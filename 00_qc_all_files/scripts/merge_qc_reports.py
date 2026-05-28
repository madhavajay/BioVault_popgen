#!/usr/bin/env python3
"""Merge per-file QC reports into the standard qc_all_files outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


TABLES = {
    "file_summary": [
        "participant_id",
        "file",
        "source",
        "detected_format",
        "status",
        "readable",
        "raw_rows",
        "normalized_rows",
        "unique_variants",
        "facet_count",
        "facet_missing_count",
        "missing_facets",
        "facets_json",
        "errors",
        "warnings",
        "message",
    ],
    "issues": [
        "participant_id",
        "file",
        "detected_format",
        "line_number",
        "severity",
        "code",
        "message",
        "line",
    ],
    "facet_values": ["participant_id", "file", "facet", "value", "is_missing"],
    "facet_summary": ["facet", "value", "count"],
}


def read_tsvs(paths: list[Path], columns: list[str]) -> pd.DataFrame:
    frames = []
    for path in paths:
        if path.exists() and path.stat().st_size > 0:
            frames.append(pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False))
    if not frames:
        return pd.DataFrame(columns=columns)
    df = pd.concat(frames, ignore_index=True)
    for col in columns:
        if col not in df.columns:
            df[col] = ""
    return df[columns]


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0).astype(int)


def numeric_value(value) -> int:
    parsed = pd.to_numeric(value, errors="coerce")
    if pd.isna(parsed):
        return 0
    return int(parsed)


def percent(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        return 0.0
    return round((numerator / denominator) * 100.0, 2)


def build_run_report(file_summary: pd.DataFrame, issues: pd.DataFrame) -> tuple[str, dict[str, int | float]]:
    files_checked = int(len(file_summary))
    errors = numeric(file_summary["errors"]) if "errors" in file_summary else pd.Series(dtype=int)
    raw_rows = numeric(file_summary["raw_rows"]) if "raw_rows" in file_summary else pd.Series(dtype=int)
    facet_missing = numeric(file_summary["facet_missing_count"]) if "facet_missing_count" in file_summary else pd.Series(dtype=int)
    files_with_errors = int((errors > 0).sum()) if not errors.empty else 0
    files_with_missing_facets = int((facet_missing > 0).sum()) if not facet_missing.empty else 0
    usable_files = files_checked - files_with_errors

    if not issues.empty:
        line_numbers = pd.to_numeric(issues["line_number"], errors="coerce").fillna(0)
        warning_rows = issues[(issues["severity"] == "WARNING") & (line_numbers > 0)][["file", "line_number"]].drop_duplicates()
        error_rows = issues[(issues["severity"] == "ERROR") & (line_numbers > 0)][["file", "line_number"]].drop_duplicates()
    else:
        warning_rows = pd.DataFrame(columns=["file", "line_number"])
        error_rows = pd.DataFrame(columns=["file", "line_number"])

    files_with_warning_rows = int(warning_rows["file"].nunique()) if not warning_rows.empty else 0
    total_warning_rows = int(len(warning_rows))

    raw_rows_total = int(raw_rows.sum()) if not raw_rows.empty else 0
    usable_rows = 0
    if not file_summary.empty:
        for row in file_summary.to_dict("records"):
            if numeric_value(row.get("errors", 0)) == 0:
                usable_rows += numeric_value(row.get("normalized_rows", 0))

    metrics = {
        "files_checked": files_checked,
        "files_with_errors": files_with_errors,
        "files_with_warning_rows": files_with_warning_rows,
        "total_warning_rows": total_warning_rows,
        "files_with_missing_facets": files_with_missing_facets,
        "usable_files_percent": percent(usable_files, files_checked),
        "usable_rows_percent": percent(usable_rows, raw_rows_total),
    }
    report = "\n".join([
        f"{metrics['files_checked']} number of files checked",
        f"{metrics['files_with_errors']} files with errors (unusable file)",
        f"{metrics['files_with_warning_rows']} files with warning rows",
        f"{metrics['total_warning_rows']} total warning rows",
        f"{metrics['files_with_missing_facets']} files with missing facets",
        f"{metrics['usable_files_percent']:.2f}% of usable files",
        f"{metrics['usable_rows_percent']:.2f}% of usable rows",
        "",
    ])
    return report, metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--file-summary", nargs="*", type=Path, default=[])
    parser.add_argument("--issues", nargs="*", type=Path, default=[])
    parser.add_argument("--facet-values", nargs="*", type=Path, default=[])
    parser.add_argument("--logs", nargs="*", type=Path, default=[])
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    file_summary = read_tsvs(args.file_summary, TABLES["file_summary"])
    issues = read_tsvs(args.issues, TABLES["issues"])
    facet_values = read_tsvs(args.facet_values, TABLES["facet_values"])

    if not facet_values.empty:
        grouped = (
            facet_values.assign(value=facet_values.apply(
                lambda row: "<MISSING>" if str(row["is_missing"]).lower() == "true" else row["value"],
                axis=1,
            ))
            .groupby(["facet", "value"], dropna=False)
            .size()
            .reset_index(name="count")
            .sort_values(["facet", "value"])
        )
    else:
        grouped = pd.DataFrame(columns=TABLES["facet_summary"])

    file_summary.to_csv(args.output_dir / "file_summary.tsv", sep="\t", index=False)
    issues.to_csv(args.output_dir / "issues.tsv", sep="\t", index=False)
    issues[issues["severity"] == "ERROR"].to_csv(args.output_dir / "errors.tsv", sep="\t", index=False)
    issues[issues["severity"] == "WARNING"].to_csv(args.output_dir / "warnings.tsv", sep="\t", index=False)
    issues[issues["severity"] == "QUALITY"].to_csv(args.output_dir / "quality_issues.tsv", sep="\t", index=False)
    facet_values.to_csv(args.output_dir / "facet_values.tsv", sep="\t", index=False)
    grouped.to_csv(args.output_dir / "facet_summary.tsv", sep="\t", index=False)

    log_parts = []
    for path in args.logs:
        if path.exists():
            log_parts.append(f"===== {path} =====\n{path.read_text(encoding='utf-8', errors='replace')}")
    (args.output_dir / "qc_all_files.log").write_text("\n".join(log_parts), encoding="utf-8")

    report_text, report_metrics = build_run_report(file_summary, issues)
    (args.output_dir / "report.txt").write_text(report_text, encoding="utf-8")
    totals = {
        "files": int(len(file_summary)),
        "pass": int((file_summary["status"] == "PASS").sum()) if "status" in file_summary else 0,
        "warn": int((file_summary["status"] == "WARN").sum()) if "status" in file_summary else 0,
        "fail": int((file_summary["status"] == "FAIL").sum()) if "status" in file_summary else 0,
        "errors": int(numeric(file_summary["errors"]).sum()) if "errors" in file_summary else 0,
        "warnings": int(numeric(file_summary["warnings"]).sum()) if "warnings" in file_summary else 0,
        "quality_issues": int((issues["severity"] == "QUALITY").sum()) if "severity" in issues else 0,
        "facets": sorted(facet_values["facet"].dropna().unique().tolist()) if "facet" in facet_values else [],
        "facet_missing_values": int((facet_values["is_missing"].astype(str).str.lower() == "true").sum()) if "is_missing" in facet_values else 0,
        "report": report_metrics,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(totals, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
