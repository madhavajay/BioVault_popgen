#!/usr/bin/env python3
"""Preflight QC for BioVault genotype input files.

The goal is to inspect every supplied file and produce actionable diagnostics
without aborting on the first malformed sample. Each issue report includes the
file path, optional participant id, line number, issue code, and full raw line.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

for _parent in Path(__file__).resolve().parents:
    if (_parent / "tools" / "genotype_normalizer.py").exists():
        sys.path.insert(0, str(_parent))
        break
sys.path.append("/opt/biovault")
from tools import genotype_normalizer as geno  # noqa: E402


DDNA_COLS = ["rsid", "chrom", "pos", "gt", "gs", "baf", "lrr"]
VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}
NO_CALLS = {"", "-", "--", ".", "0", "00", "N", "NN", "NA", "N/A", "#N/A"}
VENDOR_MISSING_VALUES = {"#N/A"}
VALID_BASES = set("ACGT")
VALID_CHROMS_WITH_PAR = VALID_CHROMS | {"XY"}
INDEL_BASES = set("ID")


@dataclass
class InputRecord:
    participant_id: str
    path: Path
    source: str
    facets: dict[str, str]


@dataclass
class Issue:
    participant_id: str
    file: str
    detected_format: str
    line_number: int
    severity: str
    code: str
    message: str
    line: str


@dataclass
class FileSummary:
    participant_id: str
    file: str
    source: str
    detected_format: str
    status: str
    readable: bool
    raw_rows: int
    normalized_rows: int
    unique_variants: int
    facet_count: int
    facet_missing_count: int
    missing_facets: str
    facets_json: str
    errors: int
    warnings: int
    message: str


@dataclass
class FacetValue:
    facet: str
    value: str
    count: int


MISSING_FACET_VALUES = {"", "na", "n/a", "null", "none", "nan", "."}
FILE_ERROR_CODES = {
    "FILE_NOT_FOUND",
    "NOT_A_FILE",
    "FILE_NOT_READABLE",
    "FILE_EMPTY",
    "NO_DATA_ROWS",
    "READ_EXCEPTION",
    "SNIFF_EXCEPTION",
    "NORMALIZER_EXCEPTION",
    "UNEXPECTED_QC_EXCEPTION",
    "QC_PROCESS_EXIT",
    "ILLUMINA_MISSING_COLUMNS",
    "ILLUMINA_NO_HEADER",
    "ILLUMINA_NO_DATA_SECTION",
}
ROW_ERROR_WARNING_CODES = {
    "DDNA_FIELD_COUNT",
    "MALFORMED_ILLUMINA_ROW",
    "INVALID_CHROM",
    "INVALID_POSITION",
    "INVALID_GENOTYPE",
    "MISSING_GENOTYPE",
}


def report_line(line: str) -> str:
    text = line.rstrip("\n").rstrip("\r")
    return "".join(
        ch if ch == "\t" or ord(ch) >= 32 else f"\\x{ord(ch):02x}"
        for ch in text
    )


def is_missing_facet(value: str) -> bool:
    return str(value).strip().lower() in MISSING_FACET_VALUES


def facet_summary(rec: InputRecord) -> tuple[int, int, str, str]:
    missing = sorted(name for name, value in rec.facets.items() if is_missing_facet(value))
    facets_json = json.dumps(rec.facets, sort_keys=True, ensure_ascii=False)
    return len(rec.facets), len(missing), ",".join(missing), facets_json


def issue(
    rec: InputRecord,
    detected_format: str,
    line_number: int,
    severity: str,
    code: str,
    message: str,
    line: str = "",
) -> Issue:
    return Issue(
        participant_id=rec.participant_id,
        file=str(rec.path),
        detected_format=detected_format,
        line_number=line_number,
        severity=severity,
        code=code,
        message=message,
        line=report_line(line),
    )


def clean_chrom(value: str) -> str:
    chrom = str(value).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.upper()


def valid_gt(value: str) -> bool:
    gt = str(value).strip().upper()
    if gt in NO_CALLS:
        return True
    return len(gt) == 2 and all(base in VALID_BASES for base in gt)


def gt_class(value: str) -> str:
    gt = str(value).strip().upper()
    if gt in NO_CALLS:
        return "no_call"
    if len(gt) == 2 and all(base in VALID_BASES for base in gt):
        return "snp"
    if len(gt) == 2 and all(base in VALID_BASES | INDEL_BASES for base in gt) and any(base in INDEL_BASES for base in gt):
        return "indel"
    return "invalid"


def is_vendor_missing_value(value: str) -> bool:
    return str(value).strip().upper() in VENDOR_MISSING_VALUES


def setup_logging(output_dir: Path) -> logging.Logger:
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("qc_all_files")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    for handler in (
        logging.StreamHandler(),
        logging.FileHandler(output_dir / "qc_all_files.log", mode="w"),
    ):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def read_samplesheet(path: Path) -> list[InputRecord]:
    records: list[InputRecord] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError(f"{path}: samplesheet has no header")
        columns = {name.lower(): name for name in reader.fieldnames}
        pid_col = columns.get("participant_id") or columns.get("sample_id") or columns.get("id")
        file_col = columns.get("genotype_file") or columns.get("path") or columns.get("file")
        if not file_col:
            raise ValueError(f"{path}: expected genotype_file/path/file column")
        facet_cols = [name for name in reader.fieldnames if name not in {pid_col, file_col}]
        base = path.parent
        for idx, row in enumerate(reader, start=2):
            raw_file = (row.get(file_col) or "").strip()
            raw_pid = (row.get(pid_col) or "") if pid_col else ""
            pid = raw_pid.strip() or f"row_{idx}"
            fpath = Path(raw_file)
            if raw_file and not fpath.is_absolute():
                fpath = base / fpath
            facets = {name: (row.get(name) or "").strip() for name in facet_cols}
            records.append(InputRecord(pid, fpath, f"samplesheet:{path}:{idx}", facets))
    return records


def discover_inputs(paths: Iterable[Path], samplesheet: Path | None) -> list[InputRecord]:
    records: list[InputRecord] = []
    if samplesheet:
        records.extend(read_samplesheet(samplesheet))
    for path in paths:
        if path.is_file():
            records.append(InputRecord(path.stem, path, f"file:{path}", {}))
        elif path.is_dir():
            for txt in sorted(path.rglob("*.txt")):
                pid = txt.parent.name if txt.parent != path else txt.stem
                records.append(InputRecord(pid, txt, f"dir:{path}", {}))
        else:
            records.append(InputRecord(path.stem or str(path), path, f"missing-input:{path}", {}))
    records.sort(key=lambda rec: (str(rec.path), rec.participant_id))
    return records


def validate_ddna(rec: InputRecord) -> tuple[int, list[Issue]]:
    issues: list[Issue] = []
    raw_rows = 0
    seen: dict[str, tuple[str, str, str, int, str]] = {}
    try:
        with rec.path.open("r", encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                if "\x00" in line:
                    issues.append(issue(rec, "ddna", line_number, "ERROR", "NUL_BYTE", "line contains NUL byte; file may be binary or corrupt", line))
                fields = line.rstrip("\n").rstrip("\r").split("\t")
                lower = [field.strip().lstrip("\ufeff").lower() for field in fields]
                if lower[:7] == DDNA_COLS:
                    if fields[0].startswith("\ufeff"):
                        issues.append(issue(rec, "ddna", line_number, "WARNING", "UTF8_BOM", "header starts with UTF-8 BOM", line))
                    if len(fields) < 7:
                        issues.append(issue(rec, "ddna", line_number, "ERROR", "DDNA_HEADER_FIELD_COUNT", "DDNA header is missing required fields", line))
                    continue
                raw_rows += 1
                if len(fields) < 7:
                    if len(fields) == 1 and "," in fields[0]:
                        issues.append(issue(rec, "ddna", line_number, "ERROR", "POSSIBLE_CSV", "expected tab-delimited DDNA, but line looks comma-delimited", line))
                    issues.append(issue(rec, "ddna", line_number, "ERROR", "DDNA_FIELD_COUNT", f"expected at least 7 tab-separated fields, found {len(fields)}", line))
                    continue
                rsid, chrom, pos, gt, gs, baf, lrr = [field.strip() for field in fields[:7]]
                if is_vendor_missing_value(gt):
                    issues.append(issue(rec, "ddna", line_number, "QUALITY", "DDNA_VENDOR_NO_CALL_ROW", "DDNA row has vendor #N/A genotype and is skipped", line))
                    continue
                rsid_missing = rsid in {"", ".", "-"}
                chrom_clean = clean_chrom(chrom)
                if chrom_clean == "0":
                    issues.append(issue(rec, "ddna", line_number, "QUALITY", "UNMAPPED_CHROM", "chromosome is 0/unmapped; downstream SNP analyses will ignore this row", line))
                elif chrom_clean not in VALID_CHROMS_WITH_PAR:
                    issues.append(issue(rec, "ddna", line_number, "ERROR", "INVALID_CHROM", f"invalid chromosome {chrom!r}", line))
                try:
                    pos_i = int(pos)
                except ValueError:
                    issues.append(issue(rec, "ddna", line_number, "ERROR", "INVALID_POSITION", f"position is not a positive integer: {pos!r}", line))
                else:
                    if pos_i <= 0:
                        issues.append(issue(rec, "ddna", line_number, "QUALITY", "UNMAPPED_POSITION", f"position is {pos_i}; downstream SNP analyses will ignore this row", line))
                if not gt:
                    issues.append(issue(rec, "ddna", line_number, "ERROR", "MISSING_GENOTYPE", "genotype is empty; downstream code may see NaN/float", line))
                else:
                    genotype_class = gt_class(gt)
                    if genotype_class == "invalid":
                        issues.append(issue(rec, "ddna", line_number, "ERROR", "INVALID_GENOTYPE", f"genotype must be A/C/G/T SNP alleles, indel I/D alleles, or a no-call; got {gt!r}", line))
                for col_name, value in (("gs", gs), ("baf", baf), ("lrr", lrr)):
                    if value == "":
                        issues.append(issue(rec, "ddna", line_number, "WARNING", f"MISSING_{col_name.upper()}", f"{col_name} is empty", line))
                        continue
                    try:
                        numeric = float(value)
                    except ValueError:
                        issues.append(issue(rec, "ddna", line_number, "WARNING", f"INVALID_{col_name.upper()}", f"{col_name} is not numeric: {value!r}", line))
                        continue
                    if not math.isfinite(numeric):
                        if gt_class(gt) == "no_call" and col_name in {"baf", "lrr"}:
                            continue
                        issues.append(issue(rec, "ddna", line_number, "WARNING", f"NONFINITE_{col_name.upper()}", f"{col_name} is missing/non-finite: {value!r}; downstream treats this as missing", line))
                if not rsid_missing:
                    signature = (clean_chrom(chrom), pos, gt.upper(), line_number, line)
                    previous = seen.get(rsid)
                    if previous:
                        prev_chrom, prev_pos, prev_gt, prev_line_number, prev_line = previous
                        if (prev_chrom, prev_pos, prev_gt) == signature[:3]:
                            issues.append(issue(rec, "ddna", line_number, "WARNING", "DUPLICATE_RSID", f"duplicate rsid also seen at line {prev_line_number}", line))
                        else:
                            severity = "QUALITY" if gt_class(gt) == "indel" or gt_class(prev_gt) == "indel" else "ERROR"
                            issues.append(issue(rec, "ddna", line_number, severity, "CONFLICTING_DUPLICATE_RSID", f"rsid conflicts with line {prev_line_number}: {prev_line.rstrip()}", line))
                    else:
                        seen[rsid] = signature
    except Exception as exc:
        issues.append(issue(rec, "ddna", 0, "ERROR", "READ_EXCEPTION", f"{type(exc).__name__}: {exc}"))
    if rec.path.exists() and rec.path.stat().st_size == 0:
        issues.append(issue(rec, "ddna", 0, "ERROR", "FILE_EMPTY", "file is zero bytes"))
    elif raw_rows == 0 and not any(item.code.startswith("READ_") for item in issues):
        issues.append(issue(rec, "ddna", 0, "ERROR", "NO_DATA_ROWS", "file has no genotype data rows"))
    return raw_rows, issues


def validate_illumina(rec: InputRecord) -> tuple[int, list[Issue]]:
    issues: list[Issue] = []
    raw_rows = 0
    required = ["SNP Name", "Chr", "Position", "Allele1 - Plus", "Allele2 - Plus"]
    seen_probes: dict[str, int] = {}
    try:
        with rec.path.open("r", encoding="utf-8", errors="replace") as handle:
            in_data = False
            header: list[str] | None = None
            idx: dict[str, int] = {}
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not in_data:
                    if stripped == "[Data]":
                        in_data = True
                    continue
                if header is None:
                    if not stripped:
                        continue
                    header = line.rstrip("\n").rstrip("\r").split("\t")
                    idx = {name: i for i, name in enumerate(header)}
                    missing = [name for name in required if name not in idx]
                    if missing:
                        issues.append(issue(rec, "illumina", line_number, "ERROR", "ILLUMINA_MISSING_COLUMNS", f"missing required columns: {', '.join(missing)}", line))
                    continue
                if not stripped:
                    continue
                raw_rows += 1
                fields = line.rstrip("\n").rstrip("\r").split("\t")
                if len(fields) < len(header):
                    issues.append(issue(rec, "illumina", line_number, "ERROR", "ILLUMINA_FIELD_COUNT", f"expected at least {len(header)} fields from header, found {len(fields)}", line))
                    continue
                if any(name not in idx for name in required):
                    continue
                probe = fields[idx["SNP Name"]].strip()
                chrom = fields[idx["Chr"]].strip()
                pos = fields[idx["Position"]].strip()
                a1 = fields[idx["Allele1 - Plus"]].strip()
                a2 = fields[idx["Allele2 - Plus"]].strip()
                if not probe:
                    issues.append(issue(rec, "illumina", line_number, "ERROR", "MISSING_PROBE_ID", "SNP Name is empty", line))
                elif probe in seen_probes:
                    issues.append(issue(rec, "illumina", line_number, "WARNING", "DUPLICATE_PROBE", f"SNP Name also seen at line {seen_probes[probe]}", line))
                else:
                    seen_probes[probe] = line_number
                chrom_clean = clean_chrom(chrom)
                if chrom_clean == "0":
                    pass
                elif chrom_clean not in VALID_CHROMS_WITH_PAR:
                    issues.append(issue(rec, "illumina", line_number, "ERROR", "INVALID_CHROM", f"invalid chromosome {chrom!r}", line))
                try:
                    pos_i = int(pos)
                except ValueError:
                    issues.append(issue(rec, "illumina", line_number, "ERROR", "INVALID_POSITION", f"position is not a positive integer: {pos!r}", line))
                else:
                    if pos_i <= 0:
                        pass
                gt = f"{a1}{a2}"
                if not a1 or not a2:
                    issues.append(issue(rec, "illumina", line_number, "ERROR", "MISSING_ALLELE", "Allele1/Allele2 is empty", line))
                else:
                    genotype_class = gt_class(gt)
                    if genotype_class == "invalid":
                        issues.append(issue(rec, "illumina", line_number, "ERROR", "INVALID_GENOTYPE", f"alleles must be A/C/G/T, indel I/D, or no-call, got {a1!r}/{a2!r}", line))
            if not in_data:
                issues.append(issue(rec, "illumina", 0, "ERROR", "ILLUMINA_NO_DATA_SECTION", "missing [Data] section"))
            elif header is None:
                issues.append(issue(rec, "illumina", 0, "ERROR", "ILLUMINA_NO_HEADER", "missing data header after [Data]"))
            elif raw_rows == 0:
                issues.append(issue(rec, "illumina", 0, "ERROR", "NO_DATA_ROWS", "file has no genotype data rows"))
    except Exception as exc:
        issues.append(issue(rec, "illumina", 0, "ERROR", "READ_EXCEPTION", f"{type(exc).__name__}: {exc}"))
    return raw_rows, issues


def count_raw_rows_fast(rec: InputRecord, detected: str) -> int:
    raw_rows = 0
    if detected == "illumina":
        in_data = False
        saw_header = False
        with rec.path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                stripped = line.strip()
                if not in_data:
                    if stripped == "[Data]":
                        in_data = True
                    continue
                if not stripped:
                    continue
                if not saw_header:
                    saw_header = True
                    continue
                raw_rows += 1
        return raw_rows

    with rec.path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = line.rstrip("\n").rstrip("\r").split("\t")
            lower = [field.strip().lstrip("\ufeff").lower() for field in fields]
            if lower[:7] == DDNA_COLS:
                continue
            raw_rows += 1
    return raw_rows


def read_normalizer_warnings(path: Path, rec: InputRecord, detected: str) -> list[Issue]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    out: list[Issue] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            code = (row.get("code") or "NORMALIZER_WARNING").strip()
            severity = "ERROR" if code in ROW_ERROR_WARNING_CODES else "WARNING"
            try:
                line_no = int(row.get("line_no") or 0)
            except ValueError:
                line_no = 0
            out.append(issue(
                rec,
                detected,
                line_no,
                severity,
                code,
                row.get("message") or "",
                row.get("raw_line") or "",
            ))
    return out


def normalize_file(rec: InputRecord, detected: str) -> tuple[int, int, list[Issue]]:
    old_warnings = os.environ.get("BIOVAULT_WARNINGS_TSV")
    old_fast = os.environ.get("BIOVAULT_FAST_NORMALIZE")
    warning_file = tempfile.NamedTemporaryFile(prefix="biovault_qc_warnings_", suffix=".tsv", delete=False)
    warning_path = Path(warning_file.name)
    warning_file.close()
    try:
        os.environ["BIOVAULT_WARNINGS_TSV"] = str(warning_path)
        os.environ["BIOVAULT_FAST_NORMALIZE"] = "1"
        df = geno.read_genotypes(rec.path)
        rows = int(len(df))
        unique = int(df["variant_id"].nunique(dropna=True)) if rows else 0
        issues = read_normalizer_warnings(warning_path, rec, detected)
        return rows, unique, issues
    except Exception as exc:
        return 0, 0, [issue(rec, "unknown", 0, "ERROR", "NORMALIZER_EXCEPTION", f"{type(exc).__name__}: {exc}")]
    finally:
        if old_warnings is None:
            os.environ.pop("BIOVAULT_WARNINGS_TSV", None)
        else:
            os.environ["BIOVAULT_WARNINGS_TSV"] = old_warnings
        if old_fast is None:
            os.environ.pop("BIOVAULT_FAST_NORMALIZE", None)
        else:
            os.environ["BIOVAULT_FAST_NORMALIZE"] = old_fast
        try:
            warning_path.unlink()
        except FileNotFoundError:
            pass


def qc_one(rec: InputRecord, logger: logging.Logger, full_diagnostics: bool = False) -> tuple[FileSummary, list[Issue]]:
    all_issues: list[Issue] = []
    facet_count, facet_missing_count, missing_facets, facets_json = facet_summary(rec)
    logger.info(
        "FACETS participant=%s file=%s facet_count=%s missing_count=%s missing=%s values=%s",
        rec.participant_id,
        rec.path,
        facet_count,
        facet_missing_count,
        missing_facets or "-",
        facets_json,
    )
    if not rec.path.exists():
        all_issues.append(issue(rec, "unknown", 0, "ERROR", "FILE_NOT_FOUND", "file does not exist"))
        logger.info("FAIL participant=%s file=%s format=unknown raw_rows=0 normalized_rows=0 errors=1 warnings=0 message=file_not_found",
                    rec.participant_id, rec.path)
        return FileSummary(rec.participant_id, str(rec.path), rec.source, "unknown", "FAIL", False, 0, 0, 0, facet_count, facet_missing_count, missing_facets, facets_json, 1, 0, "file not found"), all_issues
    if not rec.path.is_file():
        all_issues.append(issue(rec, "unknown", 0, "ERROR", "NOT_A_FILE", "path is not a regular file"))
        logger.info("FAIL participant=%s file=%s format=unknown raw_rows=0 normalized_rows=0 errors=1 warnings=0 message=not_a_file",
                    rec.participant_id, rec.path)
        return FileSummary(rec.participant_id, str(rec.path), rec.source, "unknown", "FAIL", False, 0, 0, 0, facet_count, facet_missing_count, missing_facets, facets_json, 1, 0, "not a file"), all_issues
    if not os.access(rec.path, os.R_OK):
        all_issues.append(issue(rec, "unknown", 0, "ERROR", "FILE_NOT_READABLE", "file is not readable by this process"))
        logger.info("FAIL participant=%s file=%s format=unknown raw_rows=0 normalized_rows=0 errors=1 warnings=0 message=not_readable",
                    rec.participant_id, rec.path)
        return FileSummary(rec.participant_id, str(rec.path), rec.source, "unknown", "FAIL", False, 0, 0, 0, facet_count, facet_missing_count, missing_facets, facets_json, 1, 0, "not readable"), all_issues

    try:
        detected = geno.sniff_format(rec.path)
    except Exception as exc:
        detected = "unknown"
        all_issues.append(issue(rec, detected, 0, "ERROR", "SNIFF_EXCEPTION", f"{type(exc).__name__}: {exc}"))

    detected = detected if detected == "illumina" else ("ddna" if detected != "unknown" else detected)
    try:
        raw_rows = count_raw_rows_fast(rec, detected)
    except Exception as exc:
        raw_rows = 0
        all_issues.append(issue(rec, detected, 0, "ERROR", "READ_EXCEPTION", f"{type(exc).__name__}: {exc}"))

    normalized_rows, unique_variants, norm_issues = normalize_file(rec, detected)
    for norm_issue in norm_issues:
        norm_issue.detected_format = detected
    all_issues.extend(norm_issues)

    run_diagnostics = (
        full_diagnostics
        or normalized_rows == 0
        or (detected == "ddna" and raw_rows != normalized_rows)
        or any(item.code in FILE_ERROR_CODES for item in all_issues)
    )
    if run_diagnostics:
        if detected == "illumina":
            raw_rows, raw_issues = validate_illumina(rec)
        else:
            raw_rows, raw_issues = validate_ddna(rec)
        all_issues.extend(raw_issues)

    errors = sum(1 for item in all_issues if item.severity == "ERROR")
    warnings = sum(1 for item in all_issues if item.severity == "WARNING")
    status = "PASS" if errors == 0 else "FAIL"
    if status == "PASS" and warnings:
        status = "WARN"
    message = "ok" if status == "PASS" else f"{errors} errors, {warnings} warnings"
    logger.info("%s participant=%s file=%s format=%s raw_rows=%s normalized_rows=%s errors=%s warnings=%s",
                status, rec.participant_id, rec.path, detected, raw_rows, normalized_rows, errors, warnings)
    return FileSummary(rec.participant_id, str(rec.path), rec.source, detected, status, True, raw_rows, normalized_rows, unique_variants, facet_count, facet_missing_count, missing_facets, facets_json, errors, warnings, message), all_issues


def build_facet_counts(records: list[InputRecord]) -> list[FacetValue]:
    counts: dict[tuple[str, str], int] = {}
    for rec in records:
        for facet, raw_value in rec.facets.items():
            value = raw_value.strip()
            if is_missing_facet(value):
                value = "<MISSING>"
            counts[(facet, value)] = counts.get((facet, value), 0) + 1
    return [
        FacetValue(facet=facet, value=value, count=count)
        for (facet, value), count in sorted(counts.items(), key=lambda item: (item[0][0], item[0][1]))
    ]


def write_facet_reports(output_dir: Path, records: list[InputRecord], logger: logging.Logger) -> None:
    per_file_rows = []
    for rec in records:
        for facet, raw_value in sorted(rec.facets.items()):
            value = raw_value.strip()
            per_file_rows.append(
                {
                    "participant_id": rec.participant_id,
                    "file": str(rec.path),
                    "facet": facet,
                    "value": value,
                    "is_missing": is_missing_facet(value),
                }
            )
    pd.DataFrame(
        per_file_rows,
        columns=["participant_id", "file", "facet", "value", "is_missing"],
    ).to_csv(output_dir / "facet_values.tsv", sep="\t", index=False)

    facet_counts = build_facet_counts(records)
    pd.DataFrame([asdict(row) for row in facet_counts], columns=["facet", "value", "count"]).to_csv(
        output_dir / "facet_summary.tsv",
        sep="\t",
        index=False,
    )
    for row in facet_counts:
        logger.info("FACET_VALUE facet=%s value=%s count=%s", row.facet, row.value, row.count)


def percent(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        return 0.0
    return round((numerator / denominator) * 100.0, 2)


def split_error_frames(issue_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if issue_df.empty:
        return issue_df.copy(), issue_df.copy()
    errors = issue_df[issue_df["severity"] == "ERROR"].copy()
    line_numbers = pd.to_numeric(errors["line_number"], errors="coerce").fillna(0)
    file_mask = (line_numbers <= 0) | errors["code"].isin(FILE_ERROR_CODES)
    return errors[file_mask].copy(), errors[~file_mask].copy()


def build_run_report(summaries: list[FileSummary], issue_df: pd.DataFrame) -> tuple[str, dict[str, int | float]]:
    files_checked = len(summaries)
    file_errors, row_errors = split_error_frames(issue_df)
    files_with_file_errors = int(file_errors["file"].nunique()) if not file_errors.empty else 0
    files_with_row_errors = int(row_errors["file"].nunique()) if not row_errors.empty else 0
    total_row_error_rows = int(len(row_errors[["file", "line_number"]].drop_duplicates())) if not row_errors.empty else 0
    files_with_missing_facets = sum(1 for row in summaries if row.facet_missing_count > 0)
    usable_files = files_checked - files_with_file_errors

    warning_rows_df = issue_df[
        (issue_df["severity"] == "WARNING") & (pd.to_numeric(issue_df["line_number"], errors="coerce").fillna(0) > 0)
    ] if not issue_df.empty else pd.DataFrame(columns=issue_df.columns)
    warning_row_keys = warning_rows_df[["file", "line_number"]].drop_duplicates() if not warning_rows_df.empty else warning_rows_df
    files_with_warning_rows = int(warning_row_keys["file"].nunique()) if not warning_row_keys.empty else 0
    total_warning_rows = int(len(warning_row_keys))

    raw_rows_total = sum(int(row.raw_rows) for row in summaries)
    file_error_files = set(file_errors["file"].tolist()) if not file_errors.empty else set()
    usable_rows = sum(int(row.normalized_rows) for row in summaries if row.file not in file_error_files)

    metrics = {
        "files_checked": files_checked,
        "files_with_file_errors": files_with_file_errors,
        "files_with_row_errors": files_with_row_errors,
        "total_row_error_rows": total_row_error_rows,
        "files_with_warning_rows": files_with_warning_rows,
        "total_warning_rows": total_warning_rows,
        "files_with_missing_facets": files_with_missing_facets,
        "usable_files_percent": percent(usable_files, files_checked),
        "usable_rows_percent": percent(usable_rows, raw_rows_total),
    }
    report = "\n".join([
        f"{metrics['files_checked']} number of files checked",
        f"{metrics['files_with_file_errors']} files with file errors (unusable file)",
        f"{metrics['files_with_row_errors']} files with row errors (bad rows skipped)",
        f"{metrics['total_row_error_rows']} total row error rows",
        f"{metrics['files_with_warning_rows']} files with warning rows",
        f"{metrics['total_warning_rows']} total warning rows",
        f"{metrics['files_with_missing_facets']} files with missing facets",
        f"{metrics['usable_files_percent']:.2f}% of usable files",
        f"{metrics['usable_rows_percent']:.2f}% of usable rows",
        "",
    ])
    return report, metrics


def write_reports(output_dir: Path, records: list[InputRecord], summaries: list[FileSummary], issues: list[Issue], logger: logging.Logger) -> None:
    pd.DataFrame([asdict(row) for row in summaries]).to_csv(output_dir / "file_summary.tsv", sep="\t", index=False)
    issue_cols = ["participant_id", "file", "detected_format", "line_number", "severity", "code", "message", "line"]
    issue_rows = [asdict(row) for row in issues]
    issue_df = pd.DataFrame(issue_rows, columns=issue_cols)
    issue_df.to_csv(output_dir / "issues.tsv", sep="\t", index=False)
    issue_df[issue_df["severity"] == "ERROR"].to_csv(output_dir / "errors.tsv", sep="\t", index=False)
    file_errors, row_errors = split_error_frames(issue_df)
    file_errors.to_csv(output_dir / "file_errors.tsv", sep="\t", index=False)
    row_errors.to_csv(output_dir / "row_errors.tsv", sep="\t", index=False)
    issue_df[issue_df["severity"] == "WARNING"].to_csv(output_dir / "warnings.tsv", sep="\t", index=False)
    issue_df[issue_df["severity"] == "QUALITY"].to_csv(output_dir / "quality_issues.tsv", sep="\t", index=False)
    write_facet_reports(output_dir, records, logger)
    report_text, report_metrics = build_run_report(summaries, issue_df)
    (output_dir / "report.txt").write_text(report_text, encoding="utf-8")
    totals = {
        "files": len(summaries),
        "pass": sum(1 for row in summaries if row.status == "PASS"),
        "warn": sum(1 for row in summaries if row.status == "WARN"),
        "fail": sum(1 for row in summaries if row.status == "FAIL"),
        "errors": sum(row.errors for row in summaries),
        "warnings": sum(row.warnings for row in summaries),
        "quality_issues": int((issue_df["severity"] == "QUALITY").sum()),
        "facets": sorted({facet for rec in records for facet in rec.facets}),
        "facet_missing_values": sum(1 for rec in records for value in rec.facets.values() if is_missing_facet(value)),
        "report": report_metrics,
    }
    (output_dir / "summary.json").write_text(json.dumps(totals, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path, help="Files or directories to scan recursively for .txt genotype files")
    parser.add_argument("--samplesheet", type=Path, help="CSV/TSV containing genotype_file and optional participant_id columns")
    parser.add_argument("--output-dir", type=Path, default=Path("00_qc_all_files/results"))
    parser.add_argument("--fail-on-issues", action="store_true", help="Exit 1 when any file has an ERROR")
    parser.add_argument("--full-diagnostics", action="store_true", help="Run the slower line-by-line validator for every readable file")
    args = parser.parse_args()

    logger = setup_logging(args.output_dir)
    try:
        records = discover_inputs(args.inputs, args.samplesheet)
    except Exception as exc:
        logger.exception("Failed to discover inputs: %s", exc)
        return 2
    if not records:
        logger.error("No input files found")
        return 2

    summaries: list[FileSummary] = []
    issues: list[Issue] = []
    logger.info("Scanning %d input records", len(records))
    for rec in records:
        try:
            summary, rec_issues = qc_one(rec, logger, full_diagnostics=args.full_diagnostics)
        except Exception as exc:
            logger.exception("Unexpected QC exception for %s", rec.path)
            rec_issues = [issue(rec, "unknown", 0, "ERROR", "UNEXPECTED_QC_EXCEPTION", f"{type(exc).__name__}: {exc}")]
            facet_count, facet_missing_count, missing_facets, facets_json = facet_summary(rec)
            summary = FileSummary(rec.participant_id, str(rec.path), rec.source, "unknown", "FAIL", rec.path.exists(), 0, 0, 0, facet_count, facet_missing_count, missing_facets, facets_json, 1, 0, "unexpected QC exception")
        summaries.append(summary)
        issues.extend(rec_issues)

    write_reports(args.output_dir, records, summaries, issues, logger)
    total_errors = sum(row.errors for row in summaries)
    total_warnings = sum(row.warnings for row in summaries)
    logger.info("QC complete: files=%d pass=%d warn=%d fail=%d errors=%d warnings=%d output=%s",
                len(summaries),
                sum(1 for row in summaries if row.status == "PASS"),
                sum(1 for row in summaries if row.status == "WARN"),
                sum(1 for row in summaries if row.status == "FAIL"),
                total_errors,
                total_warnings,
                args.output_dir)
    return 1 if args.fail_on_issues and total_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
