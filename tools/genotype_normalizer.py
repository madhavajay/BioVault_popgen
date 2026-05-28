#!/usr/bin/env python3
"""Normalize DDNA and Illumina GenomeStudio genotype files.

The parser returns one canonical row per comparable variant:

    variant_id  rsid  probe_id  chrom  pos  gt  gs  baf  lrr

`variant_id` is the stable downstream join key. For DDNA it is the rsid. For
Illumina it is the extracted/resolved rsid when available, otherwise a
`chrom:pos` locus key. Raw Illumina `SNP Name` is preserved as `probe_id`.
"""

from __future__ import annotations

import argparse
import os
import re
import sqlite3
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

CANON = ["variant_id", "rsid", "probe_id", "chrom", "pos", "gt", "gs", "baf", "lrr"]
PIPELINE_CANON = ["rsid", "chrom", "pos", "gt", "gs", "baf", "lrr"]
RS_RE = re.compile(r"rs\d+", re.IGNORECASE)
NO_CALLS = {"", "-", ".", "N", "n", "0"}


def sniff_format(path: str | Path) -> str:
    """Return 'illumina' or 'ddna' from the first meaningful line."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text in {"[Header]", "[Data]"}:
                return "illumina"
            return "ddna"
    return "ddna"


def extract_rsid(value: str) -> str:
    """Extract the first rsNNN token from an Illumina probe id."""
    match = RS_RE.search(value or "")
    if not match:
        return ""
    return "rs" + match.group(0)[2:]


def normalize_gt(a1: str, a2: str | None = None) -> str:
    if a2 is None:
        gt = str(a1).strip().upper()
        if not gt or all(base in NO_CALLS for base in gt):
            return "--"
        return gt
    a1 = str(a1).strip().upper()
    a2 = str(a2).strip().upper()
    if a1 in NO_CALLS or a2 in NO_CALLS:
        return "--"
    return f"{a1}{a2}"


def load_locus_map(path: str | Path | None) -> dict[tuple[str, int], str]:
    """Load optional chrom/pos -> rsid mappings from TSV/CSV/SQLite.

    Text files need columns compatible with `chrom`/`chromosome`, `pos`/`position`,
    and `rsid`. SQLite files are read from `grch38_non_rsids` if present, falling
    back to `rsid_reference`.
    """
    if path is None:
        return {}
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix in {".sqlite", ".sqlite3", ".db"}:
        return _load_locus_map_sqlite(path)
    sep = "," if path.suffix == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str)
    columns = {c.lower(): c for c in df.columns}
    chrom_col = columns.get("chrom") or columns.get("chromosome") or columns.get("chr")
    pos_col = columns.get("pos") or columns.get("position")
    rsid_col = columns.get("rsid") or columns.get("id")
    if not chrom_col or not pos_col or not rsid_col:
        raise ValueError(f"{path}: expected chrom/pos/rsid columns, got {list(df.columns)}")
    out: dict[tuple[str, int], str | None] = {}
    for row in df[[chrom_col, pos_col, rsid_col]].itertuples(index=False, name=None):
        chrom, pos, rsid = row
        try:
            pos_i = int(pos)
        except (TypeError, ValueError):
            continue
        rsid = str(rsid).strip()
        if chrom and pos_i > 0 and rsid:
            add_locus_mapping(out, str(chrom).strip(), pos_i, canonical_rsid(rsid))
    return {key: value for key, value in out.items() if value is not None}


def _load_locus_map_sqlite(path: Path) -> dict[tuple[str, int], str]:
    out: dict[tuple[str, int], str | None] = {}
    with sqlite3.connect(path) as con:
        tables = {
            row[0]
            for row in con.execute(
                "select name from sqlite_master where type = 'table'"
            ).fetchall()
        }
        table = "grch38_non_rsids" if "grch38_non_rsids" in tables else "rsid_reference"
        cols = {row[1].lower(): row[1] for row in con.execute(f"pragma table_info({table})")}
        chrom_col = cols.get("chromosome") or cols.get("chrom") or cols.get("chr")
        pos_col = cols.get("position") or cols.get("pos")
        rsid_col = cols.get("rsid") or cols.get("id")
        if not chrom_col or not pos_col or not rsid_col:
            raise ValueError(f"{path}:{table}: cannot find chrom/pos/rsid columns")
        query = f"select {chrom_col}, {pos_col}, {rsid_col} from {table}"
        for chrom, pos, rsid in con.execute(query):
            try:
                pos_i = int(pos)
            except (TypeError, ValueError):
                continue
            if chrom and pos_i > 0 and rsid is not None:
                add_locus_mapping(out, str(chrom), pos_i, canonical_rsid(str(rsid)))
    return {key: value for key, value in out.items() if value is not None}


def add_locus_mapping(
    out: dict[tuple[str, int], str | None], chrom: str, pos: int, rsid: str
) -> None:
    key = (chrom, pos)
    existing = out.get(key)
    if existing is None and key in out:
        return
    if existing is not None and existing != rsid:
        out[key] = None
        return
    out[key] = rsid


def canonical_rsid(value: str) -> str:
    value = str(value).strip()
    if not value:
        return ""
    return value if value.lower().startswith("rs") else f"rs{value}"


def variant_id_for(rsid: str, chrom: str, pos: int, locus_map: dict[tuple[str, int], str]) -> str:
    resolved = locus_map.get((str(chrom), int(pos)))
    if resolved:
        return resolved
    if rsid:
        return rsid
    return f"{chrom}:{int(pos)}"


def read_ddna(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=["rsid", "chrom", "pos", "gt", "gs", "baf", "lrr"],
        dtype={"rsid": str, "chrom": str, "pos": "Int64", "gt": str},
        na_filter=False,
        engine="c",
        low_memory=False,
    )
    df = df[df["rsid"] != "rsid"].copy()
    df["rsid"] = df["rsid"].map(canonical_rsid)
    df["probe_id"] = df["rsid"]
    df["variant_id"] = df["rsid"]
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    for col in ("gs", "baf", "lrr"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["gt"] = df["gt"].map(normalize_gt)
    return df[CANON].reset_index(drop=True)


def read_illumina(path: str | Path, locus_map: dict[tuple[str, int], str] | None = None) -> pd.DataFrame:
    locus_map = locus_map or {}
    rows: list[tuple[str, str, str, str, int, str, float, float, float]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        in_data = False
        idx: dict[str, int] | None = None
        for line in handle:
            line = line.rstrip("\n").rstrip("\r")
            if not in_data:
                if line.strip() == "[Data]":
                    in_data = True
                continue
            fields = line.split("\t")
            if idx is None:
                idx = {name: i for i, name in enumerate(fields)}
                required = ["SNP Name", "Chr", "Position", "Allele1 - Plus", "Allele2 - Plus"]
                missing = [name for name in required if name not in idx]
                if missing:
                    raise ValueError(f"{path}: missing Illumina columns {missing}")
                continue
            if not line.strip():
                continue
            try:
                probe_id = fields[idx["SNP Name"]].strip()
                chrom = fields[idx["Chr"]].strip()
                pos = int(fields[idx["Position"]].strip())
                a1 = fields[idx["Allele1 - Plus"]]
                a2 = fields[idx["Allele2 - Plus"]]
            except (IndexError, ValueError):
                continue
            if not chrom or chrom == "0" or pos <= 0:
                continue
            rsid = extract_rsid(probe_id)
            variant_id = variant_id_for(rsid, chrom, pos, locus_map)
            if not rsid and variant_id.startswith("rs"):
                rsid = variant_id
            rows.append((variant_id, rsid, probe_id, chrom, pos, normalize_gt(a1, a2), np.nan, np.nan, np.nan))
    df = pd.DataFrame(rows, columns=CANON)
    if df.empty:
        return df.astype({"variant_id": str, "rsid": str, "probe_id": str, "chrom": str})
    return merge_duplicate_probes(df)


def merge_duplicate_probes(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse Illumina replicate probes by variant_id/chrom/pos."""
    df = df.copy()
    keys = ["variant_id", "chrom", "pos"]
    df["_order"] = np.arange(len(df))
    df["_is_call"] = df["gt"] != "--"

    called = df[df["_is_call"]]
    if called.empty:
        conflicts = pd.MultiIndex.from_tuples([], names=keys)
    else:
        conflicts = called.groupby(keys, sort=False)["gt"].nunique()
        conflicts = conflicts[conflicts > 1].index

    chosen = (
        df.sort_values(keys + ["_is_call", "_order"], ascending=[True, True, True, False, True])
        .drop_duplicates(keys, keep="last")
        .copy()
    )
    if len(conflicts):
        chosen_idx = pd.MultiIndex.from_frame(chosen[keys])
        chosen.loc[chosen_idx.isin(conflicts), "gt"] = "--"
    out = chosen.sort_values("_order").drop(columns=["_order", "_is_call"])
    return out[CANON].reset_index(drop=True)


def read_genotypes(path: str | Path, locus_map: dict[tuple[str, int], str] | None = None) -> pd.DataFrame:
    if sniff_format(path) == "illumina":
        return read_illumina(path, locus_map=locus_map)
    return read_ddna(path)


def default_locus_map() -> dict[tuple[str, int], str]:
    env_path = os.environ.get("BIOVAULT_LOCUS_MAP")
    if env_path:
        return load_locus_map(env_path)
    default_path = Path(__file__).resolve().parent / "locus_map.tsv"
    if default_path.exists():
        return load_locus_map(default_path)
    return {}


def read_pipeline_genotypes(
    path: str | Path, locus_map: dict[tuple[str, int], str] | None = None
) -> pd.DataFrame:
    """Read DDNA/Illumina and return the historical 7-column pipeline frame.

    The `rsid` column is populated from canonical `variant_id`, so old callers
    that join on `rsid` use normalized Illumina identifiers.
    """
    if locus_map is None:
        locus_map = default_locus_map()
    df = read_genotypes(path, locus_map=locus_map)
    out = df[["variant_id", "chrom", "pos", "gt", "gs", "baf", "lrr"]].copy()
    out = out.rename(columns={"variant_id": "rsid"})
    return out[PIPELINE_CANON].reset_index(drop=True)


def summarize(df: pd.DataFrame) -> dict[str, int]:
    return {
        "rows": int(len(df)),
        "unique_variant_id": int(df["variant_id"].nunique(dropna=True)),
        "unique_rsid": int(df.loc[df["rsid"] != "", "rsid"].nunique(dropna=True)),
        "no_call": int((df["gt"] == "--").sum()),
        "locus_key": int(df["variant_id"].astype(str).str.contains(":").sum()),
    }


def compare(left: pd.DataFrame, right: pd.DataFrame) -> dict[str, int]:
    left_ids = set(left["variant_id"].astype(str))
    right_ids = set(right["variant_id"].astype(str))
    return {
        "left_rows": int(len(left)),
        "right_rows": int(len(right)),
        "left_unique": len(left_ids),
        "right_unique": len(right_ids),
        "shared_variant_id": len(left_ids & right_ids),
        "left_only": len(left_ids - right_ids),
        "right_only": len(right_ids - left_ids),
    }


def write_locus_map(paths: Iterable[str | Path], output: str | Path) -> None:
    frames = [read_genotypes(path) for path in paths]
    df = pd.concat(frames, ignore_index=True)
    df = df[df["rsid"] != ""].drop_duplicates(["chrom", "pos", "rsid"])
    df[["chrom", "pos", "rsid"]].to_csv(output, sep="\t", index=False)


def iter_illumina_probe_rows(path: str | Path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        in_data = False
        idx: dict[str, int] | None = None
        for line in handle:
            line = line.rstrip("\n").rstrip("\r")
            if not in_data:
                if line.strip() == "[Data]":
                    in_data = True
                continue
            fields = line.split("\t")
            if idx is None:
                idx = {name: i for i, name in enumerate(fields)}
                required = ["SNP Name", "SNP", "Chr", "Position"]
                missing = [name for name in required if name not in idx]
                if missing:
                    raise ValueError(f"{path}: missing Illumina columns {missing}")
                continue
            if not line.strip():
                continue
            try:
                probe_id = fields[idx["SNP Name"]].strip()
                design = fields[idx["SNP"]].strip()
                chrom = fields[idx["Chr"]].strip()
                pos = int(fields[idx["Position"]].strip())
            except (IndexError, ValueError):
                continue
            rsid = extract_rsid(probe_id)
            yield {
                "probe_id": probe_id,
                "rsid": rsid,
                "design": design,
                "chrom": chrom,
                "pos": pos,
                "has_rsid": bool(rsid),
                "plain_rsid": probe_id == rsid,
                "prefixed_or_suffixed": bool(rsid) and probe_id != rsid,
                "non_rsid": not bool(rsid),
                "duplicate_probe": "_ilmndup" in probe_id,
                "cnv_probe": "CNV" in probe_id,
                "unmapped": chrom in {"", "0"} or pos == 0,
            }


def write_probe_manifest(input_path: str | Path, output: str | Path) -> None:
    if sniff_format(input_path) != "illumina":
        raise ValueError(f"{input_path}: probe manifest requires an Illumina GSGT file")
    pd.DataFrame(iter_illumina_probe_rows(input_path)).to_csv(output, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    norm = sub.add_parser("normalize", help="Normalize one genotype file to TSV")
    norm.add_argument("input", type=Path)
    norm.add_argument("-o", "--output", type=Path)
    norm.add_argument("--locus-map", type=Path)

    comp = sub.add_parser("compare", help="Compare two genotype files after normalization")
    comp.add_argument("left", type=Path)
    comp.add_argument("right", type=Path)
    comp.add_argument("--locus-map", type=Path)

    lm = sub.add_parser("write-locus-map", help="Write chrom/pos -> rsid map from genotype files")
    lm.add_argument("inputs", nargs="+", type=Path)
    lm.add_argument("-o", "--output", type=Path, required=True)

    pm = sub.add_parser("write-probe-manifest", help="Write an Illumina probe manifest TSV")
    pm.add_argument("input", type=Path)
    pm.add_argument("-o", "--output", type=Path, required=True)

    args = parser.parse_args()
    if args.cmd == "normalize":
        locus_map = load_locus_map(args.locus_map) if args.locus_map else default_locus_map()
        df = read_genotypes(args.input, locus_map=locus_map)
        if args.output:
            df.to_csv(args.output, sep="\t", index=False)
        print(summarize(df))
    elif args.cmd == "compare":
        locus_map = load_locus_map(args.locus_map) if args.locus_map else default_locus_map()
        left = read_genotypes(args.left, locus_map=locus_map)
        right = read_genotypes(args.right, locus_map=locus_map)
        print(compare(left, right))
    elif args.cmd == "write-locus-map":
        write_locus_map(args.inputs, args.output)
        print(f"wrote {args.output}")
    elif args.cmd == "write-probe-manifest":
        write_probe_manifest(args.input, args.output)
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
