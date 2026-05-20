#!/usr/bin/env python3
"""
genoio — shared genotype-file reader. Sniffs DDNA vs Illumina GSGT and
returns a uniform frame, so the fast pipelines can consume either format.

CANONICAL COPY lives at 01_mock_data_generation/scripts/genoio.py.
Synced forks (kept byte-identical, like popset.py) live in each baked
script dir so the Docker image and standalone runs both import it:
  03_individual_level/gnomad_projection_fast/scripts/genoio.py
  03_individual_level/pca_qc_fast/scripts/genoio.py
  03_individual_level/sex_biased_admixture/scripts/genoio.py
If you edit one, copy it to the others.

DDNA (Dynamic DNA): `#`-comment header, then tab rows
    rsid  chrom  pos  genotype  gs  baf  lrr
Illumina GSGT: `[Header]` / `[Data]` sections, a tab column-header row,
genotype split across `Allele1 - Plus` / `Allele2 - Plus`. bvs synthetic
GSGT carries no GenCall/BAF/LRR, so those come back NaN.

`read_genotypes(path)` returns a DataFrame with columns:
    rsid:str  chrom:str  pos:int64  gt:str(2, upper)  gs:float  baf:float  lrr:float
No filtering is applied — callers keep their existing filter logic so the
DDNA numerical path stays byte-for-byte unchanged.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

_CANON = ["rsid", "chrom", "pos", "gt", "gs", "baf", "lrr"]


def sniff_format(path) -> str:
    """Return 'illumina' or 'ddna' from the first meaningful line."""
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("[Header]") or s.startswith("[Data]"):
                return "illumina"
            return "ddna"
    return "ddna"


def _read_ddna(path) -> pd.DataFrame:
    df = pd.read_csv(
        path, sep="\t", comment="#", header=None,
        names=["rsid", "chrom", "pos", "gt", "gs", "baf", "lrr"],
        dtype={"rsid": str, "chrom": str, "pos": "Int64", "gt": str},
        na_filter=False, engine="c",
    )
    df = df[df["rsid"] != "rsid"]
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    for c in ("gs", "baf", "lrr"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["gt"] = df["gt"].astype(str).str.upper().str.strip()
    return df[_CANON].reset_index(drop=True)


def _read_illumina(path) -> pd.DataFrame:
    rows_rsid, rows_chr, rows_pos, rows_gt = [], [], [], []
    with open(path, "r") as f:
        in_data = False
        idx = None
        for line in f:
            line = line.rstrip("\n")
            if not in_data:
                if line.startswith("[Data]"):
                    in_data = True
                continue
            if idx is None:
                cols = line.split("\t")
                idx = {c: i for i, c in enumerate(cols)}
                i_name = idx.get("SNP Name")
                i_chr = idx.get("Chr")
                i_pos = idx.get("Position")
                i_a1 = idx.get("Allele1 - Plus")
                i_a2 = idx.get("Allele2 - Plus")
                if None in (i_name, i_chr, i_pos, i_a1, i_a2):
                    raise ValueError(
                        f"{path}: unexpected GSGT columns {cols[:8]}…")
                continue
            p = line.split("\t")
            if len(p) <= i_a2:
                continue
            rows_rsid.append(p[i_name])
            rows_chr.append(p[i_chr])
            rows_pos.append(p[i_pos])
            rows_gt.append((p[i_a1] + p[i_a2]).upper())
    df = pd.DataFrame({
        "rsid": pd.array(rows_rsid, dtype="string"),
        "chrom": pd.array(rows_chr, dtype="string"),
        "pos": pd.to_numeric(pd.Series(rows_pos), errors="coerce").astype("Int64"),
        "gt": pd.array(rows_gt, dtype="string"),
        "gs": np.nan,
        "baf": np.nan,
        "lrr": np.nan,
    })
    df["rsid"] = df["rsid"].astype(str)
    df["chrom"] = df["chrom"].astype(str)
    df["gt"] = df["gt"].astype(str)
    return df[_CANON].reset_index(drop=True)


def read_genotypes(path) -> pd.DataFrame:
    """Format-detecting reader. Returns the canonical 7-column frame."""
    if sniff_format(path) == "illumina":
        return _read_illumina(path)
    return _read_ddna(path)
