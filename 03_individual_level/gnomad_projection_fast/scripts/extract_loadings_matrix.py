#!/usr/bin/env python3
"""
One-time export of the gnomAD v3.1 PCA loadings Hail Table to a numpy archive
so downstream projection code can avoid the Hail/Spark JVM stack.

Output (.npz) keys:
    chrom    (n_var,)   bytes, no "chr" prefix
    pos      (n_var,)   int64
    ref      (n_var,)   bytes, single-base
    alt      (n_var,)   bytes, single-base
    loadings (n_var, n_pcs)  float64   (matches Hail's array<float64>)
    pca_af   (n_var,)   float64       (matches Hail's float64)

Storing as float64 is required for byte-precise reproduction of Hail's
pc_project output; casting to float32 introduces ~1e-6 relative errors.

Usage:
    python extract_loadings_matrix.py --ht <path/to/.ht> --out <path/to/loadings.npz>
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ht", default=os.environ.get("LOADINGS_HT"),
                    help="path to gnomad.v3.1.pca_loadings.ht (env LOADINGS_HT)")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    if not args.ht:
        sys.exit("ERROR: provide --ht or set LOADINGS_HT")

    args.out.parent.mkdir(parents=True, exist_ok=True)

    import hail as hl

    hl.init(
        default_reference="GRCh38",
        quiet=True,
        log=str(args.out.parent / "extract_loadings_matrix.log"),
        spark_conf={"spark.driver.memory": "4g"},
    )

    ht = hl.read_table(args.ht)
    ht = ht.naive_coalesce(1)
    df = ht.to_pandas()
    print(f"Loaded {len(df):,} loadings rows")

    chrom = df["locus"].apply(lambda l: l.contig.replace("chr", "")).to_numpy()
    pos = df["locus"].apply(lambda l: l.position).to_numpy(dtype=np.int64)
    ref = df["alleles"].apply(lambda a: a[0]).to_numpy()
    alt = df["alleles"].apply(lambda a: a[1]).to_numpy()
    loadings = np.array(df["loadings"].tolist(), dtype=np.float64)
    pca_af = df["pca_af"].to_numpy(dtype=np.float64)

    print(f"loadings shape: {loadings.shape}")
    print(f"pca_af range:   [{pca_af.min():.3f}, {pca_af.max():.3f}]")

    np.savez(
        args.out,
        chrom=chrom.astype("S"),
        pos=pos,
        ref=ref.astype("S"),
        alt=alt.astype("S"),
        loadings=loadings,
        pca_af=pca_af,
    )
    size_mb = args.out.stat().st_size / 1024 / 1024
    print(f"Wrote {args.out} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
