#!/usr/bin/env python3
"""Write HGP1K samples.tsv and variants.tsv from hgp1k_dosage.npz."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


def text(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: write_hgp1k_sidecars.py <hgp1k_dosage.npz> <dest-dir>", file=sys.stderr)
        return 2

    npz_path = Path(sys.argv[1])
    dest_dir = Path(sys.argv[2])
    z = np.load(npz_path, allow_pickle=True)

    samples_path = dest_dir / "samples.tsv"
    if not samples_path.exists() or samples_path.stat().st_size == 0:
        with samples_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample_id"])
            for sample in z["samples"]:
                writer.writerow([text(sample)])

    variants_path = dest_dir / "variants.tsv"
    if not variants_path.exists() or variants_path.stat().st_size == 0:
        with variants_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["chrom", "pos", "rsid", "ref", "alt"])
            for chrom, pos, rsid, ref, alt in zip(z["chrom"], z["pos"], z["rsid"], z["ref"], z["alt"]):
                writer.writerow([text(chrom), int(pos), text(rsid), text(ref), text(alt)])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
