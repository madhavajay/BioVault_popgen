#!/usr/bin/env python3
"""
Numpy replacement for hl.experimental.pc_project. Reads study PLINK
bed/bim/fam directly and projects onto the gnomAD HGDP+1kGP PCA space using
the pre-exported loadings.npz.

This is a drop-in for `pca_project.py` (same CLI, same outputs) that skips
the ~45 s Hail/Spark JVM warm-up. The numerical contract is bit-for-bit
match against Hail's pc_project:

  1. Variants are joined on (locus, alleles) with a2_reference=True, so the
     join keeps a variant iff study's (A2, A1) == loadings (ref, alt).
  2. n_variants is the TOTAL row count of the loadings table (the source
     passed to pc_project). It is NOT the post-join count — Hail evaluates
     `loadings_source.count()` once up-front and uses that constant.
  3. dosage = count of A1 = count of alt (because A1 == loadings.alt in
     surviving rows).
  4. Per-variant standardization: (G - 2*af) / sqrt(n_variants * 2*af*(1-af))
  5. Missing genotypes contribute 0 (Hail's hl.agg.array_sum skips
     undefined entries).

Usage:
    fast_pca_project.py <study_plink_prefix> <output_dir>

Outputs (in <output_dir>):
    study_pca_projection.tsv      one row per sample: s<TAB>[scores...]
    pca_projection.png            PC1 vs PC2 scatter
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np


LOADINGS_NPZ_DEFAULT = "/opt/biovault/reference/pca_loadings/loadings.npz"


def read_bim(prefix: str):
    chrom, rsid, pos, a1, a2 = [], [], [], [], []
    with open(f"{prefix}.bim") as f:
        for line in f:
            parts = line.rstrip("\n").split()
            chrom.append(parts[0])
            rsid.append(parts[1])
            pos.append(int(parts[3]))
            a1.append(parts[4])
            a2.append(parts[5])
    return (
        np.array(chrom, dtype=object),
        np.array(rsid, dtype=object),
        np.array(pos, dtype=np.int64),
        np.array(a1, dtype=object),
        np.array(a2, dtype=object),
    )


def read_fam(prefix: str):
    sample_ids = []
    with open(f"{prefix}.fam") as f:
        for line in f:
            parts = line.rstrip("\n").split()
            sample_ids.append(parts[1])
    return sample_ids


def read_bed_dosage(prefix: str, n_variants: int, n_samples: int) -> np.ndarray:
    """Return dosage matrix shape (n_variants, n_samples), float32, NaN for missing.

    PLINK 1 bed (variant-major) packs 4 sample genotypes per byte, bits 0-1
    sample0 ... bits 6-7 sample3. Codes: 0b00 hom A1, 0b01 missing,
    0b10 het, 0b11 hom A2. With a2_reference=True, n_alt = count of A1, so
    dosage = 2,1,0 for 0b00,0b10,0b11; NaN for 0b01.
    """
    bytes_per_var = (n_samples + 3) // 4
    expected = n_variants * bytes_per_var
    with open(f"{prefix}.bed", "rb") as f:
        magic = f.read(3)
        if magic != b"\x6c\x1b\x01":
            raise ValueError(f"bad bed magic: {magic!r}")
        raw = f.read()
    if len(raw) != expected:
        raise ValueError(f"bed size mismatch: got {len(raw)}, expected {expected}")

    arr = np.frombuffer(raw, dtype=np.uint8).reshape(n_variants, bytes_per_var)
    codes = np.empty((n_variants, bytes_per_var * 4), dtype=np.uint8)
    codes[:, 0::4] = arr & 0b11
    codes[:, 1::4] = (arr >> 2) & 0b11
    codes[:, 2::4] = (arr >> 4) & 0b11
    codes[:, 3::4] = (arr >> 6) & 0b11
    codes = codes[:, :n_samples]

    dosage = np.full(codes.shape, np.nan, dtype=np.float32)
    dosage[codes == 0b00] = 2.0
    dosage[codes == 0b10] = 1.0
    dosage[codes == 0b11] = 0.0
    return dosage


def load_loadings(path: str):
    arc = np.load(path)
    chrom = np.array([s.decode() if isinstance(s, bytes) else str(s)
                      for s in arc["chrom"]], dtype=object)
    pos = arc["pos"].astype(np.int64)
    ref = np.array([s.decode() if isinstance(s, bytes) else str(s)
                    for s in arc["ref"]], dtype=object)
    alt = np.array([s.decode() if isinstance(s, bytes) else str(s)
                    for s in arc["alt"]], dtype=object)
    loadings = arc["loadings"].astype(np.float64)
    pca_af = arc["pca_af"].astype(np.float64)
    return chrom, pos, ref, alt, loadings, pca_af


def java_double_str(x: float) -> str:
    """Format a float the way Java's Double.toString does, to match Hail's
    Table.export output byte-for-byte.

    Rules (Java spec):
      - 0.0 -> "0.0"  (negative zero -> "-0.0")
      - finite, 1e-3 <= |x| < 1e7: decimal notation, minimum 1 digit after dot
      - otherwise: scientific d[.ddd]E[-]dd  with one digit before the dot
      - shortest decimal representation that round-trips to the same double
    """
    if x == 0.0:
        return "-0.0" if (1.0 / x) < 0 else "0.0"
    if not np.isfinite(x):
        return "NaN" if np.isnan(x) else ("Infinity" if x > 0 else "-Infinity")
    sign = "-" if x < 0 else ""
    ax = abs(x)
    # Python's repr() yields the shortest round-trippable decimal for a double.
    r = repr(ax)
    if "e" in r or "E" in r:
        mant, exp = r.lower().split("e")
        e = int(exp)
    else:
        # decimal form; convert to mantissa + exponent
        if "." in r:
            intpart, frac = r.split(".")
        else:
            intpart, frac = r, ""
        if intpart == "0":
            # 0.000XYZ -> shift to single leading digit
            stripped = frac.lstrip("0")
            shift = len(frac) - len(stripped)
            digits = stripped if stripped else "0"
            e = -(shift + 1)
            mant = digits[0] + ("." + digits[1:] if len(digits) > 1 else "")
        else:
            e = len(intpart) - 1
            digits = intpart + frac
            mant = digits[0] + ("." + digits[1:] if len(digits) > 1 else "")
    # Strip trailing zeros / dot in mantissa
    if "." in mant:
        mant = mant.rstrip("0").rstrip(".")
    # Apply Java's threshold rule using exponent and original magnitude
    if -3 <= e < 7:
        # decimal output
        if "." not in mant:
            mant_digits = mant
            mant_int_len = len(mant_digits)
        else:
            int_p, frac_p = mant.split(".")
            mant_digits = int_p + frac_p
            mant_int_len = len(int_p)
        # We have number = mant_digits * 10^(e - (mant_int_len - 1))
        # Shift to standard decimal:
        total_int_digits = e + 1  # digits before decimal
        if total_int_digits <= 0:
            # leading zeros: 0.<zeros><digits>
            s = "0." + "0" * (-total_int_digits) + mant_digits
        elif total_int_digits >= len(mant_digits):
            s = mant_digits + "0" * (total_int_digits - len(mant_digits)) + ".0"
        else:
            s = mant_digits[:total_int_digits] + "." + mant_digits[total_int_digits:]
        return sign + s
    else:
        # scientific
        if "." not in mant:
            mant = mant + ".0" if mant.isdigit() and len(mant) == 1 else mant
        # Java omits trailing zeros but always has a digit after the dot when present.
        # If mantissa has no dot (single digit), Java writes e.g. "1.0E10" not "1E10".
        if "." not in mant:
            mant = mant + ".0"
        return f"{sign}{mant}E{e}"


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    prefix = sys.argv[1]
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    loadings_path = os.environ.get("LOADINGS_NPZ", LOADINGS_NPZ_DEFAULT)
    if not os.path.exists(loadings_path):
        sys.exit(f"ERROR: loadings.npz not found at {loadings_path}; "
                 "run extract_loadings_matrix.py first.")

    t0 = time.time()
    chrom, rsid, pos, a1, a2 = read_bim(prefix)
    sample_ids = read_fam(prefix)
    n_var_study = len(chrom)
    n_samples = len(sample_ids)
    print(f"[fast_pca_project] study: {n_var_study:,} variants x {n_samples} samples")

    dosage = read_bed_dosage(prefix, n_var_study, n_samples)
    t_read = time.time() - t0
    print(f"[fast_pca_project] read bed/bim/fam in {t_read:.2f}s")

    t0 = time.time()
    L_chrom, L_pos, L_ref, L_alt, L_load, L_af = load_loadings(loadings_path)
    print(f"[fast_pca_project] loadings: {len(L_pos):,} variants x {L_load.shape[1]} PCs")

    # Index loadings by (chrom, pos, ref, alt). Study chrom may have "chr"
    # prefix in some import paths; strip for both sides to be safe.
    def norm_chrom(arr):
        return np.array([c[3:] if isinstance(c, str) and c.startswith("chr") else c
                         for c in arr], dtype=object)
    s_chrom = norm_chrom(chrom)
    l_chrom = norm_chrom(L_chrom)
    l_index = {(c, int(p), r, a): i for i, (c, p, r, a) in
               enumerate(zip(l_chrom, L_pos, L_ref, L_alt))}

    # Hail import_plink with a2_reference=True puts alleles=[A2, A1] in the
    # row key. Loadings key is [ref, alt]. Join succeeds iff A2==ref AND
    # A1==alt.
    keep_study, keep_load = [], []
    for i in range(n_var_study):
        idx = l_index.get((s_chrom[i], int(pos[i]), a2[i], a1[i]))
        if idx is None:
            continue
        af = L_af[idx]
        if not (0.0 < af < 1.0):
            continue
        keep_study.append(i)
        keep_load.append(idx)
    keep_study = np.asarray(keep_study, dtype=np.int64)
    keep_load = np.asarray(keep_load, dtype=np.int64)
    n_joined = int(keep_study.size)
    n_var_loadings_total = int(len(L_pos))  # Hail uses loadings.count(), not the join size
    t_join = time.time() - t0
    print(f"[fast_pca_project] joined {n_joined:,} of {n_var_loadings_total:,} "
          f"loadings variants in {t_join:.2f}s")

    if n_joined == 0:
        sys.exit("ERROR: zero variants survived the join — check alleles.")

    t0 = time.time()
    G = dosage[keep_study, :].astype(np.float64)         # (n_joined, n_samples)
    af = L_af[keep_load]                                  # (n_joined,)
    L = L_load[keep_load]                                 # (n_joined, n_pcs)

    # Hail's pc_project normalizes by the TOTAL loadings rowcount, not the
    # post-join count. See hail/experimental/pca.py: n_variants = loadings_source.count().
    sd = np.sqrt(n_var_loadings_total * 2.0 * af * (1.0 - af))  # (n_joined,)
    L_scaled = L / sd[:, None]

    dev = G - 2.0 * af[:, None]
    np.copyto(dev, 0.0, where=~np.isfinite(dev))
    scores = dev.T @ L_scaled                             # (n_samples, n_pcs)
    t_proj = time.time() - t0
    print(f"[fast_pca_project] projection in {t_proj:.3f}s")

    out_tsv = out_dir / "study_pca_projection.tsv"
    with out_tsv.open("w") as f:
        f.write("s\tscores\n")
        for sid, row in zip(sample_ids, scores):
            f.write(f"{sid}\t[{','.join(java_double_str(v) for v in row.tolist())}]\n")
    print(f"[fast_pca_project] wrote {out_tsv}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(scores[:, 0], scores[:, 1], c="black", marker="*", s=200,
                   edgecolors="gold", linewidth=0.8)
        for i, sid in enumerate(sample_ids):
            ax.annotate(str(sid), (scores[i, 0], scores[i, 1]),
                        fontsize=8, alpha=0.7, xytext=(4, 4),
                        textcoords="offset points")
        ax.set_xlabel("PC1 (gnomAD HGDP+1kGP space)")
        ax.set_ylabel("PC2 (gnomAD HGDP+1kGP space)")
        ax.set_title("Study samples projected onto gnomAD reference PCs")
        ax.grid(True, alpha=0.2)
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.axvline(0, color="grey", linewidth=0.5)
        plot_path = out_dir / "pca_projection.png"
        plt.tight_layout()
        plt.savefig(plot_path, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"[fast_pca_project] plot -> {plot_path}")
    except Exception as e:
        print(f"[fast_pca_project] plot skipped: {e}")


if __name__ == "__main__":
    main()
