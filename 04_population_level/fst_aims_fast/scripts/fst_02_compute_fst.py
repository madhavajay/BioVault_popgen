"""
FST step 2 (flow fork) - pairwise Weir & Cockerham 1984 FST.

Verbatim fork of 04_population_level/fst_islands/scripts/02_compute_fst.py.
The original is already population-agnostic (it reads whatever columns are in
merged_allele_freq.tsv), so the only change is BV_WORK_DIR so the flow can
control the working tree. See the original for the full WC84 derivation.
"""

import itertools
import logging
import os
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

BASE_DIR = Path(os.environ.get("BV_WORK_DIR", Path(__file__).resolve().parents[1]))
MERGED = BASE_DIR / "data" / "merged"
FST_DIR = BASE_DIR / "data" / "fst"
LOG_DIR = BASE_DIR / "logs"

SAVE_PER_SNP = False

LOG_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "02_compute_fst.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


def wc84_pairwise_components(
    p1: np.ndarray, n1: np.ndarray, p2: np.ndarray, n2: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = 2.0
    n_total = n1 + n2
    n_bar = n_total / r
    p_bar = (n1 * p1 + n2 * p2) / n_total
    nc = (n_total - (n1**2 + n2**2) / n_total) / (r - 1)
    s2 = (n1 * (p1 - p_bar) ** 2 + n2 * (p2 - p_bar) ** 2) / ((r - 1) * n_bar)
    h1 = np.where(n1 > 1, 2 * n1 * p1 * (1 - p1) / (2 * n1 - 1), np.nan)
    h2 = np.where(n2 > 1, 2 * n2 * p2 * (1 - p2) / (2 * n2 - 1), np.nan)
    h_bar = (h1 + h2) / r
    inner = p_bar * (1 - p_bar) - ((r - 1) / r) * s2 - h_bar / 4.0
    a = (n_bar / nc) * (s2 - inner / (n_bar - 1))
    b = (n_bar / (n_bar - 1)) * (
        p_bar * (1 - p_bar)
        - ((r - 1) / r) * s2
        - (2 * n_bar - 1) * h_bar / (4 * n_bar)
    )
    c = h_bar / 2.0
    return a, b, c


def genome_wide_fst(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    valid = np.isfinite(a) & np.isfinite(b) & np.isfinite(c)
    a_sum = a[valid].sum()
    denom = (a[valid] + b[valid] + c[valid]).sum()
    return float(a_sum / denom) if denom != 0 else float("nan")


def main() -> None:
    FST_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading merged matrices ...")
    freq = pd.read_csv(MERGED / "merged_allele_freq.tsv", sep="\t", index_col="locus_key")
    n_df = pd.read_csv(MERGED / "merged_allele_number.tsv", sep="\t", index_col="locus_key")

    populations = list(freq.columns)
    n_pops = len(populations)
    log.info(f"  {n_pops} populations x {len(freq):,} SNPs")

    freq_np = freq.values.astype(float)
    n_np = n_df.values.astype(float)

    fst_matrix = pd.DataFrame(
        np.zeros((n_pops, n_pops)), index=populations, columns=populations
    )

    pairs = list(itertools.combinations(range(n_pops), 2))
    log.info(f"Computing {len(pairs)} pairwise FST values ...")

    if SAVE_PER_SNP:
        per_snp_dir = FST_DIR / "fst_per_snp"
        per_snp_dir.mkdir(exist_ok=True)

    for i, j in pairs:
        pop_i, pop_j = populations[i], populations[j]
        p1, n1 = freq_np[:, i], n_np[:, i]
        p2, n2 = freq_np[:, j], n_np[:, j]
        a, b, c = wc84_pairwise_components(p1, n1, p2, n2)
        fst_val = genome_wide_fst(a, b, c)
        fst_matrix.loc[pop_i, pop_j] = fst_val
        fst_matrix.loc[pop_j, pop_i] = fst_val
        log.info(f"  FST({pop_i} vs {pop_j}) = {fst_val:.6f}")
        if SAVE_PER_SNP:
            per_snp = pd.DataFrame(
                {
                    "rsid": freq.index,
                    "a": a,
                    "b": b,
                    "c": c,
                    "fst_snp": np.where(a + b + c > 0, a / (a + b + c), np.nan),
                }
            )
            per_snp.to_csv(per_snp_dir / f"fst_{pop_i}_vs_{pop_j}.tsv", sep="\t", index=False)

    for p in populations:
        fst_matrix.loc[p, p] = 0.0

    out_path = FST_DIR / "fst_matrix.tsv"
    fst_matrix.to_csv(out_path, sep="\t", float_format="%.6f")
    log.info(f"\nPairwise FST matrix saved -> {out_path}")
    log.info("\n" + fst_matrix.round(5).to_string())
    log.info("\nFST step 2 complete.")


if __name__ == "__main__":
    main()
