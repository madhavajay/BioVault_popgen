#!/usr/bin/env python3
"""
sex_biased_admixture_fast — byte-identical fast sibling of
../sex_biased_admixture.

The only bottleneck in the original is `load_all_samples` reading ~100
genotype files serially with pandas. Everything downstream
(`assign_sex`, `chrom_stats_for_sample`, `nmf_components` with fixed
`random_state=42`, `build_results_table`, `plot_figure`) is deterministic.

So this script imports the original module verbatim and *only* replaces
the loader with a parallel one (multiprocessing fork pool, same
`read_txt` per file). Identical per-sample frames in → identical
`sex_bias_results.tsv` out (NMF is seeded; all downstream uses
`sorted(samples)`, so parallel collection order is irrelevant).

`results/sex_bias_results.tsv` is guaranteed byte-for-byte identical.
The figure is rendered from the same DataFrame by the same code; pixel-
identical, though matplotlib embeds a timestamp so the PDF/PNG bytes
match only when SOURCE_DATE_EPOCH is pinned (the runner sets it).

Outputs (in sex_biased_admixture_fast/):
  results/sex_bias_results.tsv
  plots/figure4_sex_biased_admixture.{pdf,png}
  logs/sex_biased_admixture.log
"""

from __future__ import annotations

import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

ORIG_SCRIPTS = (Path(__file__).resolve().parents[2]
                / "sex_biased_admixture" / "scripts")
sys.path.insert(0, str(ORIG_SCRIPTS))

import sex_biased_admixture as sba  # noqa: E402  reuse ALL analysis math

# Redirect outputs into the _fast tree (original module hardcodes its own).
_BASE = Path(__file__).resolve().parents[1]
sba.RESULTS_DIR = _BASE / "results"
sba.PLOTS_DIR = _BASE / "plots"
sba.LOGS_DIR = _BASE / "logs"
sba.ERRORS_TSV = sba.LOGS_DIR / "errors.tsv"
sba.WARNINGS_TSV = sba.LOGS_DIR / "warnings.tsv"
for _d in (sba.RESULTS_DIR, sba.PLOTS_DIR, sba.LOGS_DIR):
    _d.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("BIOVAULT_FAST_NORMALIZE", "1")
os.environ.setdefault("BIOVAULT_WARNINGS_TSV", str(sba.WARNINGS_TSV))

# Honor BIOVAULT_DATA_DIR so the step runner / flow can point at a sample
# subset (mirrors pca_qc_fast). assign_sex already reads BIOVAULT_SEX_MAPPING
# or <data_dir>/sex_mapping.tsv (the `sex` facet, same as `country`).
_dd = os.environ.get("BIOVAULT_DATA_DIR")
if _dd:
    sba.DATA_DIR = Path(_dd).resolve()

log = sba.log


def _read_one(args):
    sid, txt = args
    try:
        return sid, sba.read_txt(txt), None
    except Exception as exc:
        return sid, None, str(exc).replace("\t", " ").replace("\n", " ")


def load_all_samples_parallel(data_dir: Path, workers: int) -> dict:
    dirs = sorted(d for d in data_dir.iterdir()
                  if d.is_dir() and d.name.isdigit())
    tasks = []
    errors: list[dict[str, str]] = []
    for d in dirs:
        txts = list(d.glob("*.txt"))
        if not txts:
            log.warning(f"  No .txt file in {d}; skipping")
            errors.append({
                "participant_id": d.name,
                "file": str(d),
                "code": "NO_TXT_FILE",
                "message": "sample directory contains no .txt genotype file",
            })
            continue
        tasks.append((d.name, txts[0]))
    out: dict = {}
    if workers > 1 and len(tasks) > 1:
        ctx = mp.get_context("fork")
        with ctx.Pool(min(workers, len(tasks))) as pool:
            for sid, df, err in pool.imap_unordered(_read_one, tasks, chunksize=1):
                if err:
                    path = dict(tasks).get(sid, "")
                    log.error(f"  Skipping {sid}: {err}")
                    errors.append({
                        "participant_id": sid,
                        "file": str(path),
                        "code": "PARSE_FAILED",
                        "message": err,
                    })
                else:
                    out[sid] = df
    else:
        for t in tasks:
            sid, df, err = _read_one(t)
            if err:
                log.error(f"  Skipping {sid}: {err}")
                errors.append({
                    "participant_id": sid,
                    "file": str(t[1]),
                    "code": "PARSE_FAILED",
                    "message": err,
                })
            else:
                out[sid] = df
    sba.write_errors(errors)
    return out


def main():
    workers = int(os.environ.get("BV_WORKERS",
                                 max(1, (os.cpu_count() or 4))))
    log.info("=== Sex-biased admixture (fast) ===")
    log.info(f"Data dir : {sba.DATA_DIR}  | workers={workers}")
    if not sba.DATA_DIR.exists():
        log.error(f"Data directory not found: {sba.DATA_DIR}")
        sys.exit(1)

    t0 = time.time()
    samples = load_all_samples_parallel(sba.DATA_DIR, workers)
    if not samples:
        log.error("No samples found.")
        sys.exit(1)
    log.info(f"Loaded {len(samples)} samples in {time.time()-t0:.1f}s "
             f"(parallel)")

    df = sba.build_results_table(samples)

    out_tsv = sba.RESULTS_DIR / "sex_bias_results.tsv"
    df.to_csv(out_tsv, sep="\t", index=False, float_format="%.4f")
    log.info(f"Results → {out_tsv}")
    log.info("\n" + df[[
        "sample_id", "sex", "auto_het", "x_het",
        "auto_c1", "x_c1", "x_minus_auto_c1",
    ]].to_string(index=False))

    sba.plot_figure(df)
    log.info("Done.")


if __name__ == "__main__":
    main()
