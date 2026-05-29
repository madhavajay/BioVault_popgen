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
import shutil
import sys
import time
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

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

BATCH_SIZE = int(os.environ.get("BV_BATCH_SIZE", "100"))
NMF_MIN_AUTO_SIGNALS = int(os.environ.get("BV_NMF_MIN_AUTO_SIGNALS", "10000"))
NMF_MIN_X_SIGNALS = int(os.environ.get("BV_NMF_MIN_X_SIGNALS", "10000"))


def _rsid_num(rsid) -> int:
    text = str(rsid)
    if len(text) > 2 and text.startswith("rs") and text[2:].isdigit():
        return int(text[2:])
    return -1


def _finite_mean(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(arr.mean()) if arr.size else np.nan


def _compact_stats(df: pd.DataFrame, sample_id: str) -> dict:
    row = {"sample_id": sample_id}
    auto_snps = 0
    auto_hets = []
    auto_lrrs = []
    for chrom in sba.AUTOSOME_CHROMS:
        sub = df[df["chrom"] == chrom]
        auto_snps += len(sub)
        if len(sub):
            auto_hets.append(float(sub["is_het"].mean()))
            auto_lrrs.append(float(sub["lrr"].mean()))
    x = df[df["chrom"] == sba.X_CHROM]
    row.update({
        "auto_snps": auto_snps,
        "x_snps": len(x),
        "auto_het": _finite_mean(auto_hets),
        "x_het": float(x["is_het"].mean()) if len(x) else np.nan,
        "auto_lrr_mean": _finite_mean(auto_lrrs),
        "x_lrr_mean": float(x["lrr"].mean()) if len(x) else np.nan,
    })
    return row


def _signals_for_region(df: pd.DataFrame, chroms: list[str], *, autosome: bool):
    sub = df.loc[df["chrom"].isin(chroms), ["rsid", "chrom", "pos", "genotype"]].copy()
    if sub.empty:
        return (
            np.array([], dtype=np.uint32),
            np.array([], dtype=np.int8),
            np.array([], dtype=object),
            np.array([], dtype=np.float64),
        )
    sub["rsid"] = sub["rsid"].astype(str)
    sub = sub[sub["rsid"].str.match(r"^rs\d+$", na=False)]
    sub = sub.drop_duplicates("rsid")
    sub["rsid_num"] = sub["rsid"].map(_rsid_num).astype(np.int64)
    sub = sub[sub["rsid_num"] > 0]
    if autosome:
        # Match the intent of SUBSAMPLE_STEP without needing the whole cohort
        # panel in memory first.
        sub = sub[(sub["rsid_num"] % sba.SUBSAMPLE_STEP) == 0]
    if sub.empty:
        return (
            np.array([], dtype=np.uint32),
            np.array([], dtype=np.int8),
            np.array([], dtype=object),
            np.array([], dtype=np.float64),
        )
    sig = sub["genotype"].map(sba._genotype_signal)
    called = sig.notna()
    sub = sub.loc[called]
    sig = sig.loc[called]
    return (
        sub["rsid_num"].to_numpy(dtype=np.uint32),
        sig.to_numpy(dtype=np.int8),
        sub["chrom"].to_numpy(dtype=object),
        pd.to_numeric(sub["pos"], errors="coerce").to_numpy(dtype=np.float64),
    )


def _read_one_compact(args):
    sid, txt, cache_dir = args
    try:
        df = sba.read_txt(txt)
        stats = _compact_stats(df, sid)
        auto_ids, auto_sig, auto_chrom, auto_pos = _signals_for_region(
            df, sba.AUTOSOME_CHROMS, autosome=True
        )
        x_ids, x_sig, x_chrom, x_pos = _signals_for_region(
            df, [sba.X_CHROM], autosome=False
        )
        cache_path = Path(cache_dir) / f"{sid}.npz"
        np.savez(
            cache_path,
            auto_ids=auto_ids,
            auto_sig=auto_sig,
            x_ids=x_ids,
            x_sig=x_sig,
        )
        return sid, str(cache_path), stats, {
            "autosomes": (auto_ids, auto_sig, auto_chrom, auto_pos),
            "x": (x_ids, x_sig, x_chrom, x_pos),
        }, None
    except Exception as exc:
        return sid, None, None, None, str(exc).replace("\t", " ").replace("\n", " ")


def _add_region_counts(region_counts: dict, region_meta: dict, payload) -> None:
    ids, sig, chroms, poss = payload
    for vid, val, chrom, pos in zip(ids, sig, chroms, poss):
        key = int(vid)
        count, total, total_sq = region_counts.get(key, (0, 0.0, 0.0))
        v = float(val)
        region_counts[key] = (count + 1, total + v, total_sq + v * v)
        if key not in region_meta:
            region_meta[key] = (str(chrom), float(pos) if np.isfinite(pos) else np.nan)


def _discover_tasks(data_dir: Path):
    tasks = []
    errors: list[dict[str, str]] = []
    for d in sorted(p for p in data_dir.iterdir() if p.is_dir()):
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
    return tasks, errors


def _chunks(items, size: int):
    for start in range(0, len(items), size):
        yield start // size + 1, items[start:start + size]


def parse_compact_batches(data_dir: Path, workers: int, batch_size: int):
    cache_dir = sba.RESULTS_DIR / "compact_cache"
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    tasks, errors = _discover_tasks(data_dir)
    sample_cache: dict[str, Path] = {}
    stats_rows: dict[str, dict] = {}
    counts = {"autosomes": {}, "x": {}}
    meta = {"autosomes": {}, "x": {}}

    total = len(tasks)
    t0 = time.time()
    completed = 0
    ok = 0
    total_batches = (total + batch_size - 1) // batch_size if batch_size else 1
    log.info(
        "Compact parse: %d samples, batch_size=%d, parallel_batches=1, workers_per_batch=%d",
        total,
        batch_size,
        workers,
    )
    for batch_no, batch in _chunks(tasks, batch_size):
        log.info("Starting compact batch %d/%d (%d samples)", batch_no, total_batches, len(batch))
        batch_args = [(sid, txt, str(cache_dir)) for sid, txt in batch]
        if workers > 1 and len(batch_args) > 1:
            ctx = mp.get_context("fork")
            with ctx.Pool(min(workers, len(batch_args))) as pool:
                iterator = pool.imap_unordered(_read_one_compact, batch_args, chunksize=1)
                for sid, cache_path, stats, payload, err in iterator:
                    completed, ok = _handle_compact_result(
                        sid, cache_path, stats, payload, err, dict(tasks),
                        sample_cache, stats_rows, counts, meta, errors,
                        completed, ok, total, t0,
                    )
        else:
            for args in batch_args:
                sid, cache_path, stats, payload, err = _read_one_compact(args)
                completed, ok = _handle_compact_result(
                    sid, cache_path, stats, payload, err, dict(tasks),
                    sample_cache, stats_rows, counts, meta, errors,
                    completed, ok, total, t0,
                )
        log.info("Finished compact batch %d/%d", batch_no, total_batches)

    sba.write_errors(errors)
    return sample_cache, stats_rows, counts, meta


def write_compact_batch(data_dir: Path, workers: int, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_cache, stats_rows, counts, meta = parse_compact_batches(
        data_dir, workers=workers, batch_size=max(1, 10**9)
    )
    cache_out = out_dir / "cache"
    if cache_out.exists():
        shutil.rmtree(cache_out)
    cache_out.mkdir(parents=True)

    manifest_rows = []
    for sid, cache_path in sorted(sample_cache.items()):
        dst = cache_out / f"{sid}.npz"
        shutil.copy2(cache_path, dst)
        manifest_rows.append({"sample_id": sid, "cache_file": f"cache/{sid}.npz"})
    pd.DataFrame(manifest_rows, columns=["sample_id", "cache_file"]).to_csv(
        out_dir / "manifest.tsv", sep="\t", index=False
    )
    pd.DataFrame([stats_rows[sid] for sid in sorted(stats_rows)]).to_csv(
        out_dir / "stats.tsv", sep="\t", index=False
    )
    for region in ("autosomes", "x"):
        rows = []
        for vid, (called, total, total_sq) in counts[region].items():
            chrom, pos = meta[region].get(vid, ("", np.nan))
            rows.append({
                "rsid_num": int(vid),
                "chrom": chrom,
                "pos": pos,
                "called": int(called),
                "sum": float(total),
                "sum_sq": float(total_sq),
            })
        pd.DataFrame(rows, columns=["rsid_num", "chrom", "pos", "called", "sum", "sum_sq"]).to_csv(
            out_dir / f"counts_{region}.tsv", sep="\t", index=False
        )
    if sba.ERRORS_TSV.exists():
        shutil.copy2(sba.ERRORS_TSV, out_dir / "errors.tsv")
    else:
        pd.DataFrame(columns=["participant_id", "file", "severity", "code", "message"]).to_csv(
            out_dir / "errors.tsv", sep="\t", index=False
        )
    log.info("Wrote compact batch -> %s (%d usable samples)", out_dir, len(sample_cache))


def read_compact_batches(batch_dirs: list[Path]):
    sample_cache: dict[str, Path] = {}
    stats_rows: dict[str, dict] = {}
    counts = {"autosomes": {}, "x": {}}
    meta = {"autosomes": {}, "x": {}}
    error_frames = []

    for batch_dir in batch_dirs:
        manifest = pd.read_csv(batch_dir / "manifest.tsv", sep="\t", dtype=str)
        for row in manifest.itertuples(index=False):
            sample_cache[row.sample_id] = batch_dir / row.cache_file

        stats = pd.read_csv(batch_dir / "stats.tsv", sep="\t")
        for row in stats.to_dict("records"):
            stats_rows[str(row["sample_id"])] = row

        for region in ("autosomes", "x"):
            path = batch_dir / f"counts_{region}.tsv"
            if not path.exists() or path.stat().st_size == 0:
                continue
            df = pd.read_csv(path, sep="\t")
            for row in df.itertuples(index=False):
                vid = int(row.rsid_num)
                old_called, old_sum, old_sum_sq = counts[region].get(vid, (0, 0.0, 0.0))
                counts[region][vid] = (
                    old_called + int(row.called),
                    old_sum + float(row.sum),
                    old_sum_sq + float(row.sum_sq),
                )
                if vid not in meta[region]:
                    meta[region][vid] = (str(row.chrom), float(row.pos))

        err_path = batch_dir / "errors.tsv"
        if err_path.exists() and err_path.stat().st_size:
            err = pd.read_csv(err_path, sep="\t", dtype=str)
            if not err.empty:
                error_frames.append(err)

    if error_frames:
        errors = pd.concat(error_frames, ignore_index=True).drop_duplicates()
    else:
        errors = pd.DataFrame(columns=["participant_id", "file", "severity", "code", "message"])
    errors.to_csv(sba.ERRORS_TSV, sep="\t", index=False)
    return sample_cache, stats_rows, counts, meta


def _handle_compact_result(
    sid,
    cache_path,
    stats,
    payload,
    err,
    task_paths,
    sample_cache,
    stats_rows,
    counts,
    meta,
    errors,
    completed,
    ok,
    total,
    t0,
):
    completed += 1
    if err:
        path = task_paths.get(sid, "")
        log.error(f"  Skipping {sid}: {err}")
        errors.append({
            "participant_id": sid,
            "file": str(path),
            "code": "PARSE_FAILED",
            "message": err,
        })
    else:
        ok += 1
        sample_cache[sid] = Path(cache_path)
        stats_rows[sid] = stats
        _add_region_counts(counts["autosomes"], meta["autosomes"], payload["autosomes"])
        _add_region_counts(counts["x"], meta["x"], payload["x"])
    if completed % 50 == 0 or completed == total:
        elapsed = max(time.time() - t0, 1e-6)
        rate = completed / elapsed
        eta = (total - completed) / max(rate, 1e-6)
        log.info(
            "Parsed compact %d/%d samples (ok=%d skipped=%d, %.1f files/s, ETA %.0fs)",
            completed,
            total,
            ok,
            len(errors),
            rate,
            eta,
        )
    return completed, ok


def _variant_order(region: str, counts: dict, meta: dict, n_samples: int):
    rows = []
    for vid, (called, total, total_sq) in counts.items():
        call_rate = min(1.0, called / n_samples) if n_samples else 0.0
        mean = total / called if called else 0.0
        variance = (total_sq / called - mean * mean) if called else 0.0
        chrom, pos = meta.get(vid, ("", np.nan))
        rows.append((vid, chrom, pos, called, max(0, n_samples - called), call_rate, variance))
    rows.sort(key=lambda r: (sba._chrom_sort_key(r[1]), r[2] if np.isfinite(r[2]) else 1e18, r[0]))
    df = pd.DataFrame(
        rows,
        columns=["rsid_num", "chrom", "pos", "called_samples", "missing_samples", "call_rate", "variance"],
    )
    if df.empty:
        df = pd.DataFrame(columns=["rsid", "chrom", "pos", "called_samples", "missing_samples", "call_rate", "variance", "filter"])
        path = sba.RESULTS_DIR / f"nmf_variant_filter_{region}.tsv"
        df.to_csv(path, sep="\t", index=False)
        return []
    df.insert(0, "rsid", "rs" + df["rsid_num"].astype(str))
    keep = (df["call_rate"] >= sba.NMF_MIN_CALL_RATE) & (df["variance"] > 0.0)
    df["filter"] = np.where(
        df["call_rate"] < sba.NMF_MIN_CALL_RATE,
        "low_call_rate",
        np.where(df["variance"] <= 0.0, "zero_variance", "pass"),
    )
    path = sba.RESULTS_DIR / f"nmf_variant_filter_{region}.tsv"
    df.loc[~keep, ["rsid", "chrom", "pos", "called_samples", "missing_samples", "call_rate", "variance", "filter"]].to_csv(
        path, sep="\t", index=False
    )
    log.info(
        "    NMF variant filter (%s): %d/%d pass (min_call_rate=%.2f); dropped variants -> %s",
        region,
        int(keep.sum()),
        len(df),
        sba.NMF_MIN_CALL_RATE,
        path,
    )
    return df.loc[keep, "rsid_num"].astype(int).tolist()


def _nmf_sample_filter(sample_cache: dict[str, Path]) -> list[str]:
    rows = []
    eligible = []
    for sid in sorted(sample_cache):
        with np.load(sample_cache[sid]) as data:
            auto_n = int(len(data["auto_ids"]))
            x_n = int(len(data["x_ids"]))
        reason = "pass"
        if auto_n < NMF_MIN_AUTO_SIGNALS:
            reason = "low_autosome_signal"
        elif x_n < NMF_MIN_X_SIGNALS:
            reason = "low_x_signal"
        else:
            eligible.append(sid)
        rows.append({
            "sample_id": sid,
            "auto_compact_snps": auto_n,
            "x_compact_snps": x_n,
            "filter": reason,
        })

    path = sba.RESULTS_DIR / "nmf_sample_filter.tsv"
    pd.DataFrame(rows, columns=["sample_id", "auto_compact_snps", "x_compact_snps", "filter"]).to_csv(
        path, sep="\t", index=False
    )
    dropped = len(rows) - len(eligible)
    log.info(
        "    NMF sample filter: keeping %d/%d samples "
        "(min_auto_signals=%d, min_x_signals=%d); sample filter -> %s",
        len(eligible),
        len(rows),
        NMF_MIN_AUTO_SIGNALS,
        NMF_MIN_X_SIGNALS,
        path,
    )
    if dropped:
        log.warning("    NMF sample filter dropped %d low-signal samples before variant call-rate filtering", dropped)
    return eligible


def _nmf_from_cache(region: str, sample_ids: list[str], sample_cache: dict[str, Path], variant_ids: list[int]):
    log.info("    NMF: %d compact SNPs for %s", len(variant_ids), region)
    if not variant_ids:
        log.warning("    NMF skipped for %s: no variants survived call-rate and variance filters", region)
        return pd.DataFrame(np.nan, index=sample_ids, columns=["c1", "c2", "c3"])
    col = {vid: idx for idx, vid in enumerate(variant_ids)}
    mat = np.full((len(sample_ids), len(variant_ids)), 0.5, dtype=np.float32)
    ids_key = "auto_ids" if region == "autosomes" else "x_ids"
    sig_key = "auto_sig" if region == "autosomes" else "x_sig"
    for row_idx, sid in enumerate(sample_ids):
        with np.load(sample_cache[sid]) as data:
            ids = data[ids_key]
            sig = data[sig_key]
            if len(ids) == 0:
                continue
            cols = np.fromiter((col.get(int(v), -1) for v in ids), dtype=np.int64, count=len(ids))
            keep = cols >= 0
            if keep.any():
                mat[row_idx, cols[keep]] = sig[keep].astype(np.float32)
    model = sba.NMF(n_components=3, random_state=sba.NMF_SEED, max_iter=1000)
    W = model.fit_transform(mat)
    row_s = W.sum(axis=1, keepdims=True)
    row_s = np.where(row_s == 0, 1.0, row_s)
    W_norm = W / row_s
    log.info("    NMF reconstruction error: %.4f", model.reconstruction_err_)
    return pd.DataFrame(W_norm, index=sample_ids, columns=["c1", "c2", "c3"])


def build_results_table_compact(sample_cache, stats_rows, counts, meta) -> pd.DataFrame:
    sample_ids = sorted(sample_cache)
    nmf_sample_ids = _nmf_sample_filter(sample_cache)
    sex_map = sba.assign_sex(sample_ids)
    log.info("Running compact NMF on autosomes …")
    auto_variants = _variant_order("autosomes", counts["autosomes"], meta["autosomes"], len(nmf_sample_ids))
    auto_nmf = pd.DataFrame(np.nan, index=sample_ids, columns=["c1", "c2", "c3"])
    auto_nmf.update(_nmf_from_cache("autosomes", nmf_sample_ids, sample_cache, auto_variants))
    log.info("Running compact NMF on X chromosome …")
    x_variants = _variant_order("x", counts["x"], meta["x"], len(nmf_sample_ids))
    x_nmf = pd.DataFrame(np.nan, index=sample_ids, columns=["c1", "c2", "c3"])
    x_nmf.update(_nmf_from_cache("x", nmf_sample_ids, sample_cache, x_variants))

    rows = []
    for sid in sample_ids:
        stats = stats_rows[sid]
        rows.append({
            "sample_id": sid,
            "sex": sex_map[sid],
            "auto_snps": stats["auto_snps"],
            "x_snps": stats["x_snps"],
            "auto_het": round(stats["auto_het"], 4),
            "x_het": round(stats["x_het"], 4),
            "auto_lrr_mean": round(stats["auto_lrr_mean"], 4),
            "x_lrr_mean": round(stats["x_lrr_mean"], 4),
            "auto_c1": round(auto_nmf.loc[sid, "c1"], 4),
            "auto_c2": round(auto_nmf.loc[sid, "c2"], 4),
            "auto_c3": round(auto_nmf.loc[sid, "c3"], 4),
            "x_c1": round(x_nmf.loc[sid, "c1"], 4),
            "x_c2": round(x_nmf.loc[sid, "c2"], 4),
            "x_c3": round(x_nmf.loc[sid, "c3"], 4),
            "x_minus_auto_c1": round(x_nmf.loc[sid, "c1"] - auto_nmf.loc[sid, "c1"], 4),
        })
    return pd.DataFrame(rows)


def _read_one(args):
    sid, txt = args
    try:
        return sid, sba.read_txt(txt), None
    except Exception as exc:
        return sid, None, str(exc).replace("\t", " ").replace("\n", " ")


def load_all_samples_parallel(data_dir: Path, workers: int) -> dict:
    dirs = sorted(d for d in data_dir.iterdir() if d.is_dir())
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
    t0 = time.time()
    completed = 0
    if workers > 1 and len(tasks) > 1:
        ctx = mp.get_context("fork")
        with ctx.Pool(min(workers, len(tasks))) as pool:
            for sid, df, err in pool.imap_unordered(_read_one, tasks, chunksize=1):
                completed += 1
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
                if completed % 50 == 0 or completed == len(tasks):
                    elapsed = max(time.time() - t0, 1e-6)
                    rate = completed / elapsed
                    eta = (len(tasks) - completed) / max(rate, 1e-6)
                    log.info(
                        "Loaded %d/%d samples (ok=%d skipped=%d, %.1f files/s, ETA %.0fs)",
                        completed,
                        len(tasks),
                        len(out),
                        len(errors),
                        rate,
                        eta,
                    )
    else:
        for t in tasks:
            sid, df, err = _read_one(t)
            completed += 1
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
            if completed % 50 == 0 or completed == len(tasks):
                elapsed = max(time.time() - t0, 1e-6)
                rate = completed / elapsed
                eta = (len(tasks) - completed) / max(rate, 1e-6)
                log.info(
                    "Loaded %d/%d samples (ok=%d skipped=%d, %.1f files/s, ETA %.0fs)",
                    completed,
                    len(tasks),
                    len(out),
                    len(errors),
                    rate,
                    eta,
                )
    sba.write_errors(errors)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("run", "parse-batch", "finalize"),
        default=os.environ.get("BV_MODE", "run"),
    )
    parser.add_argument("--batch-out", type=Path, default=Path("compact_batch"))
    parser.add_argument("batch_dirs", nargs="*", type=Path)
    args = parser.parse_args()

    workers = int(os.environ.get("BV_WORKERS",
                                 "8"))
    workers = max(1, workers)
    log.info("=== Sex-biased admixture (fast) ===")
    log.info(
        "Mode: %s | Data dir: %s | batch_size=%d | workers_per_batch=%d | parallel_batches=1",
        args.mode,
        sba.DATA_DIR,
        BATCH_SIZE,
        workers,
    )

    t0 = time.time()
    if args.mode == "parse-batch":
        if not sba.DATA_DIR.exists():
            log.error(f"Data directory not found: {sba.DATA_DIR}")
            sys.exit(1)
        write_compact_batch(sba.DATA_DIR, workers=workers, out_dir=args.batch_out)
        return

    if args.mode == "finalize":
        if not args.batch_dirs:
            log.error("No compact batch directories provided for finalize")
            sys.exit(1)
        sample_cache, stats_rows, counts, meta = read_compact_batches(args.batch_dirs)
    else:
        if not sba.DATA_DIR.exists():
            log.error(f"Data directory not found: {sba.DATA_DIR}")
            sys.exit(1)
        sample_cache, stats_rows, counts, meta = parse_compact_batches(
            sba.DATA_DIR, workers=workers, batch_size=BATCH_SIZE
        )

    if not sample_cache:
        log.error("No samples found.")
        sys.exit(1)
    log.info(f"Loaded compact cache for {len(sample_cache)} samples in {time.time()-t0:.1f}s")

    df = build_results_table_compact(sample_cache, stats_rows, counts, meta)

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
