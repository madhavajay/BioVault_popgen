#!/usr/bin/env python3
"""Freeze a reproducible 1000 Genomes reference subsample for HGP1K-anchored
ADMIXTURE (AFR / EUR / SAS).

Selection policy (deterministic):
  * source: data/1kgp_high_coverage/20130606_g1k_3202_samples_ped_population.txt
  * unrelated founders only (FatherID == 0 and MotherID == 0)
  * N per superpopulation, drawn from the founder pool sorted by SampleID
  * stdlib random.Random(seed).sample -> stable across machines/Python builds

Re-running with the same --seed and --per-pop reproduces reference_samples.tsv
byte-for-byte. The committed reference_samples.tsv is the source of truth for
the flow; this script only regenerates it.
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

SUPERPOPS = ("AFR", "EUR", "SAS")


def main() -> None:
    here = Path(__file__).resolve().parent
    repo_root = here.parents[2]
    default_ped = (
        repo_root
        / "data/1kgp_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
    )

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ped", type=Path, default=default_ped)
    ap.add_argument("--per-pop", type=int, default=300)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, default=here / "reference_samples.tsv")
    ap.add_argument(
        "--include-related",
        action="store_true",
        help="Draw from all samples, not just unrelated founders.",
    )
    args = ap.parse_args()

    pool: dict[str, list[dict[str, str]]] = {sp: [] for sp in SUPERPOPS}
    with args.ped.open() as handle:
        reader = csv.DictReader(handle, delimiter=" ")
        for row in reader:
            sp = row["Superpopulation"].strip()
            if sp not in pool:
                continue
            is_founder = row["FatherID"].strip() == "0" and row["MotherID"].strip() == "0"
            if not args.include_related and not is_founder:
                continue
            pool[sp].append(
                {
                    "sample_id": row["SampleID"].strip(),
                    "superpopulation": sp,
                    "population": row["Population"].strip(),
                    "sex": row["Sex"].strip(),
                }
            )

    selected: list[dict[str, str]] = []
    for sp in SUPERPOPS:
        candidates = sorted(pool[sp], key=lambda r: r["sample_id"])
        if len(candidates) < args.per_pop:
            raise SystemExit(
                f"ERROR: only {len(candidates)} {sp} candidates, need {args.per_pop}"
            )
        rng = random.Random(f"{args.seed}:{sp}")
        picked = rng.sample(candidates, args.per_pop)
        picked.sort(key=lambda r: r["sample_id"])
        selected.extend(picked)
        print(f"{sp}: {args.per_pop}/{len(candidates)} founders selected")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "superpopulation", "population", "sex"])
        for r in selected:
            writer.writerow(
                [r["sample_id"], r["superpopulation"], r["population"], r["sex"]]
            )
    print(f"Wrote {len(selected)} reference samples -> {args.out}")


if __name__ == "__main__":
    main()
