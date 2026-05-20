#!/usr/bin/env python3
"""
Assign each mock participant to one of the 6 Caribbean islands.

Algorithm: every island gets a guaranteed floor of N participants (default
100); whatever's left over is distributed by independent uniform random
draws across the islands. With 1000 participants and a floor of 100 per
island, that leaves 400 random draws — so each island ends up at roughly
167 ± 8 (binomial std), enough variation to feel random while keeping the
minimum.

Outputs <output_dir>/island_mapping.tsv with columns:
    participant_id   island
"""

import argparse
import random
from pathlib import Path

ISLANDS = ["Bahamas", "Barbados", "Bermuda", "BVI", "StLucia", "TT"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output_dir", type=Path,
                   help="dir whose numeric subdirs are the participant ids")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--min", type=int, default=100, dest="min_per_island",
                   help="minimum number of participants per island")
    args = p.parse_args()

    pids = sorted(
        d.name for d in args.output_dir.iterdir()
        if d.is_dir() and d.name.isdigit()
    )
    n = len(pids)
    n_islands = len(ISLANDS)
    floor = args.min_per_island * n_islands
    if n < floor:
        raise SystemExit(
            f"min ({args.min_per_island}) × islands ({n_islands}) = {floor} "
            f"> available participants ({n}). Lower --min or generate more samples."
        )

    rng = random.Random(args.seed)
    rng.shuffle(pids)

    assignments = []
    for i, pid in enumerate(pids[:floor]):
        assignments.append((pid, ISLANDS[i // args.min_per_island]))
    for pid in pids[floor:]:
        assignments.append((pid, rng.choice(ISLANDS)))

    assignments.sort()  # stable on pid

    out = args.output_dir / "island_mapping.tsv"
    counts = {island: 0 for island in ISLANDS}
    with out.open("w") as f:
        f.write("participant_id\tisland\n")
        for pid, island in assignments:
            f.write(f"{pid}\t{island}\n")
            counts[island] += 1

    print(f"Wrote {out}  ({n} participants)")
    for island in ISLANDS:
        print(f"  {island:10s}: {counts[island]}")


if __name__ == "__main__":
    main()
