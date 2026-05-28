#!/usr/bin/env python3
"""
Generate a synthetic cohort with deliberately injected, known-truth biases
so every downstream analysis has a positive control.

Plan
----
1. Build a cohort spec: N participants, each assigned an island
   (~equal across the 6), a sex (uniform 55/45 M/F within every island so
   the sex signal can't confound island FST), and a file format
   (50/50 ddna/illumina so both readers are exercised).
2. Group participants into (island, sex, format) cells. For each cell,
   build a bvs overlay (variants-file JSON) that forces:
     - this island's structure + projection + singleton blocks  -> hom ALT
     - every other island's same blocks                          -> hom REF
     - sex_block / singleton_sex: Male -> hom ALT, Female -> hom REF
   Everything else stays at --alt-frequency 0.5 (neutral background).
3. Run `bvs synthetic` once per cell (the overlay applies to all files in
   the cell), collect the generated participant ids, and emit:
     output/<id>/<id>_*.txt
     output/cohort_spec.tsv      (participant_id, sex, island, format)
     output/island_mapping.tsv   (participant_id, island)
     output/sex_mapping.tsv      (participant_id, sex)

Usage:
    generate_biased_cohort.py --count 100 --bias-blocks bias_blocks.json \
        --out-dir 01_mock_data_generation/output [--seed 100]
"""

from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
from pathlib import Path

ISLANDS = ["Bahamas", "Barbados", "Bermuda", "BVI", "StLucia", "TT"]
SEXES = ["Male", "Female"]
FORMATS = ["ddna", "illumina"]
MALE_FRACTION = 0.55  # uniform within every island

DDNA_NAME = "{id}_X_X_GSAv3-DTC_GRCh38-{month}-{day}-{year}.txt"
ILLU_NAME = "{id}_GSAv3.txt"


def assign_spec(n: int, seed: int):
    """Return list of (island, sex, format) per participant (ids assigned
    later from bvs output)."""
    rng = random.Random(seed)
    # Even island split (remainder spread over the first islands).
    per = n // len(ISLANDS)
    rem = n % len(ISLANDS)
    plan = []
    for i, isl in enumerate(ISLANDS):
        cnt = per + (1 if i < rem else 0)
        n_male = round(cnt * MALE_FRACTION)
        sexes = ["Male"] * n_male + ["Female"] * (cnt - n_male)
        rng.shuffle(sexes)
        for j, sx in enumerate(sexes):
            fmt = FORMATS[j % 2]  # deterministic 50/50 within island
            plan.append((isl, sx, fmt))
    rng.shuffle(plan)
    return plan


def build_overlay(blocks: dict, island: str, sex: str) -> dict:
    """variants-file payload for one (island, sex) cell."""
    variants = []

    def _entry(v, gt):
        return {"rsid": v["rsid"], "chromosome": v["chromosome"],
                "position": v["position"], "reference": v["reference"],
                "alternate": v["alternate"], "genotypes": [gt]}

    def hom(v, allele):
        a = v[allele]
        return _entry(v, a + a)

    def het(v):
        return _entry(v, v["reference"] + v["alternate"])

    for isl in ISLANDS:
        tgt = isl == island
        for key in (f"island_structure[{isl}]", f"projection[{isl}]",
                    f"singleton_island[{isl}]"):
            for v in blocks.get(key, []):
                variants.append(hom(v, "alternate" if tgt else "reference"))

    # Sex control = X-hemizygosity mimicry: Male -> homozygous (low X het),
    # Female -> heterozygous (normal diploid X het). Drives panel d.
    male = sex == "Male"
    for key in ("sex_block", "singleton_sex"):
        for v in blocks.get(key, []):
            variants.append(hom(v, "reference") if male else het(v))

    return {"bias": {"variants": variants}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--count", type=int, required=True)
    ap.add_argument("--bias-blocks", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=100)
    ap.add_argument("--alt-frequency", type=float, default=0.5)
    ap.add_argument("--bvs", default="bvs")
    ap.add_argument("--from-spec", type=Path, default=None,
                    help="regenerate the exact participants in an existing "
                         "cohort_spec.tsv, preserving IDs + island/sex/format")
    ap.add_argument("--bvs-image", default=None,
                    help="run bvs via this docker image instead of host bvs "
                         "(e.g. ghcr.io/openmined/biosynth:0.1.24); the "
                         "out-dir is bind-mounted at the same absolute path")
    ap.add_argument("--biallelic", action="store_true",
                    help="pass --biallelic to bvs synthetic so every rsID is "
                         "constrained to one cohort-level allele pair")
    ap.add_argument("--preserve-filenames", action="store_true",
                    help="with --from-spec, overwrite the existing genotype "
                         "filename for each participant instead of regenerating "
                         "date placeholders")
    args = ap.parse_args()

    blocks = json.loads(args.bias_blocks.read_text())
    blocks.pop("_meta", None)

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    work = out / ".overlays"
    work.mkdir(exist_ok=True)

    def bvs_cmd(extra: list[str]) -> list[str]:
        if args.bvs_image:
            o = str(out.resolve())
            return ["docker", "run", "--rm", "-v", f"{o}:{o}",
                    "--entrypoint", "bvs", args.bvs_image] + extra
        return [args.bvs] + extra

    def synthetic_args(extra: list[str]) -> list[str]:
        cmd = ["synthetic"]
        if args.biallelic:
            cmd.append("--biallelic")
        return cmd + extra

    def existing_genotype_name(pid: str) -> str | None:
        files = sorted((out / pid).glob("*.txt"))
        if not files:
            return None
        return files[0].name

    def generated_genotype_path(pid: str, expected: Path) -> Path:
        if expected.is_file():
            return expected
        files = sorted((out / pid).glob("*.txt"))
        if len(files) == 1:
            return files[0]
        raise FileNotFoundError(
            f"could not resolve generated genotype for {pid}; expected {expected}"
        )

    # --- regenerate-from-spec: keep exact IDs + per-participant assignment ---
    if args.from_spec:
        rows = []
        for ln in args.from_spec.read_text().splitlines()[1:]:
            if not ln.strip():
                continue
            pid, sx, isl, fmt = ln.split("\t")
            rows.append((pid, sx, isl, fmt))
        print(f"Regenerating {len(rows)} participants from {args.from_spec} "
              f"(IDs preserved){' via '+args.bvs_image if args.bvs_image else ''}")
        ov_cache: dict[tuple, Path] = {}
        for i, (pid, sx, isl, fmt) in enumerate(rows):
            key = (isl, sx)
            ovf = ov_cache.get(key)
            if ovf is None:
                ovf = work / f"overlay_{isl}_{sx}.json"
                ovf.write_text(json.dumps(build_overlay(blocks, isl, sx)))
                ov_cache[key] = ovf
            name = existing_genotype_name(pid) if args.preserve_filenames else None
            if name is None:
                name = DDNA_NAME if fmt == "ddna" else ILLU_NAME
            out_path = out.resolve() / pid / name
            cmd = bvs_cmd(synthetic_args([
                "--format", fmt,
                "--output", str(out_path),
                "--count", "1",
                "--id-min", pid, "--id-max", pid,
                "--alt-frequency", str(args.alt_frequency),
                "--seed", str(args.seed + i),
                "--variants-file", str(ovf.resolve()),
            ]))
            r = subprocess.run(cmd, capture_output=True, text=True)
            try:
                generated_path = generated_genotype_path(pid, out_path)
            except FileNotFoundError:
                generated_path = out_path
            if r.returncode != 0 or not generated_path.is_file():
                sys.exit(f"bvs failed for {pid} ({isl}/{sx}/{fmt}):\n"
                         f"{r.stderr}\n{r.stdout}")
            if (i + 1) % 10 == 0 or i + 1 == len(rows):
                print(f"  {i+1}/{len(rows)}")
        (out / "cohort_spec.tsv").write_text(
            "participant_id\tsex\tisland\tformat\n"
            + "".join(f"{p}\t{s}\t{i}\t{f}\n" for p, s, i, f in
                      sorted(rows, key=lambda r: r[0])))
        (out / "island_mapping.tsv").write_text(
            "participant_id\tisland\n"
            + "".join(f"{p}\t{i}\n" for p, _s, i, _f in
                      sorted(rows, key=lambda r: r[0])))
        (out / "sex_mapping.tsv").write_text(
            "participant_id\tsex\n"
            + "".join(f"{p}\t{s}\n" for p, s, _i, _f in
                      sorted(rows, key=lambda r: r[0])))
        print(f"Regenerated {len(rows)} participants -> {out}")
        return

    plan = assign_spec(args.count, args.seed)

    # Bucket participants into cells.
    cells: dict[tuple, int] = {}
    for isl, sx, fmt in plan:
        cells[(isl, sx, fmt)] = cells.get((isl, sx, fmt), 0) + 1

    spec_rows = []  # (pid, sex, island, format)
    cell_seed = args.seed

    for (isl, sx, fmt), cnt in sorted(cells.items()):
        if cnt == 0:
            continue
        cell_seed += 1
        overlay = build_overlay(blocks, isl, sx)
        ovf = work / f"overlay_{isl}_{sx}_{fmt}.json"
        ovf.write_text(json.dumps(overlay))

        name = DDNA_NAME if fmt == "ddna" else ILLU_NAME
        out_pat = str(out / "{id}" / name)

        before = {p.name for p in out.iterdir() if p.is_dir()
                  and p.name.isdigit()}
        cmd = [args.bvs] + synthetic_args([
            "--format", fmt,
            "--output", out_pat,
            "--count", str(cnt),
            "--alt-frequency", str(args.alt_frequency),
            "--seed", str(cell_seed),
            "--variants-file", str(ovf),
        ])
        print(f"[cell] {isl:<9} {sx:<6} {fmt:<8} n={cnt}  seed={cell_seed}")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"bvs failed for {isl}/{sx}/{fmt}:\n{r.stderr}\n{r.stdout}")
        after = {p.name for p in out.iterdir() if p.is_dir()
                 and p.name.isdigit()}
        new_ids = sorted(after - before)
        if len(new_ids) != cnt:
            print(f"  WARN: expected {cnt} new dirs, got {len(new_ids)}")
        for pid in new_ids:
            spec_rows.append((pid, sx, isl, fmt))

    spec_rows.sort(key=lambda r: r[0])
    (out / "cohort_spec.tsv").write_text(
        "participant_id\tsex\tisland\tformat\n"
        + "".join(f"{p}\t{s}\t{i}\t{f}\n" for p, s, i, f in spec_rows))
    (out / "island_mapping.tsv").write_text(
        "participant_id\tisland\n"
        + "".join(f"{p}\t{i}\n" for p, _s, i, _f in spec_rows))
    (out / "sex_mapping.tsv").write_text(
        "participant_id\tsex\n"
        + "".join(f"{p}\t{s}\n" for p, s, _i, _f in spec_rows))

    n = len(spec_rows)
    from collections import Counter
    print(f"\nGenerated {n} participants")
    print("  island:", dict(Counter(r[2] for r in spec_rows)))
    print("  sex   :", dict(Counter(r[1] for r in spec_rows)))
    print("  format:", dict(Counter(r[3] for r in spec_rows)))
    print(f"  spec  : {out/'cohort_spec.tsv'}")


if __name__ == "__main__":
    main()
