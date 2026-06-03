#!/usr/bin/env python3
"""Bump BioVault_popgen image and flow versions in known repo files."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SEMVER = re.compile(r"^\d+\.\d+\.\d+$")


def replace(path: Path, patterns: list[tuple[str, str]], dry_run: bool) -> bool:
    text = path.read_text(encoding="utf-8")
    updated = text
    for pattern, repl in patterns:
        updated = re.sub(pattern, repl, updated)
    if updated == text:
        return False
    if not dry_run:
        path.write_text(updated, encoding="utf-8")
    print(("would update " if dry_run else "updated ") + str(path.relative_to(ROOT)))
    return True


def require_semver(value: str, name: str) -> None:
    if not SEMVER.match(value):
        raise SystemExit(f"{name} must be x.y.z, got {value!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", help="BioVault_popgen version, for example 0.1.4")
    parser.add_argument(
        "--biosynth-version",
        help="Optional BioSynth version to update in flow 04 and image_versions.sh",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print files that would change")
    args = parser.parse_args()

    require_semver(args.version, "version")
    if args.biosynth_version:
        require_semver(args.biosynth_version, "biosynth-version")

    v = args.version
    bio_patterns = [
        (r"VERSION=\d+\.\d+\.\d+ \./build_docker\.sh", f"VERSION={v} ./build_docker.sh"),
        (r'VERSION="\$\{VERSION:-\d+\.\d+\.\d+\}"', f'VERSION="${{VERSION:-{v}}}"'),
        (r"ghcr\.io/madhavajay/biovault-popgen:\d+\.\d+\.\d+-fast", f"ghcr.io/madhavajay/biovault-popgen:{v}-fast"),
        (r"ghcr\.io/madhavajay/biovault-popgen:\d+\.\d+\.\d+(?!-fast)", f"ghcr.io/madhavajay/biovault-popgen:{v}"),
        (r"version: \d+\.\d+\.\d+", f"version: {v}"),
    ]

    files: set[Path] = {
        ROOT / "build_docker.sh",
        ROOT / "00_qc_all_files.sh",
        ROOT / "BIOVAULT.md",
        ROOT / "TODO.md",
        ROOT / "scripts" / "image_versions.sh",
    }
    files.update((ROOT / "flows").glob("*/flow.yaml"))
    files.update((ROOT / "flows").glob("*/module.yaml"))
    files.update((ROOT / "flows").glob("*/main.nf"))

    changed = 0
    for path in sorted(files):
        if path.exists() and replace(path, bio_patterns, args.dry_run):
            changed += 1

    if args.biosynth_version:
        bv = args.biosynth_version
        biosynth_patterns = [
            (r"ghcr\.io/openmined/biosynth:0\.1\.\d+", f"ghcr.io/openmined/biosynth:{bv}"),
            (r"`biosynth:0\.1\.\d+`", f"`biosynth:{bv}`"),
            (r"biosynth `--alt-frequency 0\.5`", "biosynth `--alt-frequency 0.5`"),
            (r"biosynth 0\.1\.\d+", f"biosynth {bv}"),
        ]
        for path in (
            ROOT / "BIOVAULT.md",
            ROOT / "TODO.md",
            ROOT / "scripts" / "image_versions.sh",
            ROOT / "flows" / "04_bv_paper_population_level" / "main.nf",
        ):
            if replace(path, biosynth_patterns, args.dry_run):
                changed += 1

    if changed == 0:
        print("no files changed")


if __name__ == "__main__":
    main()
