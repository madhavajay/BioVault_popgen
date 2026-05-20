"""
Shared population resolver for the bv_paper_fst_island_aims flow.

The flow is driven by a `country` participant facet. main.nf normalizes each
distinct facet value (trim -> lowercase -> non-alphanumeric to "_" ->
collapse/strip "_") and exports the resulting set as BV_POPULATIONS
(comma-separated). Every forked script resolves its population list from here
so the four formerly-island-hardcoded steps stay in lockstep with whatever
countries the cohort actually has.

Fail-loud contract: if BV_POPULATIONS is set, the listed populations are
authoritative. A step that needs the per-country allele-freq files asserts
each `allele_freq_<pop>.tsv` exists and is non-empty and raises listing every
missing/empty one. If BV_POPULATIONS is unset (running the source pipeline by
hand) we fall back to discovering `allele_freq_*.tsv` in the raw dir.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

_NORM_RE = re.compile(r"[^a-z0-9]+")


def normalize(value: str) -> str:
    """trim -> lowercase -> non-alphanumeric runs to '_' -> strip '_'."""
    return _NORM_RE.sub("_", value.strip().lower()).strip("_")


def display(pop: str) -> str:
    """Human label for plots: 'trinidad_and_tobago' -> 'Trinidad And Tobago'."""
    return " ".join(part.capitalize() for part in pop.split("_") if part)


def _env_populations() -> list[str] | None:
    raw = os.environ.get("BV_POPULATIONS")
    if raw is None:
        return None
    pops = [normalize(p) for p in raw.split(",") if p.strip()]
    if not pops:
        raise SystemExit("BV_POPULATIONS is set but empty after normalization")
    # de-dup, preserve first-seen order
    seen: set[str] = set()
    ordered: list[str] = []
    for p in pops:
        if p not in seen:
            seen.add(p)
            ordered.append(p)
    return ordered


def resolve_populations(raw_dir: Path) -> list[str]:
    """Return the population list, exploding on any missing/empty AF file.

    raw_dir is the directory holding `allele_freq_<pop>.tsv`.
    """
    pops = _env_populations()
    if pops is None:
        discovered = sorted(
            p.name[len("allele_freq_"):-len(".tsv")]
            for p in raw_dir.glob("allele_freq_*.tsv")
        )
        if not discovered:
            raise SystemExit(
                f"No allele_freq_*.tsv files found in {raw_dir} and "
                f"BV_POPULATIONS is unset - nothing to do"
            )
        return discovered

    missing: list[str] = []
    for pop in pops:
        f = raw_dir / f"allele_freq_{pop}.tsv"
        if not f.exists() or f.stat().st_size == 0:
            missing.append(pop)
    if missing:
        raise SystemExit(
            "Country facet expected per-population allele-frequency files that "
            "were not produced (empty/missing). Missing populations: "
            + ", ".join(missing)
            + f"\nLooked in: {raw_dir}\n"
            "This usually means every participant in that country had "
            "unparseable or empty genotype data."
        )
    return pops


def require_columns(df_columns, pops: list[str], where: str) -> None:
    """Explode if any expected population column is absent from a dataframe."""
    missing = [p for p in pops if p not in set(df_columns)]
    if missing:
        raise SystemExit(
            f"{where}: expected population columns missing: "
            + ", ".join(missing)
            + f"\nPresent columns: {list(df_columns)}"
        )
