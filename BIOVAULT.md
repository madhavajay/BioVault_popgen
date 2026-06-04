# BioVault popgen — End-to-end pipeline

A researcher with a cohort of GSAv3 genotype files (DDNA or Illumina,
one per participant) walks the following steps to go from raw uploads to
sample QC, ancestry projection, sex-biased admixture checks, and ranked
ancestry-informative SNPs (AIMs) distinguishing sub-populations in the
cohort. Each step runs as a self-contained Nextflow flow.

Production flows use the slim fast image:

```
ghcr.io/madhavajay/biovault-popgen:0.2.1-fast
```

The full image is still built as `ghcr.io/madhavajay/biovault-popgen:0.2.1` for slow/Hail
reference paths and debugging. Tags are mutable: a later CI build with
the same `VERSION` will overwrite `0.1.1` and `0.1.1-fast`. Use
`sha-<short>` / `sha-<short>-fast` tags when a frozen production image is
required.

Attribution: the original population-genetics analysis work was done by
`kkkathy211`. The fast paths and BioVault flows keep the same analysis
intent while adapting the implementation for local Docker and Nextflow
execution.

The numbered analysis flows are independent — each consumes the same
`List[GenotypeRecord]` samplesheet and none chains its outputs into
another. The order below is the recommended workflow, not a data
dependency. The only real data dependency is *inside* flow 4, where the
per-country split → FST → AIMs sub-steps chain.

```
        participants                participant facet
       (DDNA/Illumina TXTs)         country (required → 4)
                                    sex     (required → 3)
            │                            │
            ▼                            │
  ┌───────────────────────────┐          │   (1) sanity check
  │ bv_paper_pca_qc_fast      │          │   ───────────────
  │   within-cohort PCA       │          │   no reference;
  │   no reference panel      │          │   confirms data
  └───────────────────────────┘          │   parses + clusters
            │                            │
            ▼                            │
  ┌───────────────────────────┐          │   (2) ancestry
  │ bv_paper_gnomad_          │          │   projection
  │       projection_fast     │          │   ──────────
  │   project onto gnomAD     │          │   per-sample
  │   HGDP+1kGP 30-PC space   │          │   PC1..PC30
  └───────────────────────────┘          │
            │
            ▼                                (3) sex bias
  ┌───────────────────────────┐             ─────────────
  │ bv_paper_sex_biased_      │             sex facet →
  │       admixture_fast      │             autosome vs X
  │   X-hemizygosity signal   │             heterozygosity
  └───────────────────────────┘
                                         │
            ┌────────────────────────────┘
            ▼                                 (4) FST + AIMs
  ┌─────────────────────────────────────┐    ──────────────
  │ bv_paper_population_level           │    one flow, two
  │  split by country (biosynth)        │    containers:
  │   → FST (WC84 matrix + plots)       │    per-country AF →
  │   → AIMs (vs gnomAD ref + panels)   │    FST → AIMs SNPs
  └─────────────────────────────────────┘
```

## Step 0 — Data upload

Researcher hands BioVault two things:

- A per-participant directory layout of GSAv3 TXT files, one TXT per
  participant (e.g. `100383/100383_X_X_GSAv3-DTC_GRCh38-…txt`).
- A country/island label per participant, exposed to BioVault as the
  `country` participant facet.
- A sex label per participant, exposed to BioVault as the `sex`
  participant facet.

Flows consume both via a `List[GenotypeRecord]` samplesheet
(participant_id + genotype_file path + requested facets).

## Reference prep — 1000 Genomes high-coverage VCFs

Use the repository-level helpers when the local 1KGP high-coverage panel
needs to be created or refreshed. The fetch step downloads chromosomes
`1-22` and `X` plus metadata into `data/1kgp_high_coverage/` by default:

```bash
tools/fetch_1kgp_vcf.sh all
```

The downloader is resumable. It skips files whose local byte size matches
the remote file and continues partial downloads when possible. Use
`OUT_DIR=/path/to/dir` to change the download location, `PARALLEL=N` to
change concurrent downloads or aria2 split count, and `--no-verify` to
skip the final VCF sample-label verification.

After the VCFs are present, filter every chromosome to the loci in
`tools/locus_map.tsv`:

```bash
tools/filter_1kgp_by_locus_map.sh --force --jobs 6 all
```

The filtered VCFs and indexes are written to
`data/1kgp_high_coverage/filtered/`. `--force` replaces existing filtered
outputs. `--jobs` controls chromosome-level parallelism; each chromosome
also uses `THREADS` bcftools compression threads (`THREADS=4` by default).
The filter script checks each input `.tbi` index and rebuilds it if it is
missing or older than the VCF.

To rerun a subset, pass explicit chromosomes:

```bash
tools/filter_1kgp_by_locus_map.sh --force --jobs 6 3 4 5 X
```

## Step 1 — `bv_paper_pca_qc_fast` (sanity check)

**What it does.** Cohort-internal PCA over the SNPs that pass GENO,
MIND, MAF, HWE, and LD-prune filters. No reference panel, no projection
— purely "do my N samples look sensible, do they cluster?".

**Why first.** Catches data-loading bugs (wrong build, malformed TXTs,
missingness spikes) before you spend the run-time on the reference
projection step.

**Outputs**
- `pca.eigenvec`, `pca.eigenval` — cohort PC scores and variance
  explained per PC.
- `snp_info.tsv` — SNPs that survived QC.
- `pca_pc1_pc2.png`, `pca_pc3_pc4.png` — quick-look scatter plots.
- `fast_pipeline.log` — per-step timing and per-filter dropped counts.

Source: `03_individual_level/pca_qc_fast/scripts/fast_pipeline.py`
(baked into the image at
`/opt/biovault/scripts/pca_qc_fast/`).

## Step 2 — `bv_paper_gnomad_projection_fast` (ancestry projection)

**What it does.** Projects each participant onto the gnomAD v3.1
HGDP+1kGP PCA reference (30 PCs). Same I/O contract as the older
`biovault-popgen-gnomad-projection-1` flow, same numerical output, but
runs ~2× faster at N=10 and ~7× at N=1000 (see "Performance" below).

**How it differs from the slow flow under the hood**

| Stage          | Slow                                  | Fast                                          |
|----------------|---------------------------------------|-----------------------------------------------|
| DDNA → PLINK   | Python row-loop tped + plink2 --make-bed | One vectorized pass; bed/bim/fam written directly |
| PC projection  | `hl.experimental.pc_project` (~45 s JVM warm-up) | Numpy matmul (<1 s) reproducing Hail's formula |

The fast flow's `study_pca_projection.tsv` is **bit-identical** to the
slow flow's at float64 precision (max abs diff ~1e-17 across all tested
N).

**Outputs**
- `study_pca_projection.tsv` — one row per sample, 30 PC scores in
  HGDP+1kGP space.
- `qc_report.txt` — SNPs in/out, samples in/out, applied QC thresholds.
- `pca_projection.png` — PC1 vs PC2 scatter of the projected cohort.

Source: `03_individual_level/gnomad_projection_fast/scripts/`
(baked into the image at
`/opt/biovault/scripts/gnomad_projection_fast/`).

### Performance (mock data, biosynth `--alt-frequency 0.5`)

| N samples | Slow flow  | Fast flow | Speedup |
|-----------|------------|-----------|---------|
| 10        | 64 s       | 33 s      | 1.9×    |
| 50        | 107 s      | 35 s      | 3.1×    |
| 100       | ~167 s     | 44 s      | 3.8×    |
| 1000      | ~1500 s    | 206 s     | 7.3×    |

At N=1000 the dominant cost is the per-sample DDNA `pd.read_csv` (~0.5 s
per file × 1000 / 8 cores ≈ 60 s), then per-variant aggregation in
numpy, then PLINK 2's downstream QC (`--rm-dup`, `--geno`, `--mind`,
`--maf`, `--hwe`) which takes ~30 s for 396k SNPs × 1000 samples.

## Step 3 — `bv_paper_sex_biased_admixture_fast`

**What it does.** Computes per-sample autosomal and X-chromosome
heterozygosity and NMF components, then contrasts those values by the
declared `sex` participant facet. In synthetic positive-control data,
male samples have lower X heterozygosity than female samples while
autosomal heterozygosity stays comparable.

**Input.** Same `List[GenotypeRecord]` samplesheet as flows 1–2, plus a
required `sex` facet (`module.yaml` declares `required_facets: [sex]`).
The flow materializes that facet into `sex_mapping.tsv` inside the
process and passes it via `BIOVAULT_SEX_MAPPING`. The analysis never
infers sex from genotype data.

**Outputs**
- `sex_bias_results.tsv` — per-participant autosomal/X heterozygosity
  and NMF components.
- `figure4_sex_biased_admixture.png` — Figure 4 plot.
- `figure4_sex_biased_admixture.pdf` — Figure 4 PDF.
- `sex_biased_admixture.log` — sample counts and sex-facet labels loaded.

Source: `03_individual_level/sex_biased_admixture_fast/scripts/` plus
the original reusable analysis in
`03_individual_level/sex_biased_admixture/scripts/` (both baked into the
fast image). The original analysis work is credited to `kkkathy211`.

## Step 4 — `bv_paper_population_level` (population FST + AIMs)

**What it does.** A single flow that wraps per-country allele frequency
generation, FST, and AIMs into one pipeline, driven by a `country`
participant facet instead of a static `island_mapping.tsv`.

**Input.** Same `List[GenotypeRecord]` samplesheet as flows 1–2, plus a
required `country` facet (`flow.yaml` declares `required_facets:
[country]`). The desktop refuses to generate the samplesheet if any
selected participant has an empty `country`, so a hole can't reach the
flow.

**Internal stages** (two containers, orchestrated by Nextflow — the
desktop runner pre-pulls every per-process `container` it finds):

1. **Split by country** — `container ghcr.io/openmined/biosynth:0.1.31`.
   Per-participant `bvs emit-long`, then per-country
   `bvs aggregate-long` → `allele_freq_<country>.tsv`. The country label
   is the facet value normalized: trim → lowercase → non-alphanumeric
   runs → `_` → strip `_` (e.g. `"Trinidad and Tobago"` →
   `allele_freq_trinidad_and_tobago.tsv`). The Groovy normalizer in
   `main.nf` and `scripts/popset.py` are kept identical.
2. **FST** — `container ghcr.io/madhavajay/biovault-popgen:0.2.1-fast`. Load/merge
   per-country AF → pairwise Weir & Cockerham 1984 matrix → heatmap /
   dendrogram / population PCA.
3. **AIMs** — same container. Merge against the bundled gnomAD HGDP+TGP
   reference, per-population differential SNPs vs gnomAD AFR/global, and
   AFR/NFE + AFR/SAS AIMs panels with a combined PCA.

**Fail-loud, two layers.** Desktop `validate_required_facets` rejects a
samplesheet missing `country`; inside the flow the split aborts if any
country yields zero usable genotypes, and every popgen script
re-asserts its expected `allele_freq_<pop>.tsv` exists and is non-empty
(`popset.resolve_populations`).

**Outputs** (published to `params.results_dir`)
- `country_map.tsv` — participant → normalized country.
- `fst_matrix.tsv` — pairwise WC84 FST.
- `merged_allele_freq_annotated.tsv` — per-locus AF, all populations.
- `master_af_table.tsv` — populations + gnomAD global/AFR/NFE/SAS.
- `all_outliers_long.tsv` — long-format per-population differential SNPs.
- `aims_combined.tsv` — deduplicated AFR/NFE + AFR/SAS AIMs panel.
- `population_level_summary.txt` — FST matrix + master-AF summary.
- FST and AIMs plots (`*.png` / `*.pdf`).

Source: `04_population_level/fst_aims_fast/scripts/`, baked into the
image at `/opt/biovault/scripts/population_level/`. The original FST and
AIMs analysis work is credited to `kkkathy211`.

The AIMs gnomAD reference is **pre-baked at image-build time** (mirroring
the PCA-loadings pattern): `build_docker.sh` mirrors the ~80 GB HGDP+TGP
VCFs once and `build/derive_gnomad_aims_af.py` derives the small
`gnomad_af_per_locus.tsv`, which is the only file COPYed into the runtime
image (at `/opt/biovault/reference/aims/`). Runtime never touches the
VCFs.

## Step 5 — `bv_paper_pgx` (PharmCAT pharmacogenomics)

**What it does.** Converts each selected participant genotype TXT to a
compressed VCF with biosynth, runs the Docker PharmCAT pipeline on each
VCF, and aggregates per-gene PGx calls by `country` and `sex`. The
summary plots treat non-reference / non-wildtype PharmCAT diplotypes as a
PGx call-burden signal, since PharmCAT report TSVs do not include
population allele frequencies.

**Input.** Same `List[GenotypeRecord]` samplesheet as the other flows,
plus required `country` and `sex` facets.

**Internal stages** (two containers, orchestrated by Nextflow):

1. **Genotype TXT → VCF** — `container ghcr.io/openmined/biosynth:0.1.31`.
   Runs `bvs genotype-to-vcf -i <genotype.txt> --output <id>.vcf.gz --gzip`
   per participant.
2. **PharmCAT** — `container pgkb/pharmcat`. Runs `pharmcat_pipeline`
   with HTML, JSON, and calls-only TSV reports enabled.
3. **Aggregation** — adds `country` and `sex` facets to every PharmCAT
   gene row, then writes both country+sex counts and country-only counts.
4. **Visualization** — `04_population_level/pgx/plot_pgx_accumulation.py`
   writes country and country+sex heatmaps showing where non-reference PGx
   calls accumulate by gene, plus a top-gene burden bar chart.

**Outputs**
- `pgx_participant_results.tsv` — one row per participant and PharmCAT
  gene result, with `country` and `sex` columns.
- `pgx_country_sex_summary.tsv` — counts grouped by country, sex, gene,
  diplotype, and phenotype.
- `pgx_country_summary.tsv` — same counts grouped by country only,
  summing across sex.
- `pgx_participant_manifest.tsv` — successful participant/facet manifest.
- `pgx_gene_country_burden.tsv` — non-reference PGx call rate/count by
  country and gene.
- `pgx_gene_country_sex_burden.tsv` — non-reference PGx call rate/count
  by country, sex, and gene.
- `pgx_plots/*` — heatmaps and the top-gene burden chart as PNG/PDF.
- `vcfs/*.vcf.gz` — per-participant compressed VCFs.
- `pharmcat_reports/*.report.tsv`, `*.report.json`, `*.report.html` —
  per-participant PharmCAT reports.

## Declared outputs and raw-data boundary

The flows publish only declared headline artifacts from each process
root. Nextflow work directories and local CLI `results/*/work/` folders
can contain raw/intermediate files for debugging, but those are not
declared module outputs.

Declared flow outputs are:

- `01_bv_paper_pca_qc_fast`: `pca.eigenvec`, `pca.eigenval`,
  `snp_info.tsv`, `pca_pc1_pc2.png`, `pca_pc3_pc4.png`,
  `fast_pipeline.log`.
- `02_bv_paper_gnomad_projection_fast`: `study_pca_projection.tsv`,
  `qc_report.txt`, `pca_projection.png`.
- `03_bv_paper_sex_biased_admixture_fast`: `sex_bias_results.tsv`,
  `figure4_sex_biased_admixture.png`,
  `figure4_sex_biased_admixture.pdf`, `sex_biased_admixture.log`.
- `04_bv_paper_population_level`: `country_map.tsv`, `fst_matrix.tsv`,
  `merged_allele_freq_annotated.tsv`, `master_af_table.tsv`,
  `all_outliers_long.tsv`, `aims_combined.tsv`,
  `population_level_summary.txt`, and FST/AIMs plots.
- `05_bv_paper_pgx`: `pgx_participant_results.tsv`,
  `pgx_country_sex_summary.tsv`, `pgx_country_summary.tsv`,
  `pgx_participant_manifest.tsv`, `pgx_gene_country_burden.tsv`,
  `pgx_gene_country_sex_burden.tsv`, `pgx_plots/*`, `vcfs/*.vcf.gz`,
  `pharmcat_reports/*.report.tsv`, `pharmcat_reports/*.report.json`,
  `pharmcat_reports/*.report.html`, `errors.tsv`, `warnings.tsv`,
  `pgx_pipeline.log`.

None of these are raw participant genotype files or PLINK BED/BIM/FAM
intermediates. `merged_allele_freq_annotated.tsv` is an aggregate
per-locus allele-frequency table. `country_map.tsv` is participant facet
metadata, not genotype data.

## Image build

The repo now builds two runtime images from the same Dockerfile:

### Fast production image

`ghcr.io/madhavajay/biovault-popgen:0.2.1-fast` bakes:

- A smaller `biovault_popgen` conda env: Python, PLINK 2, pandas, numpy,
  scipy, scikit-learn, matplotlib-base, seaborn.
- `loadings.npz` under `/opt/biovault/reference/pca_loadings/`.
- The pre-baked AIMs gnomAD reference TSV under
  `/opt/biovault/reference/aims/gnomad_af_per_locus.tsv`.
- Fast production script trees under `/opt/biovault/scripts/`:
  `gnomad_projection_fast/`, `pca_qc_fast/`,
  `sex_biased_admixture_fast/`, `sex_biased_admixture/`, and
  `population_level/`.

The fast image intentionally omits Hail, OpenJDK, bcftools/samtools,
admixture, the expanded gnomAD Hail Table, and the slow script tree. It
is the image referenced by the BioVault flow files.

### Full debug/reference image

`ghcr.io/madhavajay/biovault-popgen:0.2.1` bakes:

- Full `biovault_popgen` conda env (PLINK/PLINK 2, bcftools, samtools,
  htslib, Hail/Spark, pandas, numpy, sklearn, plotting packages).
- gnomAD v3.1 PCA loadings HT, `loadings_variants.tsv`, and
  `loadings.npz` under `/opt/biovault/reference/pca_loadings/`.
- The pre-baked AIMs gnomAD reference TSV under
  `/opt/biovault/reference/aims/gnomad_af_per_locus.tsv`.
- Full script trees under `/opt/biovault/scripts/`:
  - `gnomad_projection/` — slow (Hail) reference implementation.
  - `gnomad_projection_fast/` — numpy-only fast implementation.
  - `pca_qc_fast/` — within-cohort QC.
  - `sex_biased_admixture/` and `sex_biased_admixture_fast/`.
  - `population_level/` — population FST + AIMs.

Build both images locally:

```bash
VERSION=0.2.1 ./build_docker.sh
```

This builds the full image first, then the fast image. It produces:

```text
ghcr.io/madhavajay/biovault-popgen:0.2.1
biovault-popgen:latest
ghcr.io/madhavajay/biovault-popgen:0.2.1-fast
biovault-popgen:fast
```

Build only the full image and skip the fast image:

```bash
VERSION=0.1.2 BUILD_FAST=0 ./build_docker.sh
```

Build only the fast target directly, assuming the local reference cache
already exists:

```bash
docker build --platform linux/amd64 \
  --target fast-runtime \
  -t ghcr.io/madhavajay/biovault-popgen:0.2.1-fast \
  -t biovault-popgen:fast \
  .
```

The HT cache lives on the host at
`.docker/reference/pca_loadings/` and is populated by `gsutil` on first
build; the AIMs AF cache lives at `.docker/reference/aims/` and is
derived from a one-time ~80 GB HGDP+TGP VCF mirror. Both caches are
idempotent — subsequent builds reuse them. A fresh clone can build if the
committed `.docker/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht.tar.gz`,
`.docker/reference/pca_loadings/loadings.npz`, and
`.docker/reference/aims/gnomad_af_per_locus.tsv` cache files are present.
Without those, the full build needs network access and `gsutil`/`uvx` to
prime the reference cache.

CI publishes the mutable version tags plus SHA tags:

```text
ghcr.io/madhavajay/biovault-popgen:0.2.1
ghcr.io/madhavajay/biovault-popgen:0.2.1-fast
ghcr.io/madhavajay/biovault-popgen:latest
ghcr.io/madhavajay/biovault-popgen:fast
ghcr.io/madhavajay/biovault-popgen:sha-<short>
ghcr.io/madhavajay/biovault-popgen:sha-<short>-fast
```

Use a SHA tag for a deployment that must never move.

## Local test scripts

The local scripts run the same script trees as the BioVault flows, but
from the repo against files in `01_mock_data_generation/output`. Use a
small participant count first so failures are quick.

Fast production image:

```bash
docker tag biovault-popgen:fast ghcr.io/madhavajay/biovault-popgen:0.2.1-fast
IMAGE=ghcr.io/madhavajay/biovault-popgen:0.2.1-fast ./03_individual_level.sh --qc 3
IMAGE=ghcr.io/madhavajay/biovault-popgen:0.2.1-fast ./03_individual_level.sh --fast 3
IMAGE=ghcr.io/madhavajay/biovault-popgen:0.2.1-fast ./03_individual_level.sh --sex 4
IMAGE=ghcr.io/madhavajay/biovault-popgen:0.2.1-fast ./04_population_level.sh --limit 1
```

Expected output roots:

- `./03_individual_level.sh --qc 3` →
  `results/pca_qc_fast/`.
- `./03_individual_level.sh --fast 3` →
  `results/gnomad_projection_fast/`.
- `./03_individual_level.sh --sex 4` →
  `results/sex_biased_admixture_fast/`.
- `./04_population_level.sh --limit 1` →
  `results/population_level/`.

Slow/reference paths using the full image:

```bash
IMAGE=ghcr.io/madhavajay/biovault-popgen:0.2.1 ./03_individual_level.sh --slow 3
IMAGE=ghcr.io/madhavajay/biovault-popgen:0.2.1 ./03_individual_level.sh --qc --slow 3
IMAGE=ghcr.io/madhavajay/biovault-popgen:0.2.1 ./03_individual_level.sh --sex --slow 4
IMAGE=ghcr.io/madhavajay/biovault-popgen:0.2.1 ./04_population_level.sh --slow --limit 1
```

Expected slow output roots:

- `./03_individual_level.sh --slow 3` →
  `results/gnomad_projection/`.
- `./03_individual_level.sh --qc --slow 3` →
  `results/pca_qc/`.
- `./03_individual_level.sh --sex --slow 4` →
  `results/sex_biased_admixture/`.
- `./04_population_level.sh --slow --limit 1` →
  `results/population_level_slow/`.

## Flow inventory

| Flow                                       | Step | Image |
|--------------------------------------------|------|-------|
| `01_bv_paper_pca_qc_fast`                  | 1    | `ghcr.io/madhavajay/biovault-popgen:0.2.1-fast` |
| `02_bv_paper_gnomad_projection_fast`       | 2    | `ghcr.io/madhavajay/biovault-popgen:0.2.1-fast` |
| `03_bv_paper_sex_biased_admixture_fast`    | 3    | `ghcr.io/madhavajay/biovault-popgen:0.2.1-fast` |
| `04_bv_paper_population_level`             | 4    | `ghcr.io/openmined/biosynth:0.1.31` + `ghcr.io/madhavajay/biovault-popgen:0.2.1-fast` |
| `05_bv_paper_pgx`                          | 5    | `ghcr.io/openmined/biosynth:0.1.31` + `pgkb/pharmcat` + `ghcr.io/madhavajay/biovault-popgen:0.2.1-fast` |

Each flow lives at `flows/<name>/` with `flow.yaml`, `module.yaml`, and
`main.nf`. Inputs are a `List[GenotypeRecord]` samplesheet; outputs are
declared via `publishDir` to `params.results_dir`.
