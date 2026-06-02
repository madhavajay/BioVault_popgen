# BioVault_popgen

Population-genetics analyses of Caribbean ancestry, organized by **input**: 
individual-level genotype analyses on one side, 
population-level allele-frequency analyses on the other. 
The repo also documents how the input data are generated and which reference panels are used.

This repository is a fork of the original population-genetics analysis
work by `kkkathy211`, adapted for BioVault flows, local Docker execution,
and reproducible synthetic-data testing.

## Directory layout

```
BioVault_popgen/
├── 01_mock_data_generation/   # how to generate synthetic GSA genotypes (biosynth)
├── 02_reference_panels/       # download scripts + .tbi indices for gnomAD / 1KGP
├── 03_individual_level/       # analyses whose input is per-individual genotypes
│   ├── pca_qc/                # genotype QC + Python/PLINK PCA
│   ├── admixture/             # 1KGP-based admixture & local ancestry
│   ├── gnomad_projection/     # PCA projection onto gnomAD HGDP+TGP space
│   └── sex_biased_admixture/  # X-chromosome vs autosome ancestry by sex
├── 04_population_level/       # analyses whose input is per-population allele frequencies
│   ├── fst_islands/           # FST between Caribbean islands
│   └── aims_differential_snps/# ancestry-informative markers + per-island differential SNPs
└── docs/                      # see the shared link I posted in WhatsApp 
```

## Workflow

```
              01_mock_data_generation/output/      02_reference_panels/
              (per-individual GSA TXT files)       (gnomAD, 1KGP — download)
                          │                                │
                ┌─────────┴─────────────┐                  │
                ▼                       ▼                  ▼
        03_individual_level/                     04_population_level/
        - pca_qc                                 - fst_islands
        - admixture        ◄─────────────────────  (uses per-country
        - gnomad_projection ◄────────────────────   allele freq TSVs)
        - sex_biased_admixture                   - aims_differential_snps
                                                   (uses gnomAD AF)
```

Numbers are dependency order, not strict sequence — `03` and `04` are independent and can be run in either order once the inputs in `01`/`02` are in place.

## Quick start

### Docker

The easiest way to get a reproducible environment for the whole repo is the
single Docker image:

```bash
./build_docker.sh
```

By default this builds `biovault-popgen:latest` for `linux/amd64`. You can
override those without editing the script:

```bash
IMAGE=biovault-popgen:dev PLATFORM=linux/amd64 ./build_docker.sh
```

The image includes the conda environment, PLINK/ADMIXTURE/Hail tooling, the GCS
connector needed by Hail, and the baked gnomAD v3.1 PCA loadings:

```bash
/opt/biovault/reference/pca_loadings/loadings_variants.tsv         # variant list (study-side pruning)
/opt/biovault/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht/  # full Hail Table (projection)
```

`build_docker.sh` keeps an inspectable host-side cache at:

```bash
.docker/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht/
.docker/reference/pca_loadings/loadings_variants.tsv
```

If the cache is already populated (Hail's `_SUCCESS` marker present and TSV
non-empty), the build copies it straight into the image — no GCS hit. If the
cache is missing, the build mirrors the public gnomAD Hail Table from
`gs://gcp-public-data--gnomad/release/3.1/pca/` once via `gsutil -m cp -r`
(or `uvx --from gsutil gsutil` if uvx is available but gsutil isn't), then
derives the TSV alongside via Hail and bakes both into the image. The
`gsutil -m` mirror is resumable — rerunning the build skips files already on
disk. Subsequent builds reuse the cache. Projection (`pca_project.py`) reads
from the baked path, so the runtime container needs no network access.

To rebuild the host cache from the remote:

```bash
FORCE_REFERENCE_CACHE=1 ./build_docker.sh
```

Run a pipeline inside the container with:

```bash
./scripts/run_in_docker.sh bash 03_individual_level/pca_qc_fast/scripts/run_pipeline.sh
./scripts/run_gnomad_projection_docker.sh
```

1. **Generate mock genotypes** (or supply your own real GSA files in the same naming scheme):
   ```bash
   cd 01_mock_data_generation/scripts
   bash generate_mock_genotypes.sh
   # outputs land in ../output/{id}/{id}_X_X_GSAv3-DTC_GRCh38-...txt
   ```

2. **Download reference panels** you need (skip the ones you don't):
   ```bash
   cd 02_reference_panels/scripts
   bash download_gnomad_v3_hgdp_tgp.sh   # for gnomad_projection
   bash download_gnomad_v3_sites.sh      # for aims_differential_snps
   bash download_1kgp_high_coverage.sh   # for admixture
   ```

### 1000 Genomes high-coverage VCF prep

The repository-level helper scripts in `tools/` download the full
1000 Genomes high-coverage phased panel and then filter it to the loci in
`tools/locus_map.tsv`.

Download all chromosomes (`1-22` and `X`) plus metadata:

```bash
tools/fetch_1kgp_vcf.sh all
```

The downloader writes to `data/1kgp_high_coverage/` by default. It is
resumable: completed files are skipped by comparing local and remote byte
sizes, and partial downloads continue when possible. Use `--no-verify` to
skip the final sample-label check, or set `OUT_DIR=/path/to/dir` to use a
different storage location.

Filter all downloaded VCFs to the project locus map:

```bash
tools/filter_1kgp_by_locus_map.sh --force --jobs 6 all
```

Filtered VCFs and `.tbi` indexes are written to
`data/1kgp_high_coverage/filtered/`. `--force` replaces existing filtered
outputs. `--jobs` controls how many chromosomes are processed concurrently;
each chromosome also uses `THREADS` bcftools compression threads
(`THREADS=4` by default), so choose `JOBS` with available CPUs and disk I/O
in mind. The filter script checks each input `.tbi` and rebuilds missing or
stale indexes before filtering.

3. **Run any analysis.** Each analysis directory is self-contained with its own `run_pipeline.sh` (or equivalent) and a `README.md` explaining inputs, outputs, and dependencies.

## What each analysis answers

| Analysis | Question | Input granularity |
|---|---|---|
| `pca_qc` | Do the samples cluster sensibly after QC? | Individual genotypes |
| `admixture` | What fraction of each individual's genome comes from AFR / EUR / NAT ancestral populations? | Individual genotypes + 1KGP reference |
| `gnomad_projection` | Where do study samples land in the gnomAD HGDP+TGP PCA space? | Individual genotypes + gnomAD HGDP+TGP |
| `sex_biased_admixture` | Is X-chromosome AFR ancestry > autosomal AFR ancestry (signature of sex-biased colonial admixture)? | Individual genotypes |
| `fst_islands` | How genetically differentiated are the Caribbean islands from each other? | Per-island allele frequencies |
| `aims_differential_snps` | Which SNPs best distinguish Caribbean islands from gnomAD reference populations? | Per-island AF vs gnomAD AF |

## Data not in this repo

The following are **not committed** (see `.gitignore`) but can be regenerated or re-downloaded:

- Synthetic GSA TXT files in `01_mock_data_generation/output/*/` — re-create via the biosynth script.
- gnomAD / 1KGP VCFs in `02_reference_panels/` — re-download via the scripts in `scripts/`.
- Large per-island allele frequency TSVs in `04_population_level/fst_islands/data/` — provided by collaborators.
- PLINK binary files (`*.bed`, `*.bim`, `*.fam`) and Hail tables (`*.ht/`) — produced by the pipelines.

## Repository convention

Each analysis directory follows the same skeleton:

```
<analysis>/
├── README.md     # what the analysis does, how to run, expected I/O
├── scripts/      # numbered scripts (01_*, 02_*, ...) + run_pipeline.sh
├── data/         # primary inputs to this analysis (gitignored if large)
├── results/      # tabular outputs (TSV, PLINK, eigenvecs, etc.)
├── plots/        # figures (PNG, PDF)
└── logs/         # pipeline run logs
```

Not every analysis has every folder 
e.g. `pca_qc` reads its inputs from `01_mock_data_generation/output/`, so it doesn't have its own `data/`.
