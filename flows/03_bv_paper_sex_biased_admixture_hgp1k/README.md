# 03_bv_paper_sex_biased_admixture_hgp1k

HGP1K-anchored **sex-biased admixture** via ADMIXTURE.

Joint ADMIXTURE of the study cohort **+ a frozen 1000 Genomes reference**
(300 AFR + 300 EUR + 300 SAS unrelated founders) at **K = 3, 4, 5**, run:

1. on the **combined genome** — the headline ancestry result, and
2. separately on **autosomes** and **X** — compared, per labelled component, to
   detect **sex-biased admixture** (a component over-represented on X vs
   autosomes indicates a female-biased contribution from that ancestry, since
   X spends ⅔ of its time in females).

This replaces the older NMF approach (`sex_biased_admixture`): ancestry now
comes from **reference-anchored ADMIXTURE** instead of unsupervised NMF.

## How components are made comparable

ADMIXTURE components are anonymous. Each is **labelled by which reference
superpopulation loads highest** on it (AFR / EUR / SAS; extra components for
K>3 get suffixed labels). Because the autosome run and the X run are labelled
against the **same reference**, their components line up — so "the AFR
component on X vs the AFR component on autosomes" is a meaningful comparison.

## Sex / X handling

Males are hemizygous on X. The **`sex` participant facet** (required) is
materialised into `sex_mapping.tsv` and applied with `plink2 --update-sex`, so
PLINK encodes male non-PAR X as **haploid (dosage 0/1)** — for both study and
reference. `--split-par` keeps the pseudo-autosomal regions diploid. Sex is
**never inferred** from the data. The X table includes a per-sex breakdown
(`mean_x_female`, `mean_x_male_haploid`); read male X estimates with the
haploid caveat in mind.

## Inputs

A samplesheet of `GenotypeRecord`s with a **`sex`** facet:

```csv
participant_id,genotype_file,sex
PC0001,/path/PC0001_..._GSAv3....txt,Female
PC0159,/path/PC0159_..._Raw Data....txt,Male
```

`sex` accepts `M/F`, `Male/Female`, or `1/2`.

## Outputs (published to `results_dir`)

| File | What |
|---|---|
| `sex_bias_x_vs_auto.tsv` | **The key result.** Per K, per ancestry: `mean_auto`, `mean_x`, `delta_x_minus_auto`, plus `mean_x_female` / `mean_x_male_haploid`. |
| `sex_bias_x_vs_auto.png/.pdf` | Bar chart of mean study ancestry, autosomes vs X, per K. |
| `component_labels.tsv` | component → ancestry-label map, per run (combined/auto/x) and K. |
| `admixture_{run}_K{k}_labeled_Q.tsv` | per-sample labelled proportions (`group` = study/reference, `superpopulation`, `sex`). |
| `admixture_{run}_K{k}.Q` / `.P` | raw ADMIXTURE matrices (`run` ∈ combined/auto/x). |
| `admixture_{run}_K{k}.log` | ADMIXTURE logs. |
| `qc_report.txt` | per-compartment variant/sample counts, common SNPs, pruned counts, QC params. |
| `plink_bed.tar.gz` | study cohort PLINK BED. |

## Pipeline

```
cohort_bed  (biosynth)          bvs cohort-bed  -> study_raw BED + sex_mapping.tsv
hgp1k_admixture (biovault-admixture)
  study_raw ─ plink1.9 clean + monomorphic ALT==REF→'0'
            ─ plink2 re-key chr:pos, --update-sex, --rm-dup
            ─ split: --chr 1-22 (auto) / --chr X --split-par (X)
  per compartment:
            ─ intersect study ∩ reference on chr:pos
            ─ drop strand-ambiguous (A/T, C/G)
            ─ plink --bmerge  (+ .missnp retry)        [study + reference]
            ─ QC: --geno .05 --mind .10 --maf .01 --hwe 1e-6
            ─ LD prune: --indep-pairwise 200 50 0.2
  combined  = concat(auto_pruned, x_pruned)
  ADMIXTURE K=3,4,5 on combined / auto / x
            ─ label components by reference superpop
            ─ X-vs-autosome comparison + plot
```

QC thresholds and LD parameters match the other pipelines (`pca_qc`,
`find_k`, `gnomad_projection`) and are overridable via env
(`BV_GENO`, `BV_MIND`, `BV_MAF`, `BV_HWE`, `BV_LD_WINDOW/STEP/R2`,
`BV_ADMIXTURE_K`, `BV_THREADS`).

## The baked reference

`reference/reference_samples.tsv` — **frozen** 900-sample subset (300 each
AFR/EUR/SAS), unrelated founders, deterministic (`seed 42`). Regenerate:

```bash
python reference/select_reference_samples.py            # --per-pop, --seed, --include-related
```

The reference **PLINK BED** is pre-generated from the filtered 1KGP VCFs
(`data/1kgp_high_coverage/filtered`, all chr incl. chrX), keyed by `chr:pos`,
pre-split into `reference_auto` / `reference_x`, and baked into the
`biovault-admixture` image. Build it (run inside an image with
bcftools+plink2, e.g. the admixture `tools` stage):

```bash
flows/03_bv_paper_sex_biased_admixture_hgp1k/reference/build_reference_bed.sh
```

### Committed reference (no VCFs needed in CI)

The generated BED is committed (split <100 MB) under
**`data/hgp1k_900_sex_bias/`** — same scheme as `data/hgp1k_dosage_split`:

```
hgp1k_900_sex_bias.tar.gz.aa   # tar.gz of reference_auto/reference_x BED + labels + sample list
hgp1k_900_sex_bias.yaml        # b3sums, reassembly instructions, and the sample-selection metadata
reference_samples.tsv          # which 900 1KGP samples were chosen (sample_id, superpop, population, sex)
```

`build_admixture_docker.sh` resolves the reference in this order:
1. reuse an existing BED in `.docker/`,
2. **reassemble from the committed shards** (`reassemble_reference.sh`) — the CI path,
3. regenerate from the filtered 1KGP VCFs (`build_reference_bed.sh`) — needs the VCFs.

So a clean checkout builds the image with **no raw 1KGP data**. Use
`FORCE_REFERENCE=1` to rebuild from VCFs and `pack_reference.sh` to re-shard.

> Note: this needs **genotype calls incl. chrX** — the projection dosage matrix
> (`hgp1k_dosage.npz`) is autosomes-only and float dosages, so it is **not**
> reusable here.
