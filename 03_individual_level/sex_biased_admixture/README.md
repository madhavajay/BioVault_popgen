# sex_biased_admixture

Sex-biased admixture analysis for Caribbean individuals.

## Background

Caribbean populations exhibit female-biased African ancestry and male-biased
European ancestry — a genetic signature of colonial-era sexual coercion.
The core signal: X-chromosome African ancestry > autosomal African ancestry,
because enslaved African women contributed disproportionately to the maternal
lineage, while Y-chromosome haplogroups skew European.

## What this analysis does

1. Reads all 10 individual GSA genotype files from `01_mock_data_generation/output/`
2. Assigns sex deterministically: sorted sample IDs → first 5 = Female, last 5 = Male
3. Computes per-chromosome statistics (SNP count, heterozygosity rate, BAF mean/variance, LRR mean) directly from the genotype data
4. Estimates ancestry-like components via **NMF (k=3)** on BAF values — autosomes (chr 1–22) and X chromosome treated separately
5. Compares Component 1 on X vs autosomes as a proxy for sex-biased ancestry
6. Produces a 4-panel figure and results table

## Inputs

| File | Description |
|---|---|
| `01_mock_data_generation/output/{id}/*.txt` | Per-individual GSA genotype files (rsid, chrom, pos, genotype, gs, baf, lrr) |

No external reference panels required.

## Outputs

| File | Description |
|---|---|
| `results/sex_bias_results.tsv` | Per-individual table (see columns below) |
| `plots/figure4_sex_biased_admixture.pdf` | 4-panel figure (vector, fonts embedded) |
| `plots/figure4_sex_biased_admixture.png` | 4-panel figure (300 dpi raster) |
| `logs/sex_biased_admixture.log` | Run log |

### Results table columns

| Column | Description |
|---|---|
| `sample_id` | Individual identifier |
| `sex` | Assigned sex (F / M) |
| `auto_snps` | SNP count on autosomes (chr 1–22) |
| `x_snps` | SNP count on X chromosome |
| `auto_het` | Heterozygosity rate — autosomes |
| `x_het` | Heterozygosity rate — X chromosome |
| `auto_lrr_mean` | Mean LRR — autosomes |
| `x_lrr_mean` | Mean LRR — X chromosome |
| `auto_c1/c2/c3` | NMF ancestry components — autosomes (sum to 1) |
| `x_c1/c2/c3` | NMF ancestry components — X chromosome (sum to 1) |
| `x_minus_auto_c1` | X − autosomal Component 1 (the sex-bias signal) |

### Figure panels

| Panel | What it shows |
|---|---|
| **a** | Scatter: Component 1 on autosomes vs X, coloured by sex. Points above the diagonal = X enriched for Component 1. |
| **b** | Lollipop: X − autosomal Component 1 per individual. Positive = X enriched; negative = autosome enriched. |
| **c** | Grouped stacked bars: NMF component proportions on autosomes (A) vs X (X) per individual, split by sex. |
| **d** | Heterozygosity rate on autosomes vs X per individual. In real data, male X het < autosomal het (hemizygous X). |

## How to run

```bash
cd 03_individual_level/sex_biased_admixture
python scripts/sex_biased_admixture.py
```

Requirements: `numpy`, `pandas`, `matplotlib`, `scikit-learn`

## ⚠ Mock data caveat

The current input was generated with `--alt-frequency 0.5` (uniform diploid
signal on all chromosomes). This has two consequences:

1. **Heterozygosity = 0**: every genotype call is homozygous (AA or BB),
   so `auto_het = x_het = 0.0` for all samples. This is correct for this
   synthetic dataset — not a bug.

2. **NMF components are noise**: all samples are genetically indistinguishable,
   so NMF finds no structure and produces arbitrary component weights that do
   not converge. The `x_minus_auto_c1` values are random rather than reflecting
   a real sex-biased signal.

**To obtain meaningful results**: replace the mock files with Carika's real GSA
data (same folder structure, same filename format). No code changes needed.

## Methodology note — why NMF instead of ADMIXTURE

ADMIXTURE requires PLINK binary files and a reference panel (e.g., 1KGP) to
anchor the ancestry components to known populations. NMF(k=3) is used here as
a reference-free alternative that runs entirely on the input BAF matrix. In
real data, the three NMF components will approximate African / European /
Native American axes; in mock data they are uninformative.

When Carika's real data is available, the recommended upgrade path is:

```
real GSA .txt → VCF (biosynth genotype-to-vcf) → PLINK BED
→ merge with 1KGP (02_reference_panels) → ADMIXTURE K=3
→ replace NMF components with supervised ADMIXTURE proportions
```
