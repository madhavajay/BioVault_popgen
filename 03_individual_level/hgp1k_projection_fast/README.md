# hgp1k_projection_fast

Joint PCA of BioVault samples and 1000 Genomes high-coverage genotypes on the
shared `tools/locus_map.tsv` loci.

This is intentionally not a gnomAD-loading projection. The gnomAD fast pipeline
projects onto a fixed PCA space where BioVault array data only overlaps a small
fraction of loading variants. This pipeline instead:

1. Reads filtered 1KGP high-coverage VCFs from `data/1kgp_high_coverage/filtered`.
2. Reads BioVault/DDNA/Illumina genotype TXT files via `tools.genotype_normalizer`.
3. Keeps biallelic A/C/G/T SNPs with unambiguous `chrom,pos` keys.
4. Writes per-locus allele frequencies for 1KGP and BioVault samples.
5. Fits PCA on the combined dosage matrix and plots both groups.

Outputs:

```text
pca_scores.tsv
study_pca_projection.tsv
allele_freqs.tsv
qc_report.txt
pca_projection.png
errors.tsv
```

Run locally:

```bash
CHROMOSOMES=1 bash scripts/run_hgp1k_projection.sh <data_dir> <working_dir> <out_dir>
```

Useful environment variables:

```text
HGP1K_VCF_DIR          filtered VCF directory
LOCUS_MAP              locus map TSV
CHROMOSOMES            all, 1, 1,2,3, etc.
HGP1K_METADATA_TSV     optional 1KGP sample metadata with sample_id plus pop labels
N_COMPONENTS           default 10
MIN_AF                 default 0.01
MAX_REF_MISSING        default 0.05
MIN_GS                 default 0.15
MAX_VARIANTS           optional testing cap
```
