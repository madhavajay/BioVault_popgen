# pca_qc_fast — Fast Genotype QC and Unsupervised PCA

Drop-in faster implementation of `../pca_qc/`.

It writes the same output paths relative to this directory:

```text
data/merged/genotype_matrix_raw.tsv
data/merged/genotype_matrix_numeric.tsv
data/merged/snp_info.tsv
data/plink/genotypes.ped
data/plink/genotypes.map
data/pca/pca.eigenvec
data/pca/pca.eigenval
plots/pca_pc1_pc2.png
plots/pca_pc3_pc4.png
logs/fast_pipeline.log
```

Run:

```bash
cd scripts
bash run_pipeline.sh
```

The implementation keeps the existing dependency set but avoids the slowest
Python row loops in the original pipeline:

- sample files are read in parallel;
- allele counting and dosage encoding are vectorized with NumPy;
- HWE filtering is vectorized;
- LD pruning keeps the original greedy window behavior but vectorizes each
  window comparison;
- PED/MAP text output is emitted with vectorized string preparation.
