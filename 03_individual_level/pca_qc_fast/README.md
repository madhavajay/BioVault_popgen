# pca_qc_fast — Fast Genotype QC and Unsupervised PCA

Drop-in faster implementation of `../pca_qc/`.

It writes the same published PCA/QC outputs relative to this directory, but the
large internal matrix is stored as compact on-disk memmaps instead of pandas TSV
matrices by default:

```text
data/merged/snp_info.tsv
data/qc/filtered_snps.tsv
data/plink/genotypes.bed
data/plink/genotypes.bim
data/plink/genotypes.fam
data/plink/genotypes.map
data/pca/pca.eigenvec
data/pca/pca.eigenval
plots/pca_pc1_pc2.png
plots/pca_pc3_pc4.png
logs/fast_pipeline.log
```

`genotypes.ped` is intentionally a placeholder in the fast path because PED is
multi-GB text at cohort scale. Set `BV_WRITE_MATRICES=1` to emit the legacy
numeric genotype TSV for debugging; the raw string matrix remains skipped in
the compact backend.

By default `BV_PCA_BACKEND=auto`: use `plink2` for QC/pruning/PCA when it is on
`PATH`, otherwise fall back to the chunked Python backend. Set
`BV_PCA_BACKEND=python` to force the Python path or `BV_PCA_BACKEND=plink` to
require PLINK.

Parsing is parallel by default (`BV_WORKERS`, `BV_PARSE_MODE=process|thread|serial`).
Parsed samples are cached under `data/work/parsed_samples/` as compact byte
arrays, then reused to build the SNP universe and fill the memmaps without a
second genotype-file parse.

Run:

```bash
cd scripts
bash run_pipeline.sh
```

The implementation keeps the existing dependency set but avoids the slowest
and largest dense-matrix operations in the original pipeline:

- every readable sample contributes to the full SNP universe before QC;
- sample parsing is parallel and cached to compact per-sample arrays;
- alleles and dosages are stored as `uint8`/`int8` memmaps;
- PLINK BED/BIM/FAM is written directly;
- call-rate, MAF, HWE, LD pruning, and PCA run in chunks;
- PCA uses a chunked sample covariance matrix instead of a dense samples x SNPs
  float matrix.
