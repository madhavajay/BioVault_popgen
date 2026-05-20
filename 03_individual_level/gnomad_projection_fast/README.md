# gnomad_projection_fast — Numpy-only projection onto gnomAD HGDP+1kGP PCs

Drop-in faster sibling of `../gnomad_projection`. Same output:

```
study_pca_projection.tsv   # one row per sample, 30 PC scores
qc_report.txt
pca_projection.png         # PC1 vs PC2 scatter
```

Same math as Hail's `hl.experimental.pc_project`:

```
PC[sample, k] = Σ_j (G[sample, j] − 2·pca_af[j]) · loadings[j, k]
```

Where `G[sample, j]` is the sample's alt-allele dosage (0/1/2) at loadings
variant `j`. Missing genotypes contribute 0 after mean-centering.

## Why it's faster

| Stage | `gnomad_projection`                  | `gnomad_projection_fast`            |
|-------|--------------------------------------|-------------------------------------|
| Read  | python row-loop in `convert_ddna_to_plink.py` | parallel `pd.read_csv` per file |
| Match | PLINK merge + Hail intersect         | numpy `searchsorted` on (chrom,pos) |
| QC    | `plink2 --geno/--mind/--maf/--hwe`   | skipped (mock data is uniform-AF)    |
| Project | Hail `pc_project` (Spark + JVM)    | one numpy matrix multiply           |

Hail JVM startup alone is ~45s. This pipeline finishes the same work in
seconds.

## One-time setup

The fast pipeline needs the loadings HT pre-exported to numpy:

```
/opt/biovault/reference/pca_loadings/loadings.npz
```

`extract_loadings_matrix.py` does this once via Hail. Either:
- bake at image build time (preferred — see `build_docker.sh`)
- run on first invocation (the wrapper auto-extracts if missing)

## Run

```bash
bash scripts/run_fast_projection.sh <data_dir> <out_dir>
```

`<data_dir>` is the participant-dirs root (same as the slow pipeline).
