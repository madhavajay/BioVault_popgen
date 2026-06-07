# Synthetic sex-biased admixture test set

Generate a cohort with a **known, planted** ancestry difference between
autosomes and X, so you can verify the flow recovers it.

Genotypes are drawn from the **real 1KGP AFR/EUR allele frequencies** (from the
baked reference), with autosomal AFR ≠ X AFR. Because ADMIXTURE is anchored to
the same reference, it recovers the planted proportions. Males are hemizygous
on X (PLINK codes them haploid).

## Generate

```bash
# 100 samples, autosomal AFR=0.5, X AFR=0.85  (expect AFR Δ_x-auto ≈ +0.35)
./generate_testset.sh <out_dir> 100

# tune the planted signal
AFR_AUTO=0.4 AFR_X=0.9 N_AUTO_SNPS=60000 FEMALE_FRAC=0.5 SEED=42 \
  ./generate_testset.sh <out_dir> 300
```

Writes `<out_dir>/samplesheet.csv` (`participant_id, genotype_file, sex`),
per-participant DDNA genotype files, and `GROUND_TRUTH.txt`.

## Run + verify

Point the flow (`flow.yaml` / desktop app) at `<out_dir>/samplesheet.csv`.
Check `sex_bias_x_vs_auto.tsv`:

```
expected:  AFR delta_x_minus_auto ≈ (AFR_X − AFR_AUTO)
           EUR delta ≈ −(AFR_X − AFR_AUTO),  SAS ≈ 0
```

Validated run (24 samples, AFR auto 0.5 / X 0.85): recovered **AFR Δ = +0.33**,
EUR Δ = −0.34, SAS ≈ 0 — i.e. the planted +0.35 signal, cleanly detected.

## Notes
- Needs the baked reference at `.docker/reference/hgp1k_admixture` (the wrapper
  reassembles it from the committed shards if missing).
- A null/negative control: set `AFR_AUTO == AFR_X` → expect all Δ ≈ 0.
- This is **not** the 23andMe/PGP heterogeneity case — synthetic files all carry
  the same SNP panel, so nothing is dropped by `--mind`.
