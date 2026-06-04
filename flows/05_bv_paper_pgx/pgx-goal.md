# PGx Goal

Update `05_bv_paper_pgx` so it runs PharmCAT end to end for BioVault genotype inputs and existing VCF inputs.

## Inputs

- DDNA genotype text files.
- Illumina genotype text files.
- 23andMe/PGP-style genotype text files supported by `bvs genotype-to-vcf`.
- Existing `.vcf`, `.vcf.gz`, `.vcf.bgz`, or `.bgz` files, which should skip BVS conversion.

Each participant record must provide:

- `participant_id`
- `genotype_file`
- `country`

For 1KGP validation, use the 1KGP subpopulation as `country`.
For PGP validation, use `USA` as `country`.

## Processing

1. Prepare a per-participant VCF.
   - Existing VCF files are copied through.
   - Genotype text files are converted with:

     ```bash
     bvs genotype-to-vcf \
       -i sample.txt \
       --sample PARTICIPANT_ID \
       --output PARTICIPANT_ID.vcf.gz \
       --gzip
     ```

2. Run PharmCAT:

   ```bash
   pharmcat_pipeline PARTICIPANT_ID.vcf.gz \
     -o pharmcat_PARTICIPANT_ID \
     -bf PARTICIPANT_ID \
     -reporterHtml \
     -reporterJson \
     -reporterCallsOnlyTsv
   ```

3. Parse PharmCAT reports.
   - Prefer `*.report.json` because it contains one source diplotype object per possible call.
   - Fall back to `*.report.tsv` when JSON is unavailable.
   - Preserve unresolved genotype sets as a single comma-separated value, including match scores:

     ```text
     *1/*2 (2), *1/*35 (2), *1/*61 (2)
     ```

## Required Outputs

- `pgx_participant_possible_genotypes.tsv`
  - `participant_id`
  - `country`
  - `gene`
  - `possible_genotypes`

- `pgx_country_gene_genotype_counts.tsv`
  - `country`
  - `gene`
  - `possible_genotypes`
  - `count`
  - `sample_count`

The count aggregation treats each unresolved PharmCAT genotype set as one categorical value. If three participants have the same unresolved set, that set receives count `3`.

## Validation Ladder

1. Run one 1KGP sample using subpopulation as country.
2. Run ten 1KGP samples using subpopulation as country.
3. Run a small PGP sample from `/Users/madhavajay/dev/snpdata/pgp` using country `USA`.
4. After the small runs pass, run the full 1KGP and PGP datasets and publish results under:

   ```text
   /Users/madhavajay/dev/BioVault_popgen/results
   ```

