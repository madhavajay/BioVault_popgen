# sex_biased_admixture_find_k

Prototype calibration stage for reported-ancestry sex-biased admixture.

The old sex-biased flow used unsupervised NMF components (`c1/c2/c3`) and could
not reliably say which component represented AFR, EUR, or SAS. This stage uses
self-reported single-ancestry participants as anchors after ADMIXTURE:

- participants reporting only African anchor the AFR component
- participants reporting only European anchor the EUR component
- participants reporting only Indian/South Asian anchor the SAS component
- mixed reporters are retained in ADMIXTURE, but are not used to label components

The output is a selected K plus a component label map. A later stage should run
the chosen model for autosomes and X and compare labeled ancestry proportions.

## Inputs

```text
--bed-prefix       PLINK binary prefix for the cohort
--ancestry-map     TSV/CSV with participant_id and ancestry column
--out-dir          output directory
```

Accepted ancestry column names:

```text
self_reported_ancestry
reported_ancestry
ancestry
ethnicity
```

Mixed labels can use delimiters such as `;`, `,`, `/`, `+`, `|`, `and`, or `&`.

## Outputs

```text
admixture_cv_errors.tsv
selected_k.tsv
component_anchor_means.tsv
component_labels.tsv
reported_ancestry_normalized.tsv
ancestry_anchor_samples.tsv
admixture_K{K}.Q
admixture_K{K}.P
admixture_K{K}_labeled_Q.tsv
admixture_k_summary.tsv
admixture_k_summary.png
find_k_report.txt
```

`component_labels.tsv` is the handoff file for the next stage.

## Experiment Mode

For quick ADMIXTURE timing runs:

```bash
BV_ADMIXTURE_K_MIN=3 \
BV_ADMIXTURE_K_MAX=6 \
BV_ADMIXTURE_REPS=3 \
BV_THREADS=8 \
python scripts/find_k.py \
  --bed-prefix /path/to/plink_prefix \
  --ancestry-map /path/to/reported_ancestry.tsv \
  --out-dir /tmp/find_k_test
```

If the BED has already been filtered/pruned and should be used directly:

```bash
python scripts/find_k.py \
  --bed-prefix /path/to/prepared_prefix \
  --ancestry-map /path/to/reported_ancestry.tsv \
  --out-dir /tmp/find_k_test \
  --k-min 3 \
  --k-max 6 \
  --reps 3 \
  --skip-plink-prep
```
