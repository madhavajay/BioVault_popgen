# 05_cnv

CNV calling workflow for Illumina SNP-array intensity data.

This workflow is split into three stages:

1. IAAP / GenomeStudio export from raw `*.idat` plus `*.bpm` plus `*.egt`.
2. Conversion into PennCNV signal text with `Log R Ratio` and `B Allele Frequency`.
3. PennCNV signal splitting, PFB creation, and CNV calling.

PennCNV is available in Docker. IAAP CLI is not redistributed here; mount a trusted Illumina IAAP CLI binary into the container. The documented IAAP path produces `.gtc` files via `gencall`; those must then be converted to LRR/BAF signal text, for example with `bcftools +gtc2vcf` or GenomeStudio.

## Build

```bash
docker build --platform linux/amd64 -t biovault-cnv:local 05_cnv
docker build --platform linux/amd64 -f 05_cnv/Dockerfile.gtc2vcf -t biovault-gtc2vcf:local 05_cnv
```

## Make GTC Files Without IAAP

The `freeseek/gtc2vcf` plugin includes `bcftools +idat2gtc`, which can call GTC files directly from IDAT+BPM+EGT:

```bash
IDAT_DIR=/path/to/idats \
OUT_DIR=/path/to/results/cnv/gtc2vcf \
bash 05_cnv/scripts/run_gtc2vcf_docker.sh
```

This writes GTC files to `OUT_DIR/gtc`.

## Run from an Existing Final Report

Use this if IAAP or GenomeStudio has already exported a Final Report containing `SNP Name`, `Chr`, `Position`, `Sample ID`, `B Allele Freq`, and `Log R Ratio`.

```bash
docker run --rm --platform linux/amd64 \
  -v /path/to/idats:/input:ro \
  -v /path/to/results:/out \
  -e STAGE=penncnv \
  -e FINAL_REPORT=/out/iaap/final_report.txt \
  biovault-cnv:local
```

Outputs are written under `/out/penncnv`.

## Run with IAAP

Mount a licensed/trusted IAAP CLI binary into `/opt/illumina/bin/iaap-cli`. The default command shape is:

```bash
CLR_ICU_VERSION_OVERRIDE="$(uconv -V | sed 's/.* //g')" LANG="en_US.UTF-8" /opt/illumina/bin/iaap-cli \
  gencall \
  /input/HumanOmniExpress-12v1_H.bpm \
  /input/HumanOmniExpress-12v1_H.egt \
  /out/iaap/gtc \
  --idat-folder /input \
  --output-gtc \
  --gender-estimate-call-rate-threshold 0.0
```

If this needs adjustment for a different IAAP build, pass `IAAP_COMMAND`. The runner expands these placeholders before execution:

- `{IAAP_BIN}`
- `{IDAT_DIR}`
- `{BPM}`
- `{EGT}`
- `{SAMPLE_MANIFEST}`
- `{OUT_DIR}`
- `{FINAL_REPORT}`

Example:

```bash
docker run --rm --platform linux/amd64 \
  -v /path/to/iaap-cli:/opt/illumina/bin/iaap-cli:ro \
  -v /path/to/idats:/input:ro \
  -v /path/to/results/cnv:/out \
  -e STAGE=iaap \
  -e IAAP_COMMAND='CLR_ICU_VERSION_OVERRIDE="$(uconv -V | sed '\''s/.* //g'\'')" LANG="en_US.UTF-8" {IAAP_BIN} gencall "{BPM}" "{EGT}" "{OUT_DIR}/gtc" --idat-folder "{IDAT_DIR}" --output-gtc --gender-estimate-call-rate-threshold 0.0' \
  biovault-cnv:local
```

Known IAAP caveats from gtc2vcf documentation:

- `LANG` should be `en_US.UTF-8`, otherwise malformed GTC files can be produced.
- IDAT filenames should be `BARCODE_POSITION_(Red|Grn).idat`; filenames with more than two underscores can fail.
- IAAP CLI v1.1 is documented by Illumina as a native Linux CLI, but the normal support download can require browser/EULA flow.

## Local Convenience Runner

```bash
IDAT_DIR=/path/to/idats \
OUT_DIR=/path/to/results/cnv \
FINAL_REPORT=/path/to/results/cnv/iaap/final_report.txt \
STAGE=penncnv \
bash 05_cnv/scripts/run_docker.sh
```

## Important Inputs

The pipeline auto-detects one `.bpm`, one `.egt`, and the chip-to-sample mapping file from `IDAT_DIR` unless `BPM`, `EGT`, or `SAMPLE_MAP` are set explicitly.

PennCNV autosomal calling uses:

- HMM: `/opt/PennCNV-1.0.5/lib/hhall.hmm`
- PFB: compiled from the cohort signal files
- Output: `/out/penncnv/cnv_calls.rawcnv`

If the Final Report does not include `Chr` and `Position`, set `SNP_POS_FILE` to a PennCNV-compatible SNP position/PFB file so `compile_pfb.pl` can annotate marker coordinates.

Set `CALL_CHRX=1` and provide `SEX_FILE` to add a separate chrX CNV call.
