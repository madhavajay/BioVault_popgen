#!/usr/bin/env bash
set -euo pipefail

COLLECTION_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
BASE_URL="${COLLECTION_URL}/working/20220422_3202_phased_SNV_INDEL_SV"
POPULATION_METADATA_URL="${COLLECTION_URL}/20130606_g1k_3202_samples_ped_population.txt"
PEDIGREE_METADATA_URL="${COLLECTION_URL}/working/1kGP.3202_samples.pedigree_info.txt"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUT_DIR="${OUT_DIR:-${REPO_ROOT}/data/1kgp_high_coverage}"
PARALLEL="${PARALLEL:-4}"
TMP_FILES=()

cleanup() {
  local file
  for file in "${TMP_FILES[@]:-}"; do
    rm -f "$file"
  done
}
trap cleanup EXIT

usage() {
  cat <<'USAGE'
Usage:
  tools/fetch_1kgp_vcf.sh [--dry-run] [--no-verify] [chr ...]

Examples:
  tools/fetch_1kgp_vcf.sh          # download chr1 VCF + TBI + metadata, then verify labels
  tools/fetch_1kgp_vcf.sh --dry-run 1
  tools/fetch_1kgp_vcf.sh 1        # same as default
  tools/fetch_1kgp_vcf.sh 1 2 X    # download selected chromosomes
  tools/fetch_1kgp_vcf.sh all      # download chr1-22 and chrX

Environment:
  OUT_DIR=/path/to/data      Output directory. Default: ./data/1kgp_high_coverage
  PARALLEL=4                 Parallel downloads for curl fallback or aria2 split count.

The downloader is resumable. It uses aria2c when available, otherwise curl -C -.
Completed files are skipped by comparing local and remote byte sizes.
USAGE
}

normalize_chr() {
  local chr="$1"
  chr="${chr#chr}"
  case "$chr" in
    [1-9]|1[0-9]|2[0-2]|X) printf '%s\n' "$chr" ;;
    *) printf 'Unsupported chromosome: %s\n' "$1" >&2; return 1 ;;
  esac
}

vcf_name_for_chr() {
  local chr="$1"
  if [[ "$chr" == "X" ]]; then
    printf '1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz\n'
  else
    printf '1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz\n' "$chr"
  fi
}

build_url_list() {
  local chr filename
  for chr in "$@"; do
    filename="$(vcf_name_for_chr "$chr")"
    printf '%s/%s\n' "$BASE_URL" "$filename"
    printf '%s/%s.tbi\n' "$BASE_URL" "$filename"
  done
  printf '%s\n' "$POPULATION_METADATA_URL"
  printf '%s\n' "$PEDIGREE_METADATA_URL"
}

build_vcf_path_list() {
  local chr filename
  for chr in "$@"; do
    filename="$(vcf_name_for_chr "$chr")"
    printf '%s/%s\n' "$OUT_DIR" "$filename"
  done
}

file_size() {
  local file="$1"
  if stat -f '%z' "$file" >/dev/null 2>&1; then
    stat -f '%z' "$file"
  else
    stat -c '%s' "$file"
  fi
}

remote_size() {
  local url="$1"
  curl --fail --location --silent --show-error --head "$url" \
    | awk 'tolower($1) == "content-length:" { size=$2 } END { gsub(/\r/, "", size); print size }'
}

write_pending_url_list() {
  local input_file="$1"
  local output_file="$2"
  local url filename local_file local_size remote_bytes

  while IFS= read -r url; do
    filename="${url##*/}"
    local_file="${OUT_DIR}/${filename}"

    if [[ -f "$local_file" ]] && command -v curl >/dev/null 2>&1; then
      local_size="$(file_size "$local_file")"
      remote_bytes="$(remote_size "$url")"

      if [[ -n "$remote_bytes" && "$local_size" == "$remote_bytes" ]]; then
        printf 'Skipping complete file: %s\n' "$filename"
        continue
      fi
    fi

    printf '%s\n' "$url" >> "$output_file"
  done < "$input_file"
}

download_with_aria2c() {
  local url_file="$1"
  aria2c \
    --continue=true \
    --auto-file-renaming=false \
    --allow-overwrite=true \
    --max-concurrent-downloads="$PARALLEL" \
    --split="$PARALLEL" \
    --max-connection-per-server="$PARALLEL" \
    --retry-wait=10 \
    --max-tries=0 \
    --dir="$OUT_DIR" \
    --input-file="$url_file"
}

download_one_with_curl() {
  local url="$1"
  local filename="${url##*/}"
  curl \
    --fail \
    --location \
    --continue-at - \
    --retry 10 \
    --retry-delay 10 \
    --retry-all-errors \
    --output "${OUT_DIR}/${filename}" \
    "$url"
}

download_with_curl() {
  local url_file="$1"
  export OUT_DIR
  export -f download_one_with_curl
  xargs -n 1 -P "$PARALLEL" bash -c 'download_one_with_curl "$0"' < "$url_file"
}

extract_vcf_samples() {
  local vcf="$1"
  if command -v bcftools >/dev/null 2>&1; then
    bcftools query --list-samples "$vcf"
  elif command -v zcat >/dev/null 2>&1; then
    zcat "$vcf" | awk '/^#CHROM/ { for (i = 10; i <= NF; i++) print $i; exit }'
  elif command -v gzcat >/dev/null 2>&1; then
    gzcat "$vcf" | awk '/^#CHROM/ { for (i = 10; i <= NF; i++) print $i; exit }'
  else
    printf 'Need bcftools, zcat, or gzcat to read VCF sample IDs from %s\n' "$vcf" >&2
    return 1
  fi
}

extract_metadata_ids() {
  local metadata_file="$1"
  awk '
    NR == 1 {
      sample_col = 1
      for (i = 1; i <= NF; i++) {
        field = tolower($i)
        if (field == "sampleid" || field == "sample" || field == "individual") {
          sample_col = i
          next
        }
      }
      if (tolower($0) ~ /sample/ || tolower($0) ~ /individual/) next
    }
    NF > 0 && $sample_col !~ /^#/ { print $sample_col }
  ' "$metadata_file" | sort -u
}

verify_labels() {
  local -a chromosomes=("$@")
  local population_file="${OUT_DIR}/${POPULATION_METADATA_URL##*/}"
  local pedigree_file="${OUT_DIR}/${PEDIGREE_METADATA_URL##*/}"

  if [[ ! -f "$population_file" ]]; then
    printf 'Missing population metadata: %s\n' "$population_file" >&2
    return 1
  fi
  if [[ ! -f "$pedigree_file" ]]; then
    printf 'Missing pedigree metadata: %s\n' "$pedigree_file" >&2
    return 1
  fi

  local pop_ids ped_ids vcf_ids missing_pop missing_ped extra_pop extra_ped
  pop_ids="$(mktemp)"
  ped_ids="$(mktemp)"
  vcf_ids="$(mktemp)"
  missing_pop="$(mktemp)"
  missing_ped="$(mktemp)"
  extra_pop="$(mktemp)"
  extra_ped="$(mktemp)"
  TMP_FILES+=("$pop_ids" "$ped_ids" "$vcf_ids" "$missing_pop" "$missing_ped" "$extra_pop" "$extra_ped")

  extract_metadata_ids "$population_file" > "$pop_ids"
  extract_metadata_ids "$pedigree_file" > "$ped_ids"

  local vcf
  while IFS= read -r vcf; do
    if [[ ! -f "$vcf" ]]; then
      printf 'Missing VCF for verification: %s\n' "$vcf" >&2
      return 1
    fi
    extract_vcf_samples "$vcf" | sort -u > "$vcf_ids"

    comm -23 "$vcf_ids" "$pop_ids" > "$missing_pop"
    comm -23 "$vcf_ids" "$ped_ids" > "$missing_ped"
    comm -13 "$vcf_ids" "$pop_ids" > "$extra_pop"
    comm -13 "$vcf_ids" "$ped_ids" > "$extra_ped"

    printf 'Verifying labels for %s\n' "${vcf##*/}"
    printf '  VCF samples: %s\n' "$(wc -l < "$vcf_ids" | tr -d ' ')"
    printf '  Population metadata IDs: %s\n' "$(wc -l < "$pop_ids" | tr -d ' ')"
    printf '  Pedigree metadata IDs: %s\n' "$(wc -l < "$ped_ids" | tr -d ' ')"

    if [[ -s "$missing_pop" || -s "$missing_ped" || -s "$extra_pop" || -s "$extra_ped" ]]; then
      printf '  Label check failed.\n' >&2
      if [[ -s "$missing_pop" ]]; then
        printf '  VCF IDs missing from population metadata, first 10:\n' >&2
        sed -n '1,10p' "$missing_pop" >&2
      fi
      if [[ -s "$missing_ped" ]]; then
        printf '  VCF IDs missing from pedigree metadata, first 10:\n' >&2
        sed -n '1,10p' "$missing_ped" >&2
      fi
      if [[ -s "$extra_pop" ]]; then
        printf '  Population metadata IDs not present in VCF, first 10:\n' >&2
        sed -n '1,10p' "$extra_pop" >&2
      fi
      if [[ -s "$extra_ped" ]]; then
        printf '  Pedigree metadata IDs not present in VCF, first 10:\n' >&2
        sed -n '1,10p' "$extra_ped" >&2
      fi
      return 1
    fi

    printf '  Label check passed: VCF samples match both metadata files exactly.\n'
  done < <(build_vcf_path_list "${chromosomes[@]}")
}

main() {
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
  fi

  local dry_run=0
  local verify=1
  local -a args=()
  local arg
  for arg in "$@"; do
    case "$arg" in
      --dry-run) dry_run=1 ;;
      --no-verify) verify=0 ;;
      *) args+=("$arg") ;;
    esac
  done

  local -a requested=("${args[@]:-1}")
  local -a chromosomes=()
  local chr normalized

  if [[ "${requested[0]}" == "all" ]]; then
    chromosomes=({1..22} X)
  else
    for chr in "${requested[@]}"; do
      normalized="$(normalize_chr "$chr")"
      chromosomes+=("$normalized")
    done
  fi

  mkdir -p "$OUT_DIR"

  local url_file
  url_file="$(mktemp)"
  TMP_FILES+=("$url_file")
  build_url_list "${chromosomes[@]}" > "$url_file"

  printf 'Downloading 1000 Genomes high coverage phased panel files to %s\n' "$OUT_DIR"
  printf 'Chromosomes: %s\n' "${chromosomes[*]}"

  if [[ "$dry_run" == "1" ]]; then
    cat "$url_file"
    exit 0
  fi

  local pending_url_file
  pending_url_file="$(mktemp)"
  TMP_FILES+=("$pending_url_file")
  write_pending_url_list "$url_file" "$pending_url_file"

  if [[ ! -s "$pending_url_file" ]]; then
    printf 'All requested files are already complete.\n'
  else
    if command -v aria2c >/dev/null 2>&1; then
      download_with_aria2c "$pending_url_file"
    elif command -v curl >/dev/null 2>&1; then
      download_with_curl "$pending_url_file"
    else
      printf 'Neither aria2c nor curl is available.\n' >&2
      exit 1
    fi
  fi

  if [[ "$verify" == "1" ]]; then
    verify_labels "${chromosomes[@]}"
  fi
}

main "$@"
