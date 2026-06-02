#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LOCUS_MAP="${LOCUS_MAP:-${SCRIPT_DIR}/locus_map.tsv}"
VCF_DIR="${VCF_DIR:-${REPO_ROOT}/data/1kgp_high_coverage}"
OUT_DIR="${OUT_DIR:-${VCF_DIR}/filtered}"
THREADS="${THREADS:-4}"
JOBS="${JOBS:-1}"

usage() {
  cat <<'USAGE'
Usage:
  tools/filter_1kgp_by_locus_map.sh [--dry-run] [--force] [--jobs N] [chr ...]

Examples:
  tools/filter_1kgp_by_locus_map.sh --dry-run 1
  tools/filter_1kgp_by_locus_map.sh 1
  tools/filter_1kgp_by_locus_map.sh 1 2 X
  tools/filter_1kgp_by_locus_map.sh all

Environment:
  LOCUS_MAP=/path/locus_map.tsv     Default: tools/locus_map.tsv
  VCF_DIR=/path/1kgp_high_coverage  Default: data/1kgp_high_coverage
  OUT_DIR=/path/output              Default: $VCF_DIR/filtered
  THREADS=4                         Compression threads for bcftools
  JOBS=1                            Chromosomes to process concurrently

The script filters 1000 Genomes high-coverage VCFs to loci in locus_map.tsv,
keeping only biallelic SNP records. The locus file must have header columns
compatible with chrom/chromosome/chr and pos/position; column order does not
matter.
It skips files with active .aria2 download markers and requires each input
VCF to have a current .tbi index, rebuilding missing or stale indexes unless
--dry-run is set.
USAGE
}

normalize_chr() {
  local chr="$1"
  chr="${chr#chr}"
  case "$chr" in
    [1-9]|1[0-9]|2[0-2]|X|Y|MT) printf '%s\n' "$chr" ;;
    M) printf 'MT\n' ;;
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

contig_for_chr() {
  local vcf="$1"
  local chr="$2"
  local id

  if bcftools view -h "$vcf" | grep -q "^##contig=<ID=chr${chr}"; then
    printf 'chr%s\n' "$chr"
    return 0
  fi
  if bcftools view -h "$vcf" | grep -q "^##contig=<ID=${chr}"; then
    printf '%s\n' "$chr"
    return 0
  fi

  id="$(bcftools view -h "$vcf" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' | head -n 1)"
  if [[ -n "$id" ]]; then
    printf '%s\n' "$id"
    return 0
  fi

  printf '%s\n' "$chr"
}

write_regions_for_chr() {
  local chr="$1"
  local contig="$2"
  local out="$3"

  awk -F '\t' -v chr="$chr" -v contig="$contig" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        name = tolower($i)
        if (name == "chrom" || name == "chromosome" || name == "chr") chrom_col = i
        if (name == "pos" || name == "position") pos_col = i
      }
      if (!chrom_col || !pos_col) {
        printf "Expected chrom/pos columns in %s\n", FILENAME > "/dev/stderr"
        exit 2
      }
      next
    }
    {
      row_chr = $chrom_col
      sub(/^chr/, "", row_chr)
      if (row_chr == "M") row_chr = "MT"
    }
    row_chr == chr && $pos_col ~ /^[0-9]+$/ && $pos_col > 0 {
      print contig "\t" $pos_col "\t" $pos_col
    }
  ' "$LOCUS_MAP" | sort -k1,1 -k2,2n -u > "$out"
}

ensure_vcf_index() {
  local vcf="$1"
  local dry_run="$2"
  local index="${vcf}.tbi"

  if [[ -f "$index" && ! "$vcf" -nt "$index" ]]; then
    return 0
  fi

  if [[ "$dry_run" == "1" ]]; then
    if [[ -f "$index" ]]; then
      printf 'Dry run: input index is stale, would run: tabix -f -p vcf %s\n' "$vcf"
    else
      printf 'Dry run: input index is missing, would run: tabix -f -p vcf %s\n' "$vcf"
    fi
    return 0
  fi

  if [[ -f "$index" ]]; then
    printf 'Input index is stale, rebuilding: %s\n' "$index"
  else
    printf 'Missing input index, building: %s\n' "$index"
  fi
  tabix -f -p vcf "$vcf"
}

filter_chr() {
  local chr="$1"
  local dry_run="$2"
  local force="$3"
  local vcf filename out_vcf regions contig n_regions

  filename="$(vcf_name_for_chr "$chr")"
  vcf="${VCF_DIR}/${filename}"
  out_vcf="${OUT_DIR}/${filename}"

  if [[ ! -f "$vcf" ]]; then
    printf 'Missing input VCF: %s\n' "$vcf" >&2
    return 1
  fi
  if [[ -f "${vcf}.aria2" ]]; then
    printf 'Input still has active aria2 marker, skipping incomplete download: %s.aria2\n' "$vcf" >&2
    return 1
  fi
  if [[ "$force" != "1" && -f "$out_vcf" ]]; then
    printf 'Output exists, use --force to replace: %s\n' "$out_vcf" >&2
    return 1
  fi

  ensure_vcf_index "$vcf" "$dry_run"

  mkdir -p "$OUT_DIR"
  regions="$(mktemp)"

  contig="$(contig_for_chr "$vcf" "$chr")"
  write_regions_for_chr "$chr" "$contig" "$regions"
  n_regions="$(wc -l < "$regions" | tr -d ' ')"
  if [[ "$n_regions" == "0" ]]; then
    printf 'No loci found in %s for chromosome %s\n' "$LOCUS_MAP" "$chr" >&2
    return 1
  fi

  printf 'Chromosome %s: %s target loci -> %s\n' "$chr" "$n_regions" "$out_vcf"

  if [[ "$dry_run" == "1" ]]; then
    printf 'Dry run: bcftools view --min-alleles 2 --max-alleles 2 --types snps --regions-file %s --targets-file %s --output-type z --threads %s --output %s %s\n' "$regions" "$regions" "$THREADS" "$out_vcf" "$vcf"
    rm -f "$regions"
    return 0
  fi

  bcftools view \
    --min-alleles 2 \
    --max-alleles 2 \
    --types snps \
    --regions-file "$regions" \
    --targets-file "$regions" \
    --output-type z \
    --threads "$THREADS" \
    --output "$out_vcf" \
    "$vcf"
  tabix -f -p vcf "$out_vcf"
  rm -f "$regions"
}

wait_for_oldest_job() {
  local status=0
  local pid="${active_pids[0]}"

  if ! wait "$pid"; then
    status=1
  fi
  active_pids=("${active_pids[@]:1}")
  return "$status"
}

main() {
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
  fi

  local dry_run=0
  local force=0
  local -a args=()
  local arg chr normalized next_arg

  while [[ "$#" -gt 0 ]]; do
    arg="$1"
    case "$arg" in
      --dry-run) dry_run=1 ;;
      --force) force=1 ;;
      --jobs)
        shift
        if [[ "$#" -eq 0 ]]; then
          printf '--jobs requires a positive integer.\n' >&2
          exit 1
        fi
        JOBS="$1"
        ;;
      --jobs=*)
        JOBS="${arg#--jobs=}"
        ;;
      *) args+=("$arg") ;;
    esac
    shift
  done

  if [[ ! "$JOBS" =~ ^[1-9][0-9]*$ ]]; then
    printf 'JOBS/--jobs must be a positive integer: %s\n' "$JOBS" >&2
    exit 1
  fi

  if [[ ! -f "$LOCUS_MAP" ]]; then
    printf 'Missing locus map: %s\n' "$LOCUS_MAP" >&2
    exit 1
  fi
  if ! command -v bcftools >/dev/null 2>&1; then
    printf 'bcftools is required but was not found on PATH.\n' >&2
    exit 1
  fi
  if ! command -v tabix >/dev/null 2>&1; then
    printf 'tabix is required but was not found on PATH.\n' >&2
    exit 1
  fi

  local -a requested=("${args[@]:-1}")
  local -a chromosomes=()
  if [[ "${requested[0]}" == "all" ]]; then
    chromosomes=({1..22} X)
  else
    for chr in "${requested[@]}"; do
      normalized="$(normalize_chr "$chr")"
      chromosomes+=("$normalized")
    done
  fi

  printf 'Locus map: %s\n' "$LOCUS_MAP"
  printf 'Input VCF dir: %s\n' "$VCF_DIR"
  printf 'Output dir: %s\n' "$OUT_DIR"
  printf 'Chromosomes: %s\n' "${chromosomes[*]}"
  printf 'Parallel jobs: %s chromosome(s), %s bcftools thread(s) each\n' "$JOBS" "$THREADS"

  if [[ "$JOBS" == "1" ]]; then
    for chr in "${chromosomes[@]}"; do
      filter_chr "$chr" "$dry_run" "$force"
    done
  else
    local -a active_pids=()
    local failed=0
    for chr in "${chromosomes[@]}"; do
      filter_chr "$chr" "$dry_run" "$force" &
      active_pids+=("$!")
      if [[ "${#active_pids[@]}" -ge "$JOBS" ]]; then
        if ! wait_for_oldest_job; then
          failed=1
        fi
      fi
    done
    while [[ "${#active_pids[@]}" -gt 0 ]]; do
      if ! wait_for_oldest_job; then
        failed=1
      fi
    done
    if [[ "$failed" != "0" ]]; then
      exit "$failed"
    fi
  fi
}

main "$@"
