// BioVault popgen: PyPGx-rs PGx flow.
//
//   prepare_vcf  (BIOSYNTH_IMAGE, default ghcr.io/openmined/biosynth:0.1.32)
//       Passes VCF inputs through, or converts genotype TXT files to VCF.
//
//   pypgx_rs_pipeline  (ghcr.io/madhavajay/pypgx-rs:v0.26.0-rs.1)
//       Runs pypgx-rs run-ngs-pipeline against each VCF.
//
//   aggregate_pypgx
//       Adds country facets and aggregates unique genotypes by country/gene,
//       plus reclassification groups from results.long.tsv or alleles.zip.

nextflow.enable.dsl=2

if (!params.containsKey('pypgx_rs_image')) {
    params.pypgx_rs_image = null
}
if (!params.containsKey('pypgx_genes')) {
    params.pypgx_genes = null
}
if (!params.containsKey('pypgx_assembly')) {
    params.pypgx_assembly = 'GRCh38'
}

def BIOSYNTH_IMAGE = System.getenv('BIOSYNTH_IMAGE') ?: 'ghcr.io/openmined/biosynth:0.1.32'
def PYPGX_RS_IMAGE = params.pypgx_rs_image ?: 'ghcr.io/madhavajay/pypgx-rs:v0.26.0-rs.1'
def POPGEN_IMAGE = System.getenv('POPGEN_IMAGE') ?: 'ghcr.io/madhavajay/biovault-popgen:0.2.5-fast'

def normalizeFacet(String raw) {
    return (raw ?: '').trim()
}

def safeId(value) {
    return value.toString().replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

workflow USER {
    take:
        context
        participants

    main:
        def records = participants.flatMap { record ->
            def validation = record.validation ?: [status: 'ok', message: '']
            if (validation.status?.toString() != 'ok') {
                println "[bv] WARNING: skipping participant ${record.participant_id}: genotype file ${validation.status} - ${validation.message}"
                return []
            }
            def facets = record.facets ?: [:]
            def country = normalizeFacet((record.country ?: facets.country)?.toString())
            if (!country) {
                println "[bv] WARNING: participant ${record.participant_id}: missing country facet; using Unknown"
                country = 'Unknown'
            }
            return [tuple(record.participant_id.toString(), country, file(record.genotype_file))]
        }

        def checked = records
            .ifEmpty {
                throw new IllegalArgumentException("No valid participants with readable genotype files remained")
            }
        // Full expected roster — aggregate_pypgx diffs this against the
        // participants that actually produced a report to mark down any that
        // were dropped/killed (e.g. an OOM SIGKILL) before writing a status.
        def expected_manifest = checked
            .map { pid, country, gf -> "${pid}\t${country}" }
            .collectFile(name: 'expected_participants.tsv', newLine: true,
                         seed: "participant_id\tcountry")
        def vcfs = prepare_vcf(checked)
        def reports = pypgx_rs_pipeline(vcfs.vcf)
        def aggregate = aggregate_pypgx(
            reports.report_dir.collect(flat: false).ifEmpty([]),
            expected_manifest
        )

    emit:
        participant_results = aggregate.participant_results
        country_summary = aggregate.country_summary
        phase_reclassification_groups_aggregate = aggregate.phase_reclassification_groups_aggregate
        failures = aggregate.failures
        errors = aggregate.errors
        warnings = aggregate.warnings
        pipeline_log = aggregate.pipeline_log
}

process prepare_vcf {
    container BIOSYNTH_IMAGE
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    tag { participant_id }
    // Fail-soft: a single unconvertible file must not abort the cohort. A
    // dropped participant produces no VCF -> no report -> aggregate_pypgx marks
    // it down via the expected-vs-reported diff.
    errorStrategy 'ignore'
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_id), val(country), path(genotype_file)

    output:
        tuple val(participant_id), val(country), path("${prefix}/vcf/${prefix}.vcf*"), emit: vcf

    script:
    prefix = safeId(participant_id)
    def inputName = genotype_file.getName()
    """
    set -euo pipefail
    mkdir -p ${shellQuote("${prefix}/vcf")} ${shellQuote("${prefix}/logs")}
    input_name=${shellQuote(inputName)}
    lower="\$(printf '%s' "\${input_name}" | tr '[:upper:]' '[:lower:]')"
    case "\${lower}" in
        *.vcf)
            cp "\${input_name}" ${shellQuote("${prefix}/vcf/${prefix}.vcf")}
            ;;
        *.vcf.gz|*.vcf.bgz|*.bgz)
            cp "\${input_name}" ${shellQuote("${prefix}/vcf/${prefix}.vcf.gz")}
            ;;
        *)
            bvs genotype-to-vcf \\
                -i "\${input_name}" \\
                --sample ${shellQuote(participant_id)} \\
                --output ${shellQuote("${prefix}/vcf/${prefix}.vcf.gz")} \\
                --missing-log ${shellQuote("${prefix}/logs/${prefix}.vcf.log")} \\
                --gzip
            ;;
    esac
    """
}

process pypgx_rs_pipeline {
    container PYPGX_RS_IMAGE
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    tag { participant_id }
    // Fail-soft: one bad file must not kill the cohort. The script catches a
    // failed pypgx run, records a status row, and exits 0. 'ignore' is the
    // backstop for an uncatchable container kill — e.g. an OOM SIGKILL (exit
    // 137) where the script gets no chance to write status; aggregate_pypgx
    // then marks the participant down via the expected-vs-reported diff.
    errorStrategy 'ignore'
    maxForks 1
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_id), val(country), path(vcf_file)

    output:
        path "pypgx_${prefix}", emit: report_dir

    script:
    prefix = safeId(participant_id)
    def genesArg = params.pypgx_genes ? "--genes ${shellQuote(params.pypgx_genes)}" : ""
    """
    set -uo pipefail
    out="pypgx_${prefix}"
    mkdir -p "\${out}/raw"
    # Metadata up-front so a failed sample stays attributable.
    {
        echo "participant_id=${participant_id}"
        echo "country=${country}"
        echo "assembly=${params.pypgx_assembly}"
        echo "genes=${params.pypgx_genes ?: ''}"
    } > "\${out}/metadata.properties"

    run_log="\${out}/run.log"
    status=ok
    reason=""
    input_name=${shellQuote(vcf_file.getName())}
    # Guarded subshell: VCF normalisation + pypgx. pipefail is inherited.
    (
        set -e
        if ! command -v bgzip >/dev/null 2>&1 || ! command -v tabix >/dev/null 2>&1; then
            echo "pypgx-rs image must provide bgzip and tabix" >&2
            exit 1
        fi
        case "\${input_name}" in
            *.gz|*.bgz)
                gzip -dc "\${input_name}" > input.raw.vcf
                ;;
            *)
                cat "\${input_name}" > input.raw.vcf
                ;;
        esac
        awk '
            /^#/ { print > "input.header.vcf"; next }
            {
                chrom = \$1
                sub(/^chr/, "", chrom)
                rank = 0
                if (chrom ~ /^[0-9]+\$/) rank = chrom + 0
                else if (chrom == "X") rank = 23
                else if (chrom == "Y") rank = 24
                else if (chrom == "M" || chrom == "MT") rank = 25
                else next
                if (rank < 1 || rank > 25) next
                print rank "\\t" \$2 "\\t" \$0
            }
        ' input.raw.vcf \\
            | sort -t \$'\\t' -k1,1n -k2,2n \\
            | cut -f3- > input.records.vcf
        cat input.header.vcf input.records.vcf | bgzip -c > input.vcf.gz
        tabix -f -p vcf input.vcf.gz

        pypgx run-ngs-pipeline \\
            --vcf input.vcf.gz \\
            --assembly ${shellQuote(params.pypgx_assembly)} \\
            --output "\${out}/raw" \\
            ${genesArg}
    ) > "\${run_log}" 2>&1
    rc=\$?
    if [ "\${rc}" -ne 0 ]; then
        status=failed
        reason="\$(tail -n 5 "\${run_log}" 2>/dev/null | tr '\\t\\n' '  ' | sed 's/  */ /g' | cut -c1-300)"
        echo "[bv] WARNING: pypgx_rs_pipeline failed for ${shellQuote(participant_id)} (rc=\${rc}): \${reason}" >&2
    fi
    { printf 'participant_id\\tstatus\\treason\\n'; \\
      printf '%s\\t%s\\t%s\\n' ${shellQuote(participant_id)} "\${status}" "\${reason}"; } \\
      > "\${out}/status.tsv"
    exit 0
    """
}

process aggregate_pypgx {
    container POPGEN_IMAGE
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path report_dirs
        path expected_manifest

    output:
        path "pypgx_participant_results.tsv", emit: participant_results
        path "pypgx_country_summary.tsv", emit: country_summary
        path "phase_reclassification_groups_aggregate.csv", emit: phase_reclassification_groups_aggregate
        path "pypgx_failures.tsv", emit: failures
        path "errors.tsv", emit: errors
        path "warnings.tsv", emit: warnings
        path "pypgx_pipeline.log", emit: pipeline_log

    script:
    """
    set -euo pipefail

    python3 - <<'PY'
import csv
import hashlib
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

def clean(value):
    return (value or "").strip()

def norm_key(name):
    return re.sub(r"[^a-z0-9]+", "_", clean(name).lower()).strip("_")

def read_meta(report_dir):
    meta = {"participant_id": report_dir.name.removeprefix("pypgx_"), "country": "Unknown"}
    path = report_dir / "metadata.properties"
    if path.exists():
        for line in path.read_text(errors="replace").splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                meta[k.strip()] = v.strip()
    return meta

def detect_delimiter(path):
    text = path.read_text(errors="replace")[:4096]
    if "\\t" in text:
        return "\\t"
    return ","

def get(row, *names):
    lookup = {norm_key(k): v for k, v in row.items()}
    for name in names:
        v = lookup.get(norm_key(name))
        if v is not None and clean(v):
            return clean(v)
    return ""

def result_rows(report_dir):
    raw = report_dir / "raw"
    candidates = []
    top = raw / "results.tsv"
    if top.exists():
        candidates.append((None, top))
    candidates.extend((p.parent.name, p) for p in sorted(raw.glob("*/results.tsv")))

    rows = []
    for inferred_gene, path in candidates:
        try:
            with path.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter=detect_delimiter(path))
                for row in reader:
                    gene = get(row, "gene", "Gene", "symbol") or inferred_gene or path.parent.name
                    genotype = get(row, "genotype", "Genotype", "diplotype", "Diplotype", "haplotypes", "Haplotypes")
                    phenotype = get(row, "phenotype", "Phenotype", "phenotype_or_error", "activity", "Activity")
                    status = get(row, "status", "Status") or ("OK" if genotype else "ERROR")
                    error = get(row, "error", "Error", "message", "Message")
                    sample = get(row, "sample", "Sample", "sample_id", "Sample ID")
                    hap1 = get(row, "Haplotype1", "haplotype1")
                    hap2 = get(row, "Haplotype2", "haplotype2")
                    if not genotype and hap1 and hap2:
                        genotype = f"{hap1.rstrip(';')}/{hap2.rstrip(';')}"
                    if not genotype and not error:
                        continue
                    rows.append({
                        "sample": sample,
                        "gene": gene,
                        "status": status,
                        "genotype": genotype,
                        "phenotype": phenotype,
                        "error": error,
                        "source_file": str(path),
                    })
        except Exception as exc:
            rows.append({
                "sample": "",
                "gene": inferred_gene or path.parent.name,
                "status": "ERROR",
                "genotype": "",
                "phenotype": "",
                "error": f"parse_error: {exc}",
                "source_file": str(path),
            })
    return rows

def result_candidates(raw_dir, suffix):
    candidates = []
    top = raw_dir / suffix
    if top.exists():
        candidates.append((None, top))
    candidates.extend((p.parent.name, p) for p in sorted(raw_dir.glob(f"*/{suffix}")))
    return candidates

def phase_raw_dirs(report_dir):
    found = {}
    for mode in ("phased", "unphased"):
        candidates = [
            report_dir / mode / "raw",
            report_dir / "raw" / mode,
            report_dir / "raw" / mode / "raw",
        ]
        for raw_dir in candidates:
            if raw_dir.exists():
                found[mode] = raw_dir
                break
    if not found and (report_dir / "raw").exists():
        found["phased"] = report_dir / "raw"
    return found

def rows_by_gene(raw_dir):
    rows = defaultdict(list)
    for inferred_gene, path in result_candidates(raw_dir, "results.tsv"):
        try:
            with path.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter=detect_delimiter(path))
                for row in reader:
                    gene = get(row, "gene", "Gene", "symbol") or inferred_gene or path.parent.name
                    genotype = get(row, "genotype", "Genotype", "diplotype", "Diplotype", "haplotypes", "Haplotypes")
                    phenotype = get(row, "phenotype", "Phenotype", "phenotype_or_error", "activity", "Activity")
                    status = get(row, "status", "Status") or ("OK" if genotype else "ERROR")
                    error = get(row, "error", "Error", "message", "Message")
                    hap1 = get(row, "Haplotype1", "haplotype1")
                    hap2 = get(row, "Haplotype2", "haplotype2")
                    if not genotype and hap1 and hap2:
                        genotype = f"{hap1.rstrip(';')}/{hap2.rstrip(';')}"
                    if not gene or (not genotype and not phenotype and not error):
                        continue
                    rows[gene].append({
                        "gene": gene,
                        "status": status,
                        "genotype": genotype,
                        "phenotype": phenotype,
                        "error": error,
                    })
        except Exception as exc:
            gene = inferred_gene or path.parent.name
            rows[gene].append({
                "gene": gene,
                "status": "ERROR",
                "genotype": "",
                "phenotype": "",
                "error": f"parse_error: {exc}",
            })
    return rows

def first_gene_row(rows, gene):
    values = rows.get(gene) or []
    if not values:
        return {}
    return values[0]

def hap_label(value):
    raw = clean(value)
    if not raw:
        return ""
    low = raw.lower()
    if low in ("1", "h1", "hap1", "haplotype1", "haplotype 1"):
        return "Haplotype1"
    if low in ("2", "h2", "hap2", "haplotype2", "haplotype 2"):
        return "Haplotype2"
    m = re.search(r"([12])", low)
    if m:
        return f"Haplotype{m.group(1)}"
    return raw.replace(" ", "")

def sort_rank(value, fallback):
    raw = clean(value)
    if not raw:
        return (fallback, "")
    try:
        return (float(raw), raw)
    except ValueError:
        return (fallback, raw)

def variant_label(row):
    direct = get(
        row,
        "variant",
        "Variant",
        "variant_data",
        "VariantData",
        "marker",
        "Marker",
        "site",
        "Site",
    )
    rsid = get(row, "rsid", "RSID", "rs_id", "dbsnp", "DbSNP")
    chrom = get(row, "chromosome", "Chromosome", "chrom", "Chr", "chr")
    pos = get(row, "position", "Position", "pos", "grch38_position", "GRCh38Position", "start")
    ref = get(row, "ref", "Ref", "reference", "Reference", "grch38_ref", "GRCh38Allele")
    alt = get(row, "alt", "Alt", "alternate", "Alternate", "variant_allele", "VariantAllele")
    if rsid and chrom and pos and ref and alt:
        return f"{rsid}@{chrom}-{pos}-{ref}-{alt}"
    if direct and rsid and "@" not in direct and direct != rsid:
        return f"{rsid}@{direct}"
    return direct or rsid

def candidate_label(row, order):
    rank = get(row, "rank", "Rank", "candidate_rank", "CandidateRank", "priority", "Priority") or str(order)
    allele = get(
        row,
        "allele",
        "Allele",
        "candidate",
        "Candidate",
        "star_allele",
        "StarAllele",
        "candidate_allele",
        "CandidateAllele",
        "haplotype_allele",
        "HaplotypeAllele",
    )
    variant = variant_label(row)
    body = allele or variant
    if allele and variant:
        if variant == allele or variant.startswith(f"{allele}@"):
            body = variant
        elif variant not in allele:
            body = f"{allele}@{variant}"
    if not body:
        body = get(row, "genotype", "Genotype", "call", "Call")
    if not body:
        body = f"row{order}"
    return f"{rank}:{body}", sort_rank(rank, order)

def split_semicolon_list(value):
    return [part.strip() for part in clean(value).split(";") if part.strip()]

def parse_variant_data(value):
    lookup = {}
    for entry in split_semicolon_list(value):
        allele, sep, rest = entry.partition(":")
        if not sep:
            lookup[allele] = ""
            continue
        variant = rest.rsplit(":", 1)[0] if ":" in rest else rest
        lookup[allele] = variant
    return lookup

def allele_candidate(allele, variant, order):
    body = clean(allele)
    variant = clean(variant)
    if variant and variant != "default":
        if variant == body or variant.startswith(f"{body}@"):
            body = variant
        else:
            body = f"{body}@{variant}"
    elif variant == "default":
        body = f"{body}@default"
    return f"{order}:{body}"

def alleles_zip_signatures(raw_dir):
    signatures = {}
    for inferred_gene, path in result_candidates(raw_dir, "alleles.zip"):
        try:
            with zipfile.ZipFile(path) as archive:
                data_names = [name for name in archive.namelist() if name.endswith("/data.tsv") or name == "data.tsv"]
                if not data_names:
                    continue
                with archive.open(data_names[0]) as raw_handle:
                    text = raw_handle.read().decode("utf-8", errors="replace").splitlines()
            reader = csv.DictReader(text, delimiter="\\t")
            for row in reader:
                gene = get(row, "gene", "Gene", "symbol") or inferred_gene or path.parent.name
                variant_lookup = parse_variant_data(get(row, "VariantData", "variant_data"))
                hap_parts = []
                hap_bodies = []
                top_parts = []
                for hap_name in ("Haplotype1", "Haplotype2"):
                    alleles = split_semicolon_list(get(row, hap_name, hap_name.lower()))
                    body_candidates = [
                        allele_candidate(allele, variant_lookup.get(allele, ""), index)
                        for index, allele in enumerate(alleles, 1)
                    ]
                    body = "+".join(body_candidates)
                    hap_parts.append(f"{hap_name}={body}")
                    hap_bodies.append(body)
                    if body_candidates:
                        top_parts.append(f"{hap_name}={body_candidates[0]}")
                signatures[gene] = {
                    "top": "/".join(top_parts),
                    "full": "/".join(hap_parts),
                }
        except Exception as exc:
            gene = inferred_gene or path.parent.name
            signatures[gene] = {
                "top": f"parse_error:{exc}",
                "full": f"parse_error:{exc}",
            }
    return signatures

def long_signatures(raw_dir):
    by_gene = defaultdict(lambda: defaultdict(list))
    for inferred_gene, path in result_candidates(raw_dir, "results.long.tsv"):
        try:
            with path.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter=detect_delimiter(path))
                for order, row in enumerate(reader, 1):
                    gene = get(row, "gene", "Gene", "symbol") or inferred_gene or path.parent.name
                    hap = hap_label(get(row, "haplotype", "Haplotype", "phase", "Phase", "hap", "Hap"))
                    if not hap:
                        hap = "Haplotype1"
                    candidate, rank = candidate_label(row, order)
                    by_gene[gene][hap].append((rank, candidate))
        except Exception as exc:
            gene = inferred_gene or path.parent.name
            by_gene[gene]["parse_error"].append(((0, ""), f"parse_error:{exc}"))

    signatures = {}
    for gene, haps in by_gene.items():
        hap_parts = []
        hap_bodies = []
        top_parts = []
        for hap, candidates in sorted(haps.items()):
            ordered = [candidate for _rank, candidate in sorted(candidates, key=lambda item: item[0])]
            body = "+".join(ordered)
            hap_parts.append(f"{hap}={body}")
            hap_bodies.append(body)
            if ordered:
                top_parts.append(f"{hap}={ordered[0]}")
        signatures[gene] = {
            "top": "/".join(top_parts),
            "full": "/".join(hap_parts),
        }
    for gene, signature in alleles_zip_signatures(raw_dir).items():
        signatures.setdefault(gene, signature)
    return signatures

def chromosome_sort_key(chromosome):
    chrom = clean(chromosome)
    chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
    if chrom.isdigit():
        return (0, int(chrom), chrom)
    if chrom.upper() == "X":
        return (0, 23, chrom)
    if chrom.upper() == "Y":
        return (0, 24, chrom)
    if chrom.upper() in ("M", "MT"):
        return (0, 25, chrom)
    return (1, 0, chrom)

def split_gt(gt):
    sep = "|" if "|" in gt else "/"
    return sep, gt.split(sep)

def canonical_biallelic_gt(alleles, alt_index, sep):
    mapped = []
    for allele in alleles:
        if allele == ".":
            mapped.append(".")
        elif allele == "0":
            mapped.append("0")
        elif allele == str(alt_index):
            mapped.append("1")
        else:
            mapped.append("0")
    if all(value == "." for value in mapped):
        return "./."
    if len(mapped) == 2 and mapped[0] == "1" and mapped[1] == "1":
        return "1/1"
    return sep.join(mapped)

def raw_allele_tokens_from_vcf_lines(lines):
    tokens = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) < 10:
            continue
        chrom, pos, _id, ref, alt_text = parts[:5]
        format_keys = parts[8].split(":")
        sample_values = parts[-1].split(":")
        try:
            gt_index = format_keys.index("GT")
            gt = sample_values[gt_index]
        except (ValueError, IndexError):
            gt = sample_values[0] if sample_values else ""
        gt = clean(gt)
        if not gt:
            continue
        if gt in ("0|0", "0/0"):
            continue
        sep, gt_alleles = split_gt(gt)
        alts = [alt.strip() for alt in alt_text.split(",") if alt.strip()]
        if not alts:
            continue
        if all(allele == "." for allele in gt_alleles):
            alt_indices = range(1, len(alts) + 1)
        else:
            alt_indices = sorted(
                {int(allele) for allele in gt_alleles if allele.isdigit() and int(allele) > 0}
            )
        for alt_index in alt_indices:
            if alt_index > len(alts):
                continue
            alt = alts[alt_index - 1]
            token_gt = canonical_biallelic_gt(gt_alleles, alt_index, sep)
            try:
                pos_key = int(pos)
            except ValueError:
                pos_key = 0
            tokens.append((
                chromosome_sort_key(chrom),
                pos_key,
                ref,
                alt,
                f"{chrom}:{pos}:{ref}>{alt}={token_gt}",
            ))
    return [token for *_sort, token in sorted(tokens, key=lambda item: item[:4])]

def consolidated_raw_alleles(raw_dir):
    raw_alleles = {}
    for inferred_gene, path in result_candidates(raw_dir, "consolidated-variants.zip"):
        try:
            with zipfile.ZipFile(path) as archive:
                data_names = [name for name in archive.namelist() if name.endswith("/data.vcf") or name == "data.vcf"]
                if not data_names:
                    continue
                with archive.open(data_names[0]) as raw_handle:
                    lines = raw_handle.read().decode("utf-8", errors="replace").splitlines()
            gene = inferred_gene or path.parent.name
            raw_alleles[gene] = ";".join(raw_allele_tokens_from_vcf_lines(lines))
        except Exception as exc:
            gene = inferred_gene or path.parent.name
            raw_alleles[gene] = f"parse_error:{exc}"
    return raw_alleles

def stable_group_id(values):
    joined = "\\x1f".join(clean(v) for v in values)
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()[:12]

participant_rows = []
phase_rows = []
errors = []
warnings = []
failures = []
reported_pids = set()
for report_dir in sorted(Path(".").glob("pypgx_*")):
    if not report_dir.is_dir():
        continue
    meta = read_meta(report_dir)
    pid = meta.get("participant_id") or report_dir.name.removeprefix("pypgx_")
    country = meta.get("country") or "Unknown"
    if pid:
        reported_pids.add(pid)

    # Per-sample fail-soft status written by pypgx_rs_pipeline.
    sample_status, sample_reason = "ok", ""
    status_path = report_dir / "status.tsv"
    if status_path.exists():
        with open(status_path, newline="") as handle:
            srow = next(csv.DictReader(handle, delimiter="\\t"), {}) or {}
        sample_status = clean(srow.get("status")) or "ok"
        sample_reason = clean(srow.get("reason"))
    if sample_status != "ok":
        failures.append({"participant_id": pid, "country": country,
                         "status": sample_status,
                         "reason": sample_reason or "pypgx_rs_pipeline failed"})
        errors.append((pid, "ERROR", "pypgx_failed", sample_reason or "pypgx_rs_pipeline failed"))
        continue

    rows = result_rows(report_dir)
    if not rows:
        for mode, raw_dir in sorted(phase_raw_dirs(report_dir).items()):
            for gene, gene_rows in rows_by_gene(raw_dir).items():
                for gene_row in gene_rows:
                    rows.append({
                        "sample": "",
                        "gene": gene,
                        "status": gene_row["status"],
                        "genotype": gene_row["genotype"],
                        "phenotype": gene_row["phenotype"],
                        "error": gene_row["error"],
                        "source_file": f"{mode}:{raw_dir}",
                    })
            if rows:
                break
    if not rows:
        errors.append((pid, "ERROR", "no_results", f"No PyPGx results.tsv found under {report_dir}/raw"))
    for row in rows:
        gene = clean(row["gene"])
        status = clean(row["status"]) or "OK"
        genotype = clean(row["genotype"])
        phenotype = clean(row["phenotype"])
        error = clean(row["error"])
        if status.upper() != "OK" or error:
            errors.append((pid, "ERROR", "pypgx_gene_error", f"{gene}: {error or status}"))
        participant_rows.append({
            "participant_id": pid,
            "country": country,
            "gene": gene,
            "status": status,
            "genotype": genotype,
            "phenotype": phenotype,
            "error": error,
        })

    compact_rows = {}
    signatures = {}
    raw_alleles_by_gene = {}
    phase_genes = set()
    for _mode, raw_dir in sorted(phase_raw_dirs(report_dir).items()):
        compact_rows = rows_by_gene(raw_dir)
        signatures = long_signatures(raw_dir)
        raw_alleles_by_gene = consolidated_raw_alleles(raw_dir)
        phase_genes.update(compact_rows.keys())
        phase_genes.update(signatures.keys())
        phase_genes.update(raw_alleles_by_gene.keys())
        break
    for gene in sorted(phase_genes):
        compact = first_gene_row(compact_rows, gene)
        signature = signatures.get(gene, {})
        phase_rows.append({
            "participant_id": pid,
            "country": country,
            "gene": gene,
            "genotype": clean(compact.get("genotype")),
            "phenotype": clean(compact.get("phenotype")),
            "raw_alleles": clean(raw_alleles_by_gene.get(gene)),
        })

participant_rows.sort(key=lambda r: (r["participant_id"], r["gene"]))
phase_rows.sort(key=lambda r: (r["participant_id"], r["gene"]))

summary_counts = Counter()
summary_samples = set()
for row in participant_rows:
    if not row["genotype"]:
        continue
    key = (row["country"], row["gene"], row["genotype"], row["phenotype"])
    summary_counts[key] += 1
    summary_samples.add((row["country"], row["gene"], row["participant_id"]))
summary_denom = Counter((country, gene) for country, gene, _pid in summary_samples)
with open("pypgx_country_summary.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["country", "gene", "genotype", "phenotype", "count", "sample_count", "frequency"])
    for (country, gene, genotype, phenotype), count in sorted(summary_counts.items()):
        d = summary_denom[(country, gene)]
        writer.writerow([country, gene, genotype, phenotype, count, d, f"{(count / d) if d else 0:.6f}"])

phase_detail_fields = [
    "participant_id",
    "country",
    "gene",
    "genotype",
    "phenotype",
    "raw_alleles",
]
with open("pypgx_participant_results.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=phase_detail_fields, delimiter="\\t")
    writer.writeheader()
    writer.writerows(phase_rows)

phase_group_key_fields = [
    "country",
    "gene",
    "genotype",
    "phenotype",
    "raw_alleles",
]
grouped = {}
for row in phase_rows:
    key = tuple(row[field] for field in phase_group_key_fields)
    if key not in grouped:
        grouped[key] = {"row": row, "participants": set()}
    grouped[key]["participants"].add(row["participant_id"])
phase_summary_samples = set(
    (row["country"], row["gene"], row["participant_id"])
    for row in phase_rows
    if row["genotype"]
)
phase_summary_denom = Counter((country, gene) for country, gene, _pid in phase_summary_samples)

phase_aggregate_fields = [
    "country",
    "gene",
    "genotype",
    "phenotype",
    "count",
    "sample_count",
    "frequency",
    "group_id",
    "raw_alleles",
]
aggregate_rows = []
for key, payload in grouped.items():
    row = payload["row"]
    group_id = stable_group_id(key)
    count = len(payload["participants"])
    sample_count = phase_summary_denom[(row["country"], row["gene"])]
    aggregate_rows.append({
        "country": row["country"],
        "gene": row["gene"],
        "genotype": row["genotype"],
        "phenotype": row["phenotype"],
        "count": count,
        "sample_count": sample_count,
        "frequency": f"{(count / sample_count) if sample_count else 0:.6f}",
        "group_id": group_id,
        "raw_alleles": row["raw_alleles"],
    })
aggregate_rows.sort(key=lambda r: (r["country"], r["gene"], r["genotype"], r["phenotype"], r["group_id"]))

with open("phase_reclassification_groups_aggregate.csv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=phase_aggregate_fields)
    writer.writeheader()
    writer.writerows(aggregate_rows)

# Mark down any expected participant that never produced a report (dropped in
# prepare_vcf, or an uncatchable OOM SIGKILL in pypgx_rs_pipeline).
expected = {}
expected_path = Path("expected_participants.tsv")
if expected_path.exists():
    with open(expected_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\\t"):
            epid = clean(row.get("participant_id"))
            if epid:
                expected[epid] = clean(row.get("country"))
for epid, ecountry in expected.items():
    if epid not in reported_pids:
        failures.append({"participant_id": epid, "country": ecountry, "status": "failed",
                         "reason": "no PyPGx output produced (task dropped or killed)"})
        errors.append((epid, "ERROR", "no_output",
                       "no PyPGx output produced (task dropped or killed)"))

with open("pypgx_failures.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["participant_id", "country", "status", "reason"], delimiter="\\t")
    writer.writeheader()
    writer.writerows(sorted(failures, key=lambda r: (r["country"], r["participant_id"])))

with open("errors.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["participant_id", "severity", "code", "message"])
    writer.writerows(errors)

with open("warnings.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["participant_id", "severity", "code", "message"])
    writer.writerows(warnings)

with open("pypgx_pipeline.log", "w") as handle:
    handle.write(f"expected_participants={len(expected)}\\n")
    handle.write(f"succeeded_participants={len(set(r['participant_id'] for r in participant_rows))}\\n")
    handle.write(f"failed_participants={len(failures)}\\n")
    handle.write(f"participants={len(set(r['participant_id'] for r in participant_rows))}\\n")
    handle.write(f"participant_rows={len(participant_rows)}\\n")
    handle.write(f"phase_reclassification_rows={len(phase_rows)}\\n")
    handle.write(f"phase_reclassification_groups={len(aggregate_rows)}\\n")
    handle.write(f"errors={len(errors)}\\n")

if failures:
    print(f"[bv] PyPGx: {len(failures)} participant(s) marked down (see pypgx_failures.tsv); "
          f"{len(set(r['participant_id'] for r in participant_rows))} succeeded", flush=True)
PY
    """
}
