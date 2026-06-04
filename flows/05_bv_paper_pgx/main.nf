// BioVault popgen: PharmCAT PGx flow.
//
//   prepare_vcf  (BIOSYNTH_IMAGE, default ghcr.io/openmined/biosynth:0.1.31)
//       Passes VCF inputs through, or converts genotype TXT files to VCF:
//       bvs genotype-to-vcf -i <genotype.txt> --output <pid>.vcf.gz --gzip
//
//   pharmcat_pipeline  (pgkb/pharmcat)
//       Runs PharmCAT's Docker pipeline on each VCF and writes per-participant
//       TSV/JSON/HTML reports.
//
//   aggregate_pgx
//       Adds country facets to every PharmCAT gene row and writes grouped
//       summaries, including unresolved possible genotype combinations.

nextflow.enable.dsl=2

def BIOSYNTH_IMAGE = System.getenv('BIOSYNTH_IMAGE') ?: 'ghcr.io/openmined/biosynth:0.1.31'
def PHARMCAT_IMAGE = System.getenv('PHARMCAT_IMAGE') ?: 'pgkb/pharmcat'

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
                println "[bv] WARNING: skipping participant ${record.participant_id}: missing required country facet"
                return []
            }
            return [tuple(record.participant_id.toString(), country, file(record.genotype_file))]
        }

        def checked = records
            .ifEmpty {
                throw new IllegalArgumentException("No valid participants with readable genotype files and country facets remained")
            }

        def vcfs = prepare_vcf(checked)
        def reports = pharmcat_pipeline(vcfs.vcf)
        def aggregate = aggregate_pgx(
            reports.report_dir.collect(flat: false),
            vcfs.vcf_file.collect(flat: false)
        )

    emit:
        participant_results = aggregate.participant_results
        participant_possible_genotypes = aggregate.participant_possible_genotypes
        country_gene_genotype_counts = aggregate.country_gene_genotype_counts
        country_summary = aggregate.country_summary
        participant_manifest = aggregate.participant_manifest
        gene_country_burden = aggregate.gene_country_burden
        pgx_plots = aggregate.pgx_plots
        converted_vcfs = aggregate.converted_vcfs
        pharmcat_reports = aggregate.pharmcat_reports
        pharmcat_json = aggregate.pharmcat_json
        pharmcat_html = aggregate.pharmcat_html
        errors = aggregate.errors
        warnings = aggregate.warnings
        pipeline_log = aggregate.pipeline_log
}

process prepare_vcf {
    container BIOSYNTH_IMAGE
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    tag { participant_id }
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_id), val(country), path(genotype_file)

    output:
        tuple val(participant_id), val(country), path("${prefix}/vcf/${prefix}.vcf*"), emit: vcf
        path "${prefix}/vcf/${prefix}.vcf*", emit: vcf_file

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

process pharmcat_pipeline {
    container PHARMCAT_IMAGE
    stageInMode 'symlink'
    tag { participant_id }
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_id), val(country), path(vcf_file)

    output:
        path "pharmcat_${prefix}", emit: report_dir

    script:
    prefix = safeId(participant_id)
    """
    set -euo pipefail
    out="pharmcat_${prefix}"
    mkdir -p "\${out}"
    pharmcat_pipeline ${shellQuote(vcf_file.getName())} \\
        -o "\${out}" \\
        -bf ${shellQuote(prefix)} \\
        -reporterHtml \\
        -reporterJson \\
        -reporterCallsOnlyTsv
    printf 'participant_id\\tcountry\\tprefix\\n' > "\${out}/metadata.tsv"
    printf '%s\\t%s\\t%s\\n' ${shellQuote(participant_id)} ${shellQuote(country)} ${shellQuote(prefix)} >> "\${out}/metadata.tsv"
    """
}

process aggregate_pgx {
    container 'ghcr.io/madhavajay/biovault-popgen:0.2.2-fast'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path report_dirs
        path vcf_files

    output:
        path "pgx_participant_results.tsv", emit: participant_results
        path "pgx_participant_possible_genotypes.tsv", emit: participant_possible_genotypes
        path "pgx_country_gene_genotype_counts.tsv", emit: country_gene_genotype_counts
        path "pgx_country_summary.tsv", emit: country_summary
        path "pgx_participant_manifest.tsv", emit: participant_manifest
        path "pgx_gene_country_burden.tsv", emit: gene_country_burden
        path "pgx_plots/*", emit: pgx_plots
        path "vcfs/*.vcf*", emit: converted_vcfs
        path "pharmcat_reports/*.report.tsv", emit: pharmcat_reports
        path "pharmcat_reports/*.report.json", emit: pharmcat_json, optional: true
        path "pharmcat_reports/*.report.html", emit: pharmcat_html, optional: true
        path "errors.tsv", emit: errors
        path "warnings.tsv", emit: warnings
        path "pgx_pipeline.log", emit: pipeline_log

    script:
    """
    set -euo pipefail
    mkdir -p vcfs pharmcat_reports
    cp *.vcf* vcfs/ 2>/dev/null || true

    python3 - <<'PY'
import csv
import glob
import json
from collections import Counter, defaultdict
from pathlib import Path

def clean(value):
    return (value or "").strip()

def read_report(path):
    with open(path, newline="") as handle:
        first = handle.readline()
        header = handle.readline()
        if not header:
            return []
        rows = csv.DictReader(handle, fieldnames=header.rstrip("\\n").split("\\t"), delimiter="\\t")
        return list(rows)

def format_possible_from_json(report_json):
    if not report_json.exists():
        return {}
    try:
        data = json.loads(report_json.read_text())
    except Exception:
        return {}
    out = {}
    genes = data.get("genes") or {}
    if not isinstance(genes, dict):
        return out
    for gene, payload in genes.items():
        if not isinstance(payload, dict):
            continue
        pieces = []
        seen = set()
        for dip in payload.get("sourceDiplotypes") or []:
            if not isinstance(dip, dict):
                continue
            label = clean(dip.get("label"))
            if not label:
                continue
            score = dip.get("matchScore")
            if score is None or score == "":
                piece = label
            else:
                piece = f"{label} ({score})"
            if piece not in seen:
                seen.add(piece)
                pieces.append(piece)
        if pieces:
            out[gene] = ", ".join(pieces)
    return out

def rows_from_json(report_json):
    if not report_json.exists():
        return []
    try:
        data = json.loads(report_json.read_text())
    except Exception:
        return []
    genes = data.get("genes") or {}
    if not isinstance(genes, dict):
        return []
    rows = []
    for gene, payload in sorted(genes.items()):
        if not isinstance(payload, dict):
            continue
        labels = []
        scores = []
        phenotypes = []
        for dip in payload.get("sourceDiplotypes") or []:
            if not isinstance(dip, dict):
                continue
            label = clean(dip.get("label"))
            if label:
                labels.append(label)
            score = dip.get("matchScore")
            if score is not None and score != "":
                scores.append(str(score))
            for phenotype in dip.get("phenotypes") or []:
                phenotype = clean(phenotype)
                if phenotype and phenotype not in phenotypes:
                    phenotypes.append(phenotype)
        rows.append({
            "Gene": gene,
            "Source Diplotype": " OR ".join(labels),
            "Phenotype": " / ".join(phenotypes),
            "Activity Score": "",
            "Recommendation Lookup Diplotype": "",
            "Recommendation Lookup Phenotype": "",
            "Recommendation Lookup Activity Score": "",
            "Missing positions": "",
            "Outside Call": "",
            "Match Score": " / ".join(scores),
        })
    return rows

def format_possible_from_tsv(row):
    source = clean(row.get("Source Diplotype"))
    score = clean(row.get("Match Score"))
    if not source:
        return "no call"
    parts = [p.strip() for p in source.split(" OR ") if p.strip()]
    scores = [s.strip() for s in score.split(" / ") if s.strip()]
    if not parts:
        return source
    if len(scores) == len(parts):
        return ", ".join(f"{part} ({scores[idx]})" for idx, part in enumerate(parts))
    if len(scores) == 1:
        return ", ".join(f"{part} ({scores[0]})" for part in parts)
    return ", ".join(parts)

def is_reference_like(value):
    text = clean(value).lower()
    if not text:
        return False
    if text in {"no call", "none", "na", "n/a"}:
        return True
    compact = " ".join(text.split())
    if compact in {"reference/reference", "*1/*1", "b (reference)/b (reference)"}:
        return True
    if compact.startswith("rs") and " reference " in compact and "/" in compact and " variant " not in compact:
        return True
    return False

participant_fields = [
    "participant_id",
    "country",
    "gene",
    "source_diplotype",
    "phenotype",
    "activity_score",
    "recommendation_lookup_diplotype",
    "recommendation_lookup_phenotype",
    "recommendation_lookup_activity_score",
    "missing_positions",
    "outside_call",
    "match_score",
]
summary_fields_country = [
    "country",
    "gene",
    "phenotype",
    "recommendation_lookup_phenotype",
    "source_diplotype",
    "count",
    "sample_count",
]

participant_rows = []
possible_rows = []
manifest_rows = []
country_counts = Counter()
country_gene_genotype_counts = Counter()
country_gene_genotype_samples = defaultdict(set)
burden_country_gene_samples = defaultdict(set)
burden_country_gene_nonref = Counter()
country_samples = defaultdict(set)
errors = []
warnings = []

for report_dir in sorted(Path(".").glob("pharmcat_*")):
    if not report_dir.is_dir():
        continue
    if report_dir.name == "pharmcat_reports":
        continue
    metadata_path = report_dir / "metadata.tsv"
    if not metadata_path.exists():
        errors.append(["", str(report_dir), "ERROR", "MISSING_METADATA", "metadata.tsv missing"])
        continue
    with open(metadata_path, newline="") as handle:
        meta = next(csv.DictReader(handle, delimiter="\\t"))
    pid = clean(meta.get("participant_id"))
    country = clean(meta.get("country"))
    prefix = clean(meta.get("prefix")) or pid
    manifest_rows.append({
        "participant_id": pid,
        "country": country,
        "prefix": prefix,
    })

    for artifact in report_dir.glob("*.report.*"):
        target = Path("pharmcat_reports") / f"{prefix}{''.join(artifact.suffixes)}"
        target.write_bytes(artifact.read_bytes())

    report_paths = sorted(report_dir.glob("*.report.tsv"))
    if not report_paths:
        errors.append([pid, str(report_dir), "ERROR", "MISSING_REPORT_TSV", "No PharmCAT report TSV produced"])
        continue
    if len(report_paths) > 1:
        warnings.append([pid, str(report_dir), "WARNING", "MULTIPLE_REPORT_TSV", f"Using {report_paths[0].name}"])
    report_path = report_paths[0]
    possible_by_gene = {}
    json_rows = []
    report_json_paths = sorted(report_dir.glob("*.report.json"))
    if report_json_paths:
        possible_by_gene = format_possible_from_json(report_json_paths[0])
        json_rows = rows_from_json(report_json_paths[0])

    report_rows = read_report(report_path)
    if not report_rows and json_rows:
        report_rows = json_rows

    for row in report_rows:
        gene = clean(row.get("Gene"))
        if not gene:
            continue
        source_diplotype = clean(row.get("Source Diplotype"))
        phenotype = clean(row.get("Phenotype"))
        recommendation_phenotype = clean(row.get("Recommendation Lookup Phenotype"))
        possible_genotypes = possible_by_gene.get(gene) or format_possible_from_tsv(row)
        participant_rows.append({
            "participant_id": pid,
            "country": country,
            "gene": gene,
            "source_diplotype": source_diplotype,
            "phenotype": phenotype,
            "activity_score": clean(row.get("Activity Score")),
            "recommendation_lookup_diplotype": clean(row.get("Recommendation Lookup Diplotype")),
            "recommendation_lookup_phenotype": recommendation_phenotype,
            "recommendation_lookup_activity_score": clean(row.get("Recommendation Lookup Activity Score")),
            "missing_positions": clean(row.get("Missing positions")),
            "outside_call": clean(row.get("Outside Call")),
            "match_score": clean(row.get("Match Score")),
        })
        possible_rows.append({
            "participant_id": pid,
            "country": country,
            "gene": gene,
            "possible_genotypes": possible_genotypes,
        })
        c_key = (country, gene, phenotype, recommendation_phenotype, source_diplotype)
        possible_key = (country, gene, possible_genotypes)
        country_counts[c_key] += 1
        country_gene_genotype_counts[possible_key] += 1
        country_samples[c_key].add(pid)
        country_gene_genotype_samples[possible_key].add(pid)
        burden_cg_key = (country, gene)
        burden_country_gene_samples[burden_cg_key].add(pid)
        if not is_reference_like(source_diplotype):
            burden_country_gene_nonref[burden_cg_key] += 1

with open("pgx_participant_manifest.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["participant_id", "country", "prefix"], delimiter="\\t")
    writer.writeheader()
    writer.writerows(sorted(manifest_rows, key=lambda r: (r["country"], r["participant_id"])))

with open("pgx_participant_results.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=participant_fields, delimiter="\\t")
    writer.writeheader()
    writer.writerows(sorted(participant_rows, key=lambda r: (r["country"], r["participant_id"], r["gene"])))

with open("pgx_participant_possible_genotypes.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["participant_id", "country", "gene", "possible_genotypes"], delimiter="\\t")
    writer.writeheader()
    writer.writerows(sorted(possible_rows, key=lambda r: (r["country"], r["participant_id"], r["gene"], r["possible_genotypes"])))

with open("pgx_country_gene_genotype_counts.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["country", "gene", "possible_genotypes", "count", "sample_count"], delimiter="\\t")
    writer.writeheader()
    for key, count in sorted(country_gene_genotype_counts.items()):
        country, gene, possible_genotypes = key
        writer.writerow({
            "country": country,
            "gene": gene,
            "possible_genotypes": possible_genotypes,
            "count": count,
            "sample_count": len(country_gene_genotype_samples[key]),
        })

with open("pgx_country_summary.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_fields_country, delimiter="\\t")
    writer.writeheader()
    for key, count in sorted(country_counts.items()):
        country, gene, phenotype, rec_phenotype, source_diplotype = key
        writer.writerow({
            "country": country,
            "gene": gene,
            "phenotype": phenotype,
            "recommendation_lookup_phenotype": rec_phenotype,
            "source_diplotype": source_diplotype,
            "count": count,
            "sample_count": len(country_samples[key]),
        })

with open("pgx_gene_country_burden.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["country", "gene", "total_calls", "non_reference_calls", "non_reference_rate"], delimiter="\\t")
    writer.writeheader()
    for key in sorted(burden_country_gene_samples):
        total = len(burden_country_gene_samples[key])
        nonref = burden_country_gene_nonref[key]
        writer.writerow({
            "country": key[0],
            "gene": key[1],
            "total_calls": total,
            "non_reference_calls": nonref,
            "non_reference_rate": f"{(nonref / total) if total else 0:.6f}",
        })

with open("errors.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["participant_id", "file", "severity", "code", "message"])
    writer.writerows(errors)

with open("warnings.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["participant_id", "file", "severity", "code", "message"])
    writer.writerows(warnings)

with open("pgx_pipeline.log", "w") as handle:
    handle.write(f"participants={len(manifest_rows)}\\n")
    handle.write(f"participant_gene_rows={len(participant_rows)}\\n")
    handle.write(f"participant_possible_genotype_rows={len(possible_rows)}\\n")
    handle.write(f"country_gene_genotype_count_rows={len(country_gene_genotype_counts)}\\n")
    handle.write(f"country_summary_rows={len(country_counts)}\\n")
    handle.write(f"errors={len(errors)}\\n")
    handle.write(f"warnings={len(warnings)}\\n")
PY

    mkdir -p pgx_plots
    cp pgx_gene_country_burden.tsv pgx_plots/pgx_gene_country_burden.tsv
    """
}
