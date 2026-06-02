// BioVault popgen: PharmCAT PGx flow.
//
//   genotype_to_vcf  (BIOSYNTH_IMAGE, default ghcr.io/openmined/biosynth:0.1.30)
//       Converts each selected genotype TXT to a compressed VCF:
//       bvs genotype-to-vcf -i <genotype.txt> --output <pid>.vcf.gz --gzip
//
//   pharmcat_pipeline  (pgkb/pharmcat)
//       Runs PharmCAT's Docker pipeline on each VCF and writes per-participant
//       TSV/JSON/HTML reports.
//
//   aggregate_pgx
//       Adds country/sex facets to every PharmCAT gene row and writes grouped
//       summaries. The country-only table is the sum across sex.

nextflow.enable.dsl=2

def BIOSYNTH_IMAGE = System.getenv('BIOSYNTH_IMAGE') ?: 'ghcr.io/openmined/biosynth:0.1.30'
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
            def sex = normalizeFacet((record.sex ?: facets.sex)?.toString())
            if (!country) {
                println "[bv] WARNING: skipping participant ${record.participant_id}: missing required country facet"
                return []
            }
            if (!sex) {
                println "[bv] WARNING: skipping participant ${record.participant_id}: missing required sex facet"
                return []
            }
            return [tuple(record.participant_id.toString(), country, sex, file(record.genotype_file))]
        }

        def checked = records
            .collect(flat: false)
            .map { items ->
                if (items.isEmpty()) {
                    throw new IllegalArgumentException("No valid participants with readable genotype files and country/sex facets remained")
                }
                println "[bv] pgx participants: ${items.size()}"
                tuple(items)
            }
            .flatMap { it[0] }

        def vcfs = genotype_to_vcf(checked)
        def reports = pharmcat_pipeline(vcfs.vcf)
        def plotScript = file("${projectDir}/../../04_population_level/pgx/plot_pgx_accumulation.py")
        def aggregate = aggregate_pgx(
            reports.report_dir.collect(flat: false),
            vcfs.vcf_file.collect(flat: false),
            plotScript
        )

    emit:
        participant_results = aggregate.participant_results
        country_sex_summary = aggregate.country_sex_summary
        country_summary = aggregate.country_summary
        participant_manifest = aggregate.participant_manifest
        gene_country_burden = aggregate.gene_country_burden
        gene_country_sex_burden = aggregate.gene_country_sex_burden
        pgx_plots = aggregate.pgx_plots
        converted_vcfs = aggregate.converted_vcfs
        pharmcat_reports = aggregate.pharmcat_reports
        pharmcat_json = aggregate.pharmcat_json
        pharmcat_html = aggregate.pharmcat_html
        errors = aggregate.errors
        warnings = aggregate.warnings
        pipeline_log = aggregate.pipeline_log
}

process genotype_to_vcf {
    container BIOSYNTH_IMAGE
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    tag { participant_id }
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_id), val(country), val(sex), path(genotype_file)

    output:
        tuple val(participant_id), val(country), val(sex), path("${prefix}/${prefix}.vcf.gz"), emit: vcf
        path "${prefix}/${prefix}.vcf.gz", emit: vcf_file

    script:
    prefix = safeId(participant_id)
    def inputName = genotype_file.getName()
    """
    set -euo pipefail
    mkdir -p ${shellQuote(prefix)}
    bvs genotype-to-vcf \\
        -i ${shellQuote(inputName)} \\
        --output ${shellQuote("${prefix}/${prefix}.vcf.gz")} \\
        --gzip
    """
}

process pharmcat_pipeline {
    container PHARMCAT_IMAGE
    stageInMode 'symlink'
    tag { participant_id }
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_id), val(country), val(sex), path(vcf_file)

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
    printf 'participant_id\\tcountry\\tsex\\tprefix\\n' > "\${out}/metadata.tsv"
    printf '%s\\t%s\\t%s\\t%s\\n' ${shellQuote(participant_id)} ${shellQuote(country)} ${shellQuote(sex)} ${shellQuote(prefix)} >> "\${out}/metadata.tsv"
    """
}

process aggregate_pgx {
    container 'ghcr.io/madhavajay/biovault-popgen:0.2.0-fast'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path report_dirs
        path vcf_files
        path plot_script

    output:
        path "pgx_participant_results.tsv", emit: participant_results
        path "pgx_country_sex_summary.tsv", emit: country_sex_summary
        path "pgx_country_summary.tsv", emit: country_summary
        path "pgx_participant_manifest.tsv", emit: participant_manifest
        path "pgx_gene_country_burden.tsv", emit: gene_country_burden
        path "pgx_gene_country_sex_burden.tsv", emit: gene_country_sex_burden
        path "pgx_plots/*", emit: pgx_plots
        path "vcfs/*.vcf.gz", emit: converted_vcfs
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
    cp *.vcf.gz vcfs/ 2>/dev/null || true

    python3 - <<'PY'
import csv
import glob
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

participant_fields = [
    "participant_id",
    "country",
    "sex",
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
summary_fields_country_sex = [
    "country",
    "sex",
    "gene",
    "phenotype",
    "recommendation_lookup_phenotype",
    "source_diplotype",
    "count",
    "sample_count",
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
manifest_rows = []
country_sex_counts = Counter()
country_counts = Counter()
country_sex_samples = defaultdict(set)
country_samples = defaultdict(set)
errors = []
warnings = []

for report_dir in sorted(Path(".").glob("pharmcat_*")):
    if not report_dir.is_dir():
        continue
    metadata_path = report_dir / "metadata.tsv"
    if not metadata_path.exists():
        errors.append(["", str(report_dir), "ERROR", "MISSING_METADATA", "metadata.tsv missing"])
        continue
    with open(metadata_path, newline="") as handle:
        meta = next(csv.DictReader(handle, delimiter="\\t"))
    pid = clean(meta.get("participant_id"))
    country = clean(meta.get("country"))
    sex = clean(meta.get("sex"))
    prefix = clean(meta.get("prefix")) or pid
    manifest_rows.append({
        "participant_id": pid,
        "country": country,
        "sex": sex,
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

    for row in read_report(report_path):
        gene = clean(row.get("Gene"))
        if not gene:
            continue
        source_diplotype = clean(row.get("Source Diplotype"))
        phenotype = clean(row.get("Phenotype"))
        recommendation_phenotype = clean(row.get("Recommendation Lookup Phenotype"))
        participant_rows.append({
            "participant_id": pid,
            "country": country,
            "sex": sex,
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
        cs_key = (country, sex, gene, phenotype, recommendation_phenotype, source_diplotype)
        c_key = (country, gene, phenotype, recommendation_phenotype, source_diplotype)
        country_sex_counts[cs_key] += 1
        country_counts[c_key] += 1
        country_sex_samples[cs_key].add(pid)
        country_samples[c_key].add(pid)

with open("pgx_participant_manifest.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["participant_id", "country", "sex", "prefix"], delimiter="\\t")
    writer.writeheader()
    writer.writerows(sorted(manifest_rows, key=lambda r: (r["country"], r["sex"], r["participant_id"])))

with open("pgx_participant_results.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=participant_fields, delimiter="\\t")
    writer.writeheader()
    writer.writerows(sorted(participant_rows, key=lambda r: (r["country"], r["sex"], r["participant_id"], r["gene"])))

with open("pgx_country_sex_summary.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_fields_country_sex, delimiter="\\t")
    writer.writeheader()
    for key, count in sorted(country_sex_counts.items()):
        country, sex, gene, phenotype, rec_phenotype, source_diplotype = key
        writer.writerow({
            "country": country,
            "sex": sex,
            "gene": gene,
            "phenotype": phenotype,
            "recommendation_lookup_phenotype": rec_phenotype,
            "source_diplotype": source_diplotype,
            "count": count,
            "sample_count": len(country_sex_samples[key]),
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
    handle.write(f"country_sex_summary_rows={len(country_sex_counts)}\\n")
    handle.write(f"country_summary_rows={len(country_counts)}\\n")
    handle.write(f"errors={len(errors)}\\n")
    handle.write(f"warnings={len(warnings)}\\n")
PY

    mkdir -p pgx_plots
    python3 ${shellQuote(plot_script.getName())} \\
        --participant-results pgx_participant_results.tsv \\
        --out-dir pgx_plots
    cp pgx_plots/pgx_gene_country_burden.tsv pgx_gene_country_burden.tsv
    cp pgx_plots/pgx_gene_country_sex_burden.tsv pgx_gene_country_sex_burden.tsv
    """
}
