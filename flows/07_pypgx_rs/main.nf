// BioVault popgen: PyPGx-rs PGx flow.
//
//   prepare_vcf  (BIOSYNTH_IMAGE, default ghcr.io/openmined/biosynth:0.1.32)
//       Passes VCF inputs through, or converts genotype TXT files to VCF.
//
//   pypgx_rs_pipeline  (ghcr.io/madhavajay/pypgx-rs:latest)
//       Runs pypgx-rs run-ngs-pipeline against each VCF.
//
//   aggregate_pypgx
//       Adds country facets and aggregates unique genotypes by country/gene.

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
def PYPGX_RS_IMAGE = params.pypgx_rs_image ?: 'ghcr.io/madhavajay/pypgx-rs:latest'
def POPGEN_IMAGE = System.getenv('POPGEN_IMAGE') ?: 'ghcr.io/madhavajay/biovault-popgen:0.2.4-fast'

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
        participant_possible_genotypes = aggregate.participant_possible_genotypes
        participant_possible_genotypes_normalized = aggregate.participant_possible_genotypes_normalized
        country_gene_genotype_counts = aggregate.country_gene_genotype_counts
        country_gene_genotype_counts_normalized = aggregate.country_gene_genotype_counts_normalized
        country_summary = aggregate.country_summary
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
        path "pypgx_participant_possible_genotypes.tsv", emit: participant_possible_genotypes
        path "pypgx_participant_possible_genotypes_normalized.tsv", emit: participant_possible_genotypes_normalized
        path "pypgx_country_gene_genotype_counts.tsv", emit: country_gene_genotype_counts
        path "pypgx_country_gene_genotype_counts_normalized.tsv", emit: country_gene_genotype_counts_normalized
        path "pypgx_country_summary.tsv", emit: country_summary
        path "pypgx_failures.tsv", emit: failures
        path "errors.tsv", emit: errors
        path "warnings.tsv", emit: warnings
        path "pypgx_pipeline.log", emit: pipeline_log

    script:
    """
    set -euo pipefail

    python3 - <<'PY'
import csv
import re
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

participant_rows = []
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

participant_rows.sort(key=lambda r: (r["participant_id"], r["gene"]))

with open("pypgx_participant_results.tsv", "w", newline="") as handle:
    fields = ["participant_id", "country", "gene", "status", "genotype", "phenotype", "error"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\\t")
    writer.writeheader()
    writer.writerows(participant_rows)

with open("pypgx_participant_possible_genotypes.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["participant_id", "country", "gene", "possible_genotypes"])
    for row in participant_rows:
        if row["genotype"]:
            writer.writerow([row["participant_id"], row["country"], row["gene"], row["genotype"]])

with open("pypgx_participant_possible_genotypes_normalized.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["participant_id", "country", "gene", "possible_genotypes"])
    for row in participant_rows:
        if row["genotype"]:
            genotype = re.sub(r"\\s*\\([^)]*\\)", "", row["genotype"]).strip()
            writer.writerow([row["participant_id"], row["country"], row["gene"], genotype])

def write_counts(path, normalized=False):
    counts = Counter()
    samples = set()
    for row in participant_rows:
        genotype = row["genotype"]
        if not genotype:
            continue
        if normalized:
            genotype = re.sub(r"\\s*\\([^)]*\\)", "", genotype).strip()
        key = (row["country"], row["gene"], genotype)
        counts[key] += 1
        samples.add((row["country"], row["gene"], row["participant_id"]))
    denom = Counter((country, gene) for country, gene, _pid in samples)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\\t")
        writer.writerow(["country", "gene", "possible_genotypes", "count", "sample_count", "frequency"])
        for (country, gene, genotype), count in sorted(counts.items()):
            d = denom[(country, gene)]
            writer.writerow([country, gene, genotype, count, d, f"{(count / d) if d else 0:.6f}"])

write_counts("pypgx_country_gene_genotype_counts.tsv", normalized=False)
write_counts("pypgx_country_gene_genotype_counts_normalized.tsv", normalized=True)

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
    handle.write(f"errors={len(errors)}\\n")

if failures:
    print(f"[bv] PyPGx: {len(failures)} participant(s) marked down (see pypgx_failures.tsv); "
          f"{len(set(r['participant_id'] for r in participant_rows))} succeeded", flush=True)
PY
    """
}
