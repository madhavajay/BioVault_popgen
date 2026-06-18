// ExVitae APOL1 risk report flow.

nextflow.enable.dsl=2

if (!params.containsKey('exvitae_image')) {
    params.exvitae_image = null
}
if (!params.containsKey('analysis_max_duration_ms')) {
    params.analysis_max_duration_ms = '30000'
}

def EXVITAE_CONTAINER = System.getenv('EXVITAE_IMAGE') ?: (params.exvitae_image ?: 'ghcr.io/madhavajay/exvitae:0.2.4')
def ASSAY_ID = 'apol1'
def OUTPUT_PREFIX = 'apol1'
def NF_PARAMS = (params.containsKey('nextflow') && params.nextflow) ? params.nextflow : [error_strategy: 'terminate', max_retries: 0]
def RESULTS_DIR = params.containsKey('results_dir') ? params.results_dir : 'results'

def countryValue(record) {
    def facets = record.facets ?: [:]
    def country = (
        record.country ?:
        record.Country ?:
        record.country_code ?:
        record.countryCode ?:
        facets.country ?:
        facets.Country ?:
        facets.country_code ?:
        facets.countryCode ?:
        'Unknown'
    ).toString().trim()
    return country ?: 'Unknown'
}

def safeId(value) {
    return value.toString().replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def envInt = { name, fallback ->
    def raw = System.getenv(name)
    if (raw && raw.isInteger()) {
        return raw.toInteger()
    }
    return fallback
}

def DEFAULT_REPORT_MAX_FORKS = envInt('EXVITAE_APOL1_MAX_FORKS', envInt('EXVITAE_REPORT_MAX_FORKS', 10))
def REPORT_MAX_FORKS = NF_PARAMS.report_max_forks ?: (NF_PARAMS.max_forks ?: DEFAULT_REPORT_MAX_FORKS)

workflow USER {
    take:
        context
        participants

    main:
        def participantItems = participants.flatMap { record ->
            def validation = record.validation ?: [status: 'ok', message: '']
            if (validation.status?.toString() != 'ok') {
                println "[bv] WARNING: skipping participant ${record.participant_id}: genotype file ${validation.status} - ${validation.message}"
                return []
            }
            def country = countryValue(record)
            if (country == 'Unknown') {
                println "[bv] WARNING: participant ${record.participant_id}: missing country facet; using Unknown"
            }
            return [tuple(record.participant_id.toString(), country, file(record.genotype_file))]
        }

        def checked = participantItems
            .ifEmpty {
                throw new IllegalArgumentException("No valid participants with readable genotype files remained")
            }
        def perParticipantReports = exvitae_report(checked)
        def report_outputs = aggregate_reports(perParticipantReports.report_dir.collect())

    emit:
        participant_reports = report_outputs.participant_reports
        country_aggregates = report_outputs.country_aggregates
        observations = report_outputs.observations
        reports = report_outputs.reports
        analysis = report_outputs.analysis
        participant_analysis = report_outputs.participant_analysis
        country_analysis_counts = report_outputs.country_analysis_counts
        country_variant_counts = report_outputs.country_variant_counts
}

process exvitae_report {
    container EXVITAE_CONTAINER
    stageInMode 'copy'
    tag { participant_id }
    errorStrategy { NF_PARAMS.error_strategy }
    maxRetries { NF_PARAMS.max_retries }
    maxForks REPORT_MAX_FORKS

    input:
        tuple val(participant_id), val(country), path(input_file)

    output:
        path "${prefix}", emit: report_dir

    script:
    prefix = safeId(participant_id)
    def inputName = input_file.getName()
    """
    set -euo pipefail
    mkdir -p ${shellQuote(prefix)}
    report_input="\${PWD}/${inputName}"
    lower="\$(printf '%s' ${shellQuote(inputName)} | tr '[:upper:]' '[:lower:]')"
    case "\${lower}" in
        *.gz|*.bgz)
            INPUT_PATH="\${PWD}/${inputName}" \\
            OUTPUT_DIR=${shellQuote(prefix)} \\
            ASSAY_ID=${shellQuote(ASSAY_ID)} \\
            ANALYSIS_MAX_DURATION_MS=${shellQuote(params.analysis_max_duration_ms)} \\
            python3 - <<'PY'
import gzip
import os
import shutil
import subprocess
import sys

input_path = os.environ["INPUT_PATH"]
output_dir = os.environ["OUTPUT_DIR"]
assay_id = os.environ["ASSAY_ID"]
analysis_max_duration_ms = os.environ["ANALYSIS_MAX_DURATION_MS"]

fd = os.memfd_create("exvitae-decompressed-input", 0)
try:
    with gzip.open(input_path, "rb") as src:
        with os.fdopen(os.dup(fd), "wb", closefd=True) as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
    os.lseek(fd, 0, os.SEEK_SET)
    result = subprocess.run(
        [
            "exvitae-report",
            assay_id,
            "--input-file",
            f"/proc/self/fd/{fd}",
            "--output-dir",
            output_dir,
            "--detect-sex",
            "--analysis-max-duration-ms",
            analysis_max_duration_ms,
        ],
        pass_fds=(fd,),
    )
    sys.exit(result.returncode)
finally:
    os.close(fd)
PY
            ;;
        *)
            exvitae-report ${shellQuote(ASSAY_ID)} \\
              --input-file "\${report_input}" \\
              --output-dir ${shellQuote(prefix)} \\
              --detect-sex \\
              --analysis-max-duration-ms ${params.analysis_max_duration_ms}
            ;;
    esac
    { printf 'participant_id\\tcountry\\n'; printf '%s\\t%s\\n' ${shellQuote(participant_id)} ${shellQuote(country)}; } \\
      > ${shellQuote("${prefix}/metadata.tsv")}
    """
}

process aggregate_reports {
    container EXVITAE_CONTAINER
    publishDir RESULTS_DIR, mode: 'copy', overwrite: true
    stageInMode 'copy'
    errorStrategy 'terminate'
    maxRetries { NF_PARAMS.max_retries }

    input:
        path report_dirs

    output:
        path "participants", emit: participant_reports
        path "observations.tsv", emit: observations
        path "reports.jsonl", emit: reports
        path "analysis.jsonl", emit: analysis
        path "${OUTPUT_PREFIX}_participant_analysis.tsv", emit: participant_analysis
        path "${OUTPUT_PREFIX}_country_analysis_counts.tsv", emit: country_analysis_counts
        path "${OUTPUT_PREFIX}_country_variant_counts.tsv", emit: country_variant_counts
        path "countries", emit: country_aggregates

    script:
    """
    set -euo pipefail

    python3 - <<'PY'
import csv
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path

OBS_HEADER = "participant_id\\tassay_id\\tassay_version\\tvariant_key\\trsid\\tassembly\\tchrom\\tpos_start\\tpos_end\\tref\\talt\\tkind\\tmatch_status\\tcoverage_status\\tcall_status\\tgenotype\\tgenotype_display\\tzygosity\\tref_count\\talt_count\\tdepth\\tgenotype_quality\\tallele_balance\\toutcome\\tevidence_type\\tevidence_raw\\tfacets\\n"

def safe_id(value):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return safe.strip("._") or "unknown"

OUTPUT_PREFIX = "${OUTPUT_PREFIX}"

def read_metadata(report_dir):
    path = report_dir / "metadata.tsv"
    metadata = {"participant_id": report_dir.name, "country": "Unknown"}
    if not path.is_file():
        return metadata
    for line in path.read_text(encoding="utf-8").splitlines()[1:2]:
        parts = line.split("\\t")
        if len(parts) > 0 and parts[0].strip():
            metadata["participant_id"] = parts[0].strip()
        if len(parts) > 1 and parts[1].strip():
            metadata["country"] = parts[1].strip()
        return metadata
    return metadata

def scalar(value):
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return str(value)

NOISY_ANALYSIS_KEYS = {
    "participant_id", "system", "confidence", "genotypes", "limitations",
    "interpretation", "notes",
}

def primary_analysis_values(row, emitted_keys):
    system = scalar(row.get("system")).strip()
    blood_type = scalar(row.get("blood_type")).strip()
    phenotype = scalar(row.get("phenotype")).strip()
    if system:
        if blood_type:
            return [(system, "blood_type", blood_type)]
        if phenotype:
            return [(system, "phenotype", phenotype)]

    keys = emitted_keys or [
        key for key, value in row.items()
        if key not in NOISY_ANALYSIS_KEYS and not isinstance(value, (dict, list))
    ]
    status_keys = [
        key for key in keys
        if key not in NOISY_ANALYSIS_KEYS and (key.endswith("_status") or key in {"status", "result", "phenotype", "blood_type"})
    ]
    candidate_keys = status_keys or [
        key for key in keys
        if key not in NOISY_ANALYSIS_KEYS
    ]
    values = []
    for key in candidate_keys:
        if key not in row:
            continue
        value = scalar(row.get(key)).strip()
        if value:
            values.append(("", scalar(key), value))
    return values[:1]

participant_analysis_rows = []
observation_rows = []

report_dirs = sorted(
    path for path in Path(".").iterdir()
    if path.is_dir() and path.name not in {"participants", "countries"}
)

Path("participants").mkdir(exist_ok=True)
Path("countries").mkdir(exist_ok=True)
reports = Path("reports.jsonl").open("w", encoding="utf-8")
analysis = Path("analysis.jsonl").open("w", encoding="utf-8")
obs = Path("observations.tsv").open("w", encoding="utf-8")
obs.write(OBS_HEADER)
country_obs_started = {}

for report_dir in report_dirs:
    name = report_dir.name
    metadata = read_metadata(report_dir)
    participant_id = metadata["participant_id"]
    country = metadata["country"]
    target = Path("participants") / name
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(report_dir, target)

    obs_path = report_dir / "observations.tsv"
    if obs_path.is_file():
        obs_lines = obs_path.read_text(encoding="utf-8").splitlines()
        for line in obs_lines[1:]:
            if line:
                obs.write(line + "\\n")
        if obs_lines:
            for row in csv.DictReader(obs_lines, delimiter="\\t"):
                observation_rows.append({
                    "participant_id": participant_id,
                    "country": country,
                    "assay_id": scalar(row.get("assay_id")),
                    "variant_key": scalar(row.get("variant_key")),
                    "rsid": scalar(row.get("rsid")),
                    "genotype_display": scalar(row.get("genotype_display")),
                    "outcome": scalar(row.get("outcome")),
                })

    for source, handle in ((report_dir / "reports.jsonl", reports), (report_dir / "analysis.jsonl", analysis)):
        if source.is_file():
            text = source.read_text(encoding="utf-8")
            if text:
                handle.write(text)
                if not text.endswith("\\n"):
                    handle.write("\\n")

    analysis_path = report_dir / "analysis.jsonl"
    if analysis_path.is_file():
        for line in analysis_path.read_text(encoding="utf-8").splitlines():
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            assay_id = scalar(payload.get("assay_id"))
            analysis_id = scalar(payload.get("analysis_id"))
            emitted_keys = [
                scalar(item.get("key"))
                for item in payload.get("emits", [])
                if isinstance(item, dict) and item.get("key") is not None
            ]
            rows = payload.get("rows") if isinstance(payload.get("rows"), list) else []
            for row in rows:
                if not isinstance(row, dict):
                    continue
                for result_scope, result_key, result_value in primary_analysis_values(row, emitted_keys):
                    participant_analysis_rows.append({
                        "participant_id": participant_id,
                        "country": country,
                        "assay_id": assay_id,
                        "analysis_id": analysis_id,
                        "result_scope": result_scope,
                        "result_key": result_key,
                        "result_value": result_value,
                    })

    country_root = Path("countries") / safe_id(country)
    country_participants = country_root / "participants"
    country_participants.mkdir(parents=True, exist_ok=True)
    country_target = country_participants / name
    if country_target.exists():
        shutil.rmtree(country_target)
    shutil.copytree(report_dir, country_target)

    country_obs = country_root / "observations.tsv"
    if not country_obs_started.get(country_root):
        country_obs.write_text(OBS_HEADER, encoding="utf-8")
        country_obs_started[country_root] = True
    if obs_path.is_file():
        with country_obs.open("a", encoding="utf-8") as handle:
            for line in obs_path.read_text(encoding="utf-8").splitlines()[1:]:
                if line:
                    handle.write(line + "\\n")
    for filename in ("reports.jsonl", "analysis.jsonl"):
        source = report_dir / filename
        if source.is_file():
            text = source.read_text(encoding="utf-8")
            if text:
                with (country_root / filename).open("a", encoding="utf-8") as handle:
                    handle.write(text)
                    if not text.endswith("\\n"):
                        handle.write("\\n")

reports.close()
analysis.close()
obs.close()

participant_fields = ["participant_id", "country", "assay_id", "analysis_id", "result_scope", "result_key", "result_value"]
with open(f"{OUTPUT_PREFIX}_participant_analysis.tsv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=participant_fields, delimiter="\\t")
    writer.writeheader()
    writer.writerows(sorted(
        participant_analysis_rows,
        key=lambda row: (row["country"], row["assay_id"], row["analysis_id"], row["result_scope"], row["result_key"], row["participant_id"], row["result_value"]),
    ))

denominators = defaultdict(set)
value_participants = defaultdict(set)
for row in participant_analysis_rows:
    base = (row["country"], row["assay_id"], row["analysis_id"], row["result_scope"], row["result_key"])
    denominators[base].add(row["participant_id"])
    value_participants[base + (row["result_value"],)].add(row["participant_id"])

with open(f"{OUTPUT_PREFIX}_country_analysis_counts.tsv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["country", "assay_id", "analysis_id", "result_scope", "result_key", "result_value", "count", "sample_count", "frequency"])
    for key, participants in sorted(value_participants.items()):
        country, assay_id, analysis_id, result_scope, result_key, result_value = key
        sample_count = len(denominators[(country, assay_id, analysis_id, result_scope, result_key)])
        count = len(participants)
        frequency = (count / sample_count) if sample_count else 0
        writer.writerow([country, assay_id, analysis_id, result_scope, result_key, result_value, count, sample_count, f"{frequency:.6f}"])

variant_denominators = defaultdict(set)
variant_participants = defaultdict(set)
for row in observation_rows:
    base = (row["country"], row["assay_id"], row["variant_key"], row["rsid"])
    variant_denominators[base].add(row["participant_id"])
    variant_participants[base + (row["genotype_display"], row["outcome"])].add(row["participant_id"])

with open(f"{OUTPUT_PREFIX}_country_variant_counts.tsv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\\t")
    writer.writerow(["country", "assay_id", "variant_key", "rsid", "genotype_display", "outcome", "count", "sample_count", "frequency"])
    for key, participants in sorted(variant_participants.items()):
        country, assay_id, variant_key, rsid, genotype_display, outcome = key
        sample_count = len(variant_denominators[(country, assay_id, variant_key, rsid)])
        count = len(participants)
        frequency = (count / sample_count) if sample_count else 0
        writer.writerow([country, assay_id, variant_key, rsid, genotype_display, outcome, count, sample_count, f"{frequency:.6f}"])
PY
    """
}
