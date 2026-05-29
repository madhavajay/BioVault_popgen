// BioVault popgen: preflight genotype file QC.
//
// This flow stages every provided genotype file, runs the line-level QC
// checker, and publishes reports. Bad genotype contents are reported in
// issues.tsv; the process itself exits 0 by default so all files can be
// inspected in one run. Infrastructure failures still use the app-provided
// Nextflow error strategy.

nextflow.enable.dsl=2

def csvCell(value) {
    def text = (value == null ? '' : value.toString())
    return '"' + text.replace('"', '""') + '"'
}

def safeDir(value) {
    return value.toString().replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def qcBatchSize() {
    def raw = params.qc_batch_size ?: 50
    def size = raw as int
    return size > 0 ? size : 50
}

def qcProgressEvery() {
    def raw = params.qc_progress_every ?: 50
    def size = raw as int
    return size > 0 ? size : 50
}

def qcWorkers() {
    // Per-batch parallel workers inside qc_all_files.py. Default 1 so it does
    // not oversubscribe cores against Nextflow's batch-level parallelism; the
    // vectorized normalizer already speeds each file. Raise via params.qc_workers
    // when running few/large batches.
    def raw = params.qc_workers ?: 1
    def size = raw as int
    return size > 0 ? size : 1
}

workflow USER {
    take:
        context
        participants

    main:
        def participantRecords = participants.map { record ->
            def facets = [:]
            if (record.facets) {
                facets.putAll(record.facets)
            }
            if (record.country != null && !facets.containsKey('country')) {
                facets.country = record.country
            }
            if (record.sex != null && !facets.containsKey('sex')) {
                facets.sex = record.sex
            }
            def validation = record.validation ?: [status: 'ok', message: '']
            def pathText = record.genotype_path ?: record.genotype_file?.toString()
            tuple(
                record.participant_id.toString(),
                pathText?.toString() ?: '',
                validation.status?.toString() ?: 'ok',
                validation.message?.toString() ?: '',
                facets.collectEntries { key, value -> [(key.toString()): (value?.toString()?.trim() ?: '')] }
            )
        }

        def counted = participantRecords
            .map { it }
            .collect(flat: false)
            .map { items ->
                println "[bv] qc_all_files participants: ${items.size()}"
                println "[bv] qc_all_files valid staged inputs: ${items.count { it[2] == 'ok' }}"
                println "[bv] qc_all_files input issues: ${items.count { it[2] != 'ok' }}"
                return items
            }

        def validRecords = counted
            .flatMap { items -> items.findAll { it[2] == 'ok' } }
            .map { pid, pathText, status, message, facets ->
                tuple(pid, file(pathText), facets)
            }
            .collect(flat: false)
            .flatMap { items ->
                def batches = items.collate(qcBatchSize())
                println "[bv] qc_all_files batches: ${batches.size()} (batch_size=${qcBatchSize()}, progress_every=${qcProgressEvery()})"
                batches.withIndex().collect { batch, idx ->
                    tuple(
                        idx + 1,
                        batches.size(),
                        batch.collect { it[0] },
                        batch.collect { it[1] },
                        batch.collect { it[2] }
                    )
                }
            }

        def invalidRecords = counted
            .flatMap { items -> items.findAll { it[2] != 'ok' } }
            .map { pid, pathText, status, message, facets ->
                tuple(pid, pathText, status, message, facets)
            }

        def validQc = qc_file_batch(validRecords)
        def invalidQc = qc_input_issue(invalidRecords)
        def qc = merge_qc_reports(validQc.report_dir.mix(invalidQc.report_dir).collect(flat: false))

    emit:
        file_summary = qc.file_summary
        issues = qc.issues
        errors = qc.errors
        file_errors = qc.file_errors
        row_errors = qc.row_errors
        warnings = qc.warnings
        quality_issues = qc.quality_issues
        summary_json = qc.summary_json
        report = qc.report
        facet_values = qc.facet_values
        facet_summary = qc.facet_summary
        pipeline_log = qc.pipeline_log
}

process qc_file_batch {
    container 'ghcr.io/madhavajay/biovault-popgen:0.2.0-fast'
    stageInMode 'symlink'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(batch_number), val(batch_total), val(participant_ids), path(genotype_files), val(facet_maps)

    output:
        path "${prefix}", emit: report_dir

    script:
    prefix = "batch_${batch_number}"
    def batchLabel = "${batch_number}/${batch_total}"
    def progressEvery = qcProgressEvery()
    def facetNames = facet_maps.collectMany { it.keySet().collect { key -> key.toString() } }.unique().sort()
    def header = ['participant_id', 'genotype_file'] + facetNames
    def rows = [header.collect { csvCell(it) }.join(',')]
    def staging = []
    participant_ids.eachWithIndex { participant_id, idx ->
        def participantDir = safeDir(participant_id)
        def genotype_file = genotype_files[idx]
        def facets = facet_maps[idx]
        def staged = "input/${participantDir}/${genotype_file.getName()}"
        def values = [participant_id, staged] + facetNames.collect { name -> facets[name] ?: '' }
        rows << values.collect { csvCell(it) }.join(',')
        staging << "mkdir -p ${shellQuote(prefix)}/input/${shellQuote(participantDir)} && ln -s ../../../${shellQuote(genotype_file.getName())} ${shellQuote(prefix)}/${shellQuote(staged)}"
    }
    def samplesheetText = rows.join("\n")
    def samplesheetShell = samplesheetText.replace("'", "'\"'\"'")
    """
    set -euo pipefail

    # passwd entry for the host UID (subprocesses can use pwd-name lookups)
    USER_ID="\$(id -u)"; GROUP_ID="\$(id -g)"
    if ! getent passwd "\${USER_ID}" >/dev/null 2>&1; then
        echo "biovault:x:\${USER_ID}:\${GROUP_ID}:biovault:/tmp:/bin/bash" >> /etc/passwd
    fi
    if ! getent group "\${GROUP_ID}" >/dev/null 2>&1; then
        echo "biovault:x:\${GROUP_ID}:" >> /etc/group
    fi
    export HOME=/tmp

    mkdir -p ${shellQuote(prefix)}/qc_output
    ${staging.join('\n    ')}

    printf '%s\\n' '${samplesheetShell}' > ${shellQuote(prefix)}/selected_participants.csv

    set +e
    python /opt/biovault/scripts/qc_all_files/qc_all_files.py \\
        --samplesheet ${shellQuote(prefix)}/selected_participants.csv \\
        --output-dir ${shellQuote(prefix)}/qc_output \\
        --batch-label ${shellQuote(batchLabel)} \\
        --progress-every ${progressEvery} \\
        --workers ${qcWorkers()}
    qc_status="\$?"
    set -e

    if [ "\${qc_status}" -ne 0 ] && [ ! -s ${shellQuote(prefix)}/qc_output/file_summary.tsv ]; then
        echo "QC script exited \${qc_status} for batch ${prefix}; creating diagnostic fallback outputs" >&2
        {
            printf 'participant_id\\tfile\\tsource\\tdetected_format\\tstatus\\treadable\\traw_rows\\tnormalized_rows\\tunique_variants\\tfacet_count\\tfacet_missing_count\\tmissing_facets\\tfacets_json\\terrors\\twarnings\\tmessage\\n'
            printf '%s\\t%s\\tnextflow\\tunknown\\tFAIL\\tTrue\\t0\\t0\\t0\\t0\\t0\\t\\t{}\\t1\\t0\\tQC script exited before writing normal reports\\n' '${prefix}' '${prefix}'
        } > ${shellQuote(prefix)}/qc_output/file_summary.tsv
        {
            printf 'participant_id\\tfile\\tdetected_format\\tline_number\\tseverity\\tcode\\tmessage\\tline\\n'
            printf '%s\\t%s\\tunknown\\t0\\tERROR\\tQC_PROCESS_EXIT\\tqc_all_files.py exited with status %s before writing normal reports\\t\\n' '${prefix}' '${prefix}' "\${qc_status}"
        } > ${shellQuote(prefix)}/qc_output/issues.tsv
        cp ${shellQuote(prefix)}/qc_output/issues.tsv ${shellQuote(prefix)}/qc_output/errors.tsv
        cp ${shellQuote(prefix)}/qc_output/issues.tsv ${shellQuote(prefix)}/qc_output/file_errors.tsv
        printf 'participant_id\\tfile\\tdetected_format\\tline_number\\tseverity\\tcode\\tmessage\\tline\\n' > ${shellQuote(prefix)}/qc_output/row_errors.tsv
        printf 'participant_id\\tfile\\tdetected_format\\tline_number\\tseverity\\tcode\\tmessage\\tline\\n' > ${shellQuote(prefix)}/qc_output/warnings.tsv
        printf 'participant_id\\tfile\\tdetected_format\\tline_number\\tseverity\\tcode\\tmessage\\tline\\n' > ${shellQuote(prefix)}/qc_output/quality_issues.tsv
        printf '{\\n  "files": 1,\\n  "pass": 0,\\n  "warn": 0,\\n  "fail": 1,\\n  "errors": 1,\\n  "warnings": 0,\\n  "quality_issues": 0,\\n  "facets": [],\\n  "facet_missing_values": 0\\n}\\n' > ${shellQuote(prefix)}/qc_output/summary.json
        printf 'participant_id\\tfile\\tfacet\\tvalue\\tis_missing\\n' > ${shellQuote(prefix)}/qc_output/facet_values.tsv
        printf 'facet\\tvalue\\tcount\\n' > ${shellQuote(prefix)}/qc_output/facet_summary.tsv
        touch ${shellQuote(prefix)}/qc_output/qc_all_files.log
    fi

    cp ${shellQuote(prefix)}/qc_output/file_summary.tsv ${shellQuote(prefix)}/file_summary.tsv
    cp ${shellQuote(prefix)}/qc_output/issues.tsv ${shellQuote(prefix)}/issues.tsv
    cp ${shellQuote(prefix)}/qc_output/errors.tsv ${shellQuote(prefix)}/errors.tsv
    cp ${shellQuote(prefix)}/qc_output/file_errors.tsv ${shellQuote(prefix)}/file_errors.tsv
    cp ${shellQuote(prefix)}/qc_output/row_errors.tsv ${shellQuote(prefix)}/row_errors.tsv
    cp ${shellQuote(prefix)}/qc_output/warnings.tsv ${shellQuote(prefix)}/warnings.tsv
    cp ${shellQuote(prefix)}/qc_output/quality_issues.tsv ${shellQuote(prefix)}/quality_issues.tsv
        cp ${shellQuote(prefix)}/qc_output/summary.json ${shellQuote(prefix)}/summary.json
        [ -f ${shellQuote(prefix)}/qc_output/report.txt ] && cp ${shellQuote(prefix)}/qc_output/report.txt ${shellQuote(prefix)}/report.txt || true
        cp ${shellQuote(prefix)}/qc_output/facet_values.tsv ${shellQuote(prefix)}/facet_values.tsv
    cp ${shellQuote(prefix)}/qc_output/facet_summary.tsv ${shellQuote(prefix)}/facet_summary.tsv
    cp ${shellQuote(prefix)}/qc_output/qc_all_files.log ${shellQuote(prefix)}/qc_all_files.log
    """
}

process qc_input_issue {
    container 'ghcr.io/madhavajay/biovault-popgen:0.2.0-fast'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_id), val(genotype_path), val(validation_status), val(validation_message), val(facets)

    output:
        path "${prefix}", emit: report_dir

    script:
    prefix = "${safeDir(participant_id)}_input_issue"
    def facetNames = facets.keySet().collect { it.toString() }.sort()
    def header = ['participant_id', 'genotype_file'] + facetNames
    def rows = [header.collect { csvCell(it) }.join(',')]
    def values = [participant_id, genotype_path] + facetNames.collect { name -> facets[name] ?: '' }
    rows << values.collect { csvCell(it) }.join(',')
    def samplesheetText = rows.join("\n")
    def samplesheetShell = samplesheetText.replace("'", "'\"'\"'")
    """
    set -euo pipefail
    mkdir -p ${shellQuote(prefix)}/qc_output
    printf '%s\\n' '${samplesheetShell}' > ${shellQuote(prefix)}/selected_participants.csv
    echo "Input validation issue for participant ${participant_id}: ${validation_status} ${validation_message}" >&2

    python /opt/biovault/scripts/qc_all_files/qc_all_files.py \\
        --samplesheet ${shellQuote(prefix)}/selected_participants.csv \\
        --output-dir ${shellQuote(prefix)}/qc_output \\
        --batch-label ${shellQuote("${prefix}/input_issue")} \\
        --progress-every ${qcProgressEvery()}

    cp ${shellQuote(prefix)}/qc_output/file_summary.tsv ${shellQuote(prefix)}/file_summary.tsv
    cp ${shellQuote(prefix)}/qc_output/issues.tsv ${shellQuote(prefix)}/issues.tsv
    cp ${shellQuote(prefix)}/qc_output/errors.tsv ${shellQuote(prefix)}/errors.tsv
    cp ${shellQuote(prefix)}/qc_output/file_errors.tsv ${shellQuote(prefix)}/file_errors.tsv
    cp ${shellQuote(prefix)}/qc_output/row_errors.tsv ${shellQuote(prefix)}/row_errors.tsv
    cp ${shellQuote(prefix)}/qc_output/warnings.tsv ${shellQuote(prefix)}/warnings.tsv
    cp ${shellQuote(prefix)}/qc_output/quality_issues.tsv ${shellQuote(prefix)}/quality_issues.tsv
    cp ${shellQuote(prefix)}/qc_output/summary.json ${shellQuote(prefix)}/summary.json
    [ -f ${shellQuote(prefix)}/qc_output/report.txt ] && cp ${shellQuote(prefix)}/qc_output/report.txt ${shellQuote(prefix)}/report.txt || true
    cp ${shellQuote(prefix)}/qc_output/facet_values.tsv ${shellQuote(prefix)}/facet_values.tsv
    cp ${shellQuote(prefix)}/qc_output/facet_summary.tsv ${shellQuote(prefix)}/facet_summary.tsv
    cp ${shellQuote(prefix)}/qc_output/qc_all_files.log ${shellQuote(prefix)}/qc_all_files.log
    """
}

process merge_qc_reports {
    container 'ghcr.io/madhavajay/biovault-popgen:0.2.0-fast'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        path report_dirs

    output:
        path "file_summary.tsv",  emit: file_summary
        path "issues.tsv",        emit: issues
        path "errors.tsv",        emit: errors
        path "file_errors.tsv",   emit: file_errors
        path "row_errors.tsv",    emit: row_errors
        path "warnings.tsv",      emit: warnings
        path "quality_issues.tsv", emit: quality_issues
        path "summary.json",      emit: summary_json
        path "report.txt",        emit: report
        path "facet_values.tsv",  emit: facet_values
        path "facet_summary.tsv", emit: facet_summary
        path "qc_all_files.log",  emit: pipeline_log

    script:
    def dirs = report_dirs instanceof List ? report_dirs : [report_dirs]
    def fileSummaries = dirs.collect { shellQuote("${it}/file_summary.tsv") }.join(' ')
    def issues = dirs.collect { shellQuote("${it}/issues.tsv") }.join(' ')
    def facetValues = dirs.collect { shellQuote("${it}/facet_values.tsv") }.join(' ')
    def logs = dirs.collect { shellQuote("${it}/qc_all_files.log") }.join(' ')
    """
    set -euo pipefail
    python /opt/biovault/scripts/qc_all_files/merge_qc_reports.py \\
        --output-dir merged \\
        --file-summary ${fileSummaries} \\
        --issues ${issues} \\
        --facet-values ${facetValues} \\
        --logs ${logs}

    cp merged/file_summary.tsv file_summary.tsv
    cp merged/issues.tsv issues.tsv
    cp merged/errors.tsv errors.tsv
    cp merged/file_errors.tsv file_errors.tsv
    cp merged/row_errors.tsv row_errors.tsv
    cp merged/warnings.tsv warnings.tsv
    cp merged/quality_issues.tsv quality_issues.tsv
    cp merged/summary.json summary.json
    cp merged/report.txt report.txt
    cp merged/facet_values.tsv facet_values.tsv
    cp merged/facet_summary.tsv facet_summary.tsv
    cp merged/qc_all_files.log qc_all_files.log
    """
}
