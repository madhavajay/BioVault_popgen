// BioVault popgen: fast (numpy) gnomAD HGDP+1kGP PCA projection.
//
// Same I/O contract as biovault-popgen-gnomad-projection-1; differs only in
// the implementation under the hood:
//   - DDNA -> PLINK bed/bim/fam done in vectorized numpy (no tped, no Python
//     row loop, no plink2 --make-bed).
//   - hl.experimental.pc_project replaced by a numpy reader + matmul that
//     reproduces Hail's formula at float64 precision.
// Output (study_pca_projection.tsv, qc_report.txt, pca_projection.png) is
// bit-identical to the slow flow's output at typical cohort sizes.

nextflow.enable.dsl=2

workflow USER {
    take:
        context
        participants

    main:
        def participantTuples = participants.flatMap { record ->
            def validation = record.validation ?: [status: 'ok', message: '']
            if (validation.status?.toString() != 'ok') {
                println "[bv] WARNING: skipping participant ${record.participant_id}: genotype file ${validation.status} - ${validation.message}"
                return []
            }
            return [tuple(record.participant_id.toString(), file(record.genotype_file))]
        }
        def collected = participantTuples
            .collect(flat: false)
            .map { items ->
                if (items.isEmpty()) {
                    throw new IllegalArgumentException("No valid genotype files remained after BioVault input validation")
                }
                tuple(
                    items.collect { it[0] },
                    items.collect { it[1] }
                )
            }
        def projection = gnomad_projection_fast(collected)

    emit:
        projection_tsv = projection.projection_tsv
        qc_report = projection.qc_report
        projection_plot = projection.projection_plot
        errors = projection.errors
        warnings = projection.warnings
}

process gnomad_projection_fast {
    container 'ghcr.io/madhavajay/biovault-popgen:0.2.6-fast'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), path(genotype_files)

    output:
        path "study_pca_projection.tsv", emit: projection_tsv
        path "qc_report.txt",            emit: qc_report
        path "pca_projection.png",       emit: projection_plot, optional: true
        path "errors.tsv",               emit: errors, optional: true
        path "warnings.tsv",             emit: warnings, optional: true

    script:
    def staging = []
    participant_ids.eachWithIndex { pid, idx ->
        def fname = genotype_files[idx].getName()
        staging << "mkdir -p input/${pid} && ln -s \"../../${fname}\" \"input/${pid}/${fname}\""
    }
    """
    set -euo pipefail

    mkdir -p input
    ${staging.join('\n    ')}

    bash /opt/biovault/scripts/gnomad_projection_fast/run_fast_projection.sh \\
        "\${PWD}/input" \\
        "\${PWD}/work" \\
        "\${PWD}"
    """
}
