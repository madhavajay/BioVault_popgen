// BioVault popgen: joint BioVault + 1000 Genomes high-coverage PCA.

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
        def projection = hgp1k_projection_fast(collected)

    emit:
        scores = projection.scores
        study_projection_tsv = projection.study_projection_tsv
        allele_freqs = projection.allele_freqs
        qc_report = projection.qc_report
        projection_plot = projection.projection_plot
        errors = projection.errors
}

process hgp1k_projection_fast {
    container 'ghcr.io/madhavajay/biovault-popgen:0.2.0-fast'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), path(genotype_files)

    output:
        path "pca_scores.tsv",           emit: scores
        path "study_pca_projection.tsv", emit: study_projection_tsv
        path "allele_freqs.tsv",         emit: allele_freqs
        path "qc_report.txt",            emit: qc_report
        path "pca_projection.png",       emit: projection_plot, optional: true
        path "errors.tsv",               emit: errors, optional: true

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

    bash /opt/biovault/scripts/hgp1k_projection_fast/run_hgp1k_projection.sh \\
        "\${PWD}/input" \\
        "\${PWD}/work" \\
        "\${PWD}"
    """
}
