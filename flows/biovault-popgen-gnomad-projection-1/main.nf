// BioVault popgen: project per-participant genotypes onto the gnomAD
// HGDP+1kGP PCA space using v3.1 loadings baked into the image.

nextflow.enable.dsl=2

workflow USER {
    take:
        context
        participants

    main:
        def participantTuples = participants.map { record ->
            tuple(record.participant_id.toString(), file(record.genotype_file))
        }
        // Single emission of two aligned lists: ids and staged file paths.
        def collected = participantTuples
            .collect(flat: false)
            .map { items ->
                tuple(
                    items.collect { it[0] },
                    items.collect { it[1] }
                )
            }
        def projection = gnomad_projection(collected)

    emit:
        projection_tsv = projection.projection_tsv
        qc_report = projection.qc_report
        projection_plot = projection.projection_plot
}

process gnomad_projection {
    container 'biovault-popgen:0.1.0'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'copy'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), path(genotype_files)

    output:
        path "study_pca_projection.tsv", emit: projection_tsv
        path "qc_report.txt",            emit: qc_report
        path "pca_projection.png",       emit: projection_plot, optional: true

    script:
    def staging = []
    participant_ids.eachWithIndex { pid, idx ->
        def fname = genotype_files[idx].getName()
        staging << "mkdir -p input/${pid} && mv '${fname}' 'input/${pid}/'"
    }
    """
    set -euo pipefail

    mkdir -p input
    ${staging.join('\n    ')}

    bash /opt/biovault/scripts/gnomad_projection/run_flow_projection.sh \\
        "\${PWD}/input" \\
        "\${PWD}/work" \\
        "\${PWD}"
    """
}
