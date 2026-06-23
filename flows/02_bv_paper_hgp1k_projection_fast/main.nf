// BioVault popgen: joint BioVault + 1000 Genomes high-coverage PCA.
//
// Structure:
//   hgp1k_projection_fast (ghcr.io/madhavajay/biovault-popgen:0.2.7-fast)
//       Runs `bvs target-study-dosage` once to build an HGP1K-aligned study
//       dosage matrix plus study AF TSV, then builds the 1KGP + BioVault PCA.

nextflow.enable.dsl=2

if (!params.containsKey('hgp1k_image')) {
    params.hgp1k_image = null
}
def emitWorkers() {
    def raw = params.emit_workers ?: 0
    def n = raw as int
    return n > 0 ? n : 0
}

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
        closest_population_labels = projection.closest_population_labels
        closest_population_summary = projection.closest_population_summary
        pca_variants_used = projection.pca_variants_used
        pca_variants_dropped = projection.pca_variants_dropped
}

process hgp1k_projection_fast {
    container "${params.hgp1k_image ?: 'ghcr.io/madhavajay/biovault-popgen:0.2.7-fast'}"
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), path(genotype_files)

    output:
        path "pca_scores.tsv",                  emit: scores
        path "study_pca_projection.tsv",        emit: study_projection_tsv
        path "allele_freqs.tsv",                emit: allele_freqs
        path "qc_report.txt",                   emit: qc_report
        path "pca_projection.png",              emit: projection_plot, optional: true
        path "errors.tsv",                      emit: errors, optional: true
        path "closest_population_labels.tsv",   emit: closest_population_labels
        path "closest_population_summary.tsv",  emit: closest_population_summary, optional: true
        path "pca_variants_used.tsv",           emit: pca_variants_used
        path "pca_variants_dropped.tsv",        emit: pca_variants_dropped

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

    THREADS="${emitWorkers()}"   # 0 = all cores (bvs auto-detects)
    HGP1K_VARIANTS_TSV="\${HGP1K_VARIANTS_TSV:-/opt/biovault/reference/hgp1k/variants.tsv}"
    echo "[hgp1k] bvs target-study-dosage: \$(find -L input -type f | wc -l | tr -d ' ') participants, threads=\${THREADS}"

    bvs target-study-dosage -i input --threads "\${THREADS}" \\
        --variants-tsv "\${HGP1K_VARIANTS_TSV}" \\
        --dosage-npy bvs_study_dosage.npy \\
        --samples-tsv bvs_study_samples.tsv \\
        --study-variants-tsv bvs_study_variants.tsv \\
        --allele-freq-tsv bvs_study_allele_freq.tsv \\
        >/dev/null 2>bvs_study_dosage.warnings.log

    test -s bvs_study_dosage.npy
    test -s bvs_study_samples.tsv
    test -s bvs_study_allele_freq.tsv

    BVS_STUDY_AF_TSV="\${PWD}/bvs_study_allele_freq.tsv" \\
    BVS_STUDY_DOSAGE_NPY="\${PWD}/bvs_study_dosage.npy" \\
    BVS_STUDY_SAMPLES_TSV="\${PWD}/bvs_study_samples.tsv" \\
    bash /opt/biovault/scripts/hgp1k_projection_fast/run_hgp1k_projection.sh \\
        "\${PWD}/input" \\
        "\${PWD}/work" \\
        "\${PWD}"
    """
}
