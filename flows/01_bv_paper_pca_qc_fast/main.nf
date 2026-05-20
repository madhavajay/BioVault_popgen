// BioVault popgen: within-cohort QC + PCA sanity check.
//
// Calls pca_qc_fast/scripts/fast_pipeline.py. That script's BASE_DIR is
// computed from its own file location, so we copy the script tree into the
// Nextflow workdir before running, then publish the artefacts out of the
// resulting `data/`, `plots/`, and `logs/` subdirs.

nextflow.enable.dsl=2

workflow USER {
    take:
        context
        participants

    main:
        def participantTuples = participants.map { record ->
            tuple(
                record.participant_id.toString(),
                file(record.genotype_file),
                (record.country?.toString()?.trim() ?: '(unset)'),
                (record.sex?.toString()?.trim() ?: '(unset)')
            )
        }
        def collected = participantTuples
            .collect(flat: false)
            .map { items ->
                // Aggregate-only facet proof (counts, never participant IDs).
                // Best-effort: no required_facets, so pca_qc still runs on a
                // bare samplesheet — facets show '(unset)' until the loader
                // carries them, then real per-value counts appear here.
                def countryCounts = items.collect { it[2] }.countBy { it }.sort()
                def sexCounts     = items.collect { it[3] }.countBy { it }.sort()
                println "[bv] participants: ${items.size()}"
                println "[bv] facet country: " +
                    countryCounts.collect { k, v -> "${k}=${v}" }.join(', ')
                println "[bv] facet sex: " +
                    sexCounts.collect { k, v -> "${k}=${v}" }.join(', ')
                tuple(
                    items.collect { it[0] },
                    items.collect { it[1] }
                )
            }
        def qc = pca_qc_fast(collected)

    emit:
        eigenvec = qc.eigenvec
        eigenval = qc.eigenval
        snp_info = qc.snp_info
        pca_pc12_plot = qc.pca_pc12_plot
        pca_pc34_plot = qc.pca_pc34_plot
        pipeline_log = qc.pipeline_log
}

process pca_qc_fast {
    container 'biovault-popgen:0.1.1'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'copy'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), path(genotype_files)

    output:
        path "pca.eigenvec",     emit: eigenvec
        path "pca.eigenval",     emit: eigenval
        path "snp_info.tsv",     emit: snp_info
        path "pca_pc1_pc2.png",  emit: pca_pc12_plot, optional: true
        path "pca_pc3_pc4.png",  emit: pca_pc34_plot, optional: true
        path "fast_pipeline.log", emit: pipeline_log, optional: true

    script:
    def staging = []
    participant_ids.eachWithIndex { pid, idx ->
        def fname = genotype_files[idx].getName()
        staging << "mkdir -p input/${pid} && mv '${fname}' 'input/${pid}/'"
    }
    """
    set -euo pipefail

    # passwd entry for the host UID (matches the slow runner's setup; needed
    # by any subprocess that uses pwd-name lookups)
    USER_ID="\$(id -u)"; GROUP_ID="\$(id -g)"
    if ! getent passwd "\${USER_ID}" >/dev/null 2>&1; then
        echo "biovault:x:\${USER_ID}:\${GROUP_ID}:biovault:/tmp:/bin/bash" >> /etc/passwd
    fi
    if ! getent group "\${GROUP_ID}" >/dev/null 2>&1; then
        echo "biovault:x:\${GROUP_ID}:" >> /etc/group
    fi
    export HOME=/tmp

    mkdir -p input
    ${staging.join('\n    ')}

    # The image bakes the script *contents* flat into
    # /opt/biovault/scripts/pca_qc_fast/ (fast_pipeline.py, genoio.py, ...).
    # fast_pipeline.py uses BASE_DIR = Path(__file__).parents[1], so it must
    # live at <base>/scripts/fast_pipeline.py for outputs to land in
    # <base>/{data,plots,logs}. Reconstruct that layout in the writable workdir.
    mkdir -p pca_qc_fast/scripts pca_qc_fast/data pca_qc_fast/plots pca_qc_fast/logs
    cp /opt/biovault/scripts/pca_qc_fast/*.py pca_qc_fast/scripts/

    source /opt/conda/etc/profile.d/conda.sh
    conda activate biovault_popgen

    export BIOVAULT_DATA_DIR="\${PWD}/input"
    python3 pca_qc_fast/scripts/fast_pipeline.py

    # Hoist final artefacts to the process root so publishDir picks them up.
    cp pca_qc_fast/data/pca/pca.eigenvec       pca.eigenvec
    cp pca_qc_fast/data/pca/pca.eigenval       pca.eigenval
    cp pca_qc_fast/data/merged/snp_info.tsv    snp_info.tsv
    [ -f pca_qc_fast/plots/pca_pc1_pc2.png ] && cp pca_qc_fast/plots/pca_pc1_pc2.png pca_pc1_pc2.png || true
    [ -f pca_qc_fast/plots/pca_pc3_pc4.png ] && cp pca_qc_fast/plots/pca_pc3_pc4.png pca_pc3_pc4.png || true
    [ -f pca_qc_fast/logs/fast_pipeline.log ] && cp pca_qc_fast/logs/fast_pipeline.log fast_pipeline.log || true
    """
}
