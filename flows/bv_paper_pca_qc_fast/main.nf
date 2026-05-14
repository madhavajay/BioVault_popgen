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
            tuple(record.participant_id.toString(), file(record.genotype_file))
        }
        def collected = participantTuples
            .collect(flat: false)
            .map { items ->
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
        log = qc.log
}

process pca_qc_fast {
    container 'biovault-popgen:0.1.0'
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
        path "fast_pipeline.log", emit: log, optional: true

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

    # pca_qc_fast's BASE_DIR = parent of the script. Copy the script tree
    # into a writable location so `data/`, `plots/`, `logs/` land here.
    cp -r /opt/biovault/scripts/pca_qc_fast pca_qc_fast
    mkdir -p pca_qc_fast/data pca_qc_fast/plots pca_qc_fast/logs

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
