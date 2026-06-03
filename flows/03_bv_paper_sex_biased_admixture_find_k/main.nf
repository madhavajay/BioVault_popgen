// BioVault popgen: reported-ancestry ADMIXTURE K search.

nextflow.enable.dsl=2

if (!params.containsKey('find_k_image')) {
    params.find_k_image = null
}
if (!params.containsKey('biosynth_image')) {
    params.biosynth_image = null
}

def BIOSYNTH_IMAGE = System.getenv('BIOSYNTH_IMAGE') ?: (params.biosynth_image ?: 'ghcr.io/openmined/biosynth:0.1.31')
def FIND_K_IMAGE = params.find_k_image ?: 'ghcr.io/madhavajay/biovault-popgen:0.2.0-fast'

def ancestryValue(record) {
    def facets = record.facets ?: [:]
    return (
        record.self_reported_ancestry ?:
        record.reported_ancestry ?:
        record.ancestry ?:
        record.ethnicity ?:
        facets.self_reported_ancestry ?:
        facets.reported_ancestry ?:
        facets.ancestry ?:
        facets.ethnicity ?:
        ''
    ).toString().trim()
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
            def ancestry = ancestryValue(record)
            if (!ancestry) {
                println "[bv] WARNING: skipping participant ${record.participant_id}: missing reported ancestry facet"
                return []
            }
            return [tuple(record.participant_id.toString(), ancestry, file(record.genotype_file))]
        }

        def collected = records
            .collect(flat: false)
            .map { items ->
                if (items.isEmpty()) {
                    throw new IllegalArgumentException("No valid participants with genotype files and reported ancestry remained")
                }
                tuple(
                    items.collect { it[0] },
                    items.collect { it[1] },
                    items.collect { it[2] }
                )
            }

        def bed = cohort_bed_reported_ancestry(collected)
        def result = admixture_find_k(bed.bed_dir, bed.ancestry_map)

    emit:
        cv_errors = result.cv_errors
        k_summary = result.k_summary
        k_summary_plot = result.k_summary_plot
        selected_k = result.selected_k
        component_anchor_means = result.component_anchor_means
        component_labels = result.component_labels
        ancestry_anchor_samples = result.ancestry_anchor_samples
        reported_ancestry_normalized = result.reported_ancestry_normalized
        find_k_report = result.find_k_report
}

process cohort_bed_reported_ancestry {
    container BIOSYNTH_IMAGE
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), val(ancestries), path(genotype_files)

    output:
        path "plink_bed", emit: bed_dir
        path "reported_ancestry.tsv", emit: ancestry_map

    script:
    def staging = []
    def ancestryRows = []
    participant_ids.eachWithIndex { pid, idx ->
        def fname = genotype_files[idx].getName()
        staging << "mkdir -p input/${pid} && ln -s \"../../${fname}\" \"input/${pid}/${fname}\""
        ancestryRows << "${pid}\t${ancestries[idx].replace('\t', ' ').replace('\n', ' ')}"
    }
    def ancestryText = ancestryRows.join('\\n')
    """
    set -euo pipefail
    mkdir -p input plink_bed
    ${staging.join('\n    ')}
    { printf 'participant_id\\tself_reported_ancestry\\n'; printf '%b\\n' "${ancestryText}"; } \\
        > reported_ancestry.tsv
    bvs cohort-bed -i input \\
        --out-prefix plink_bed/genotypes \\
        --snp-info plink_bed/snp_info.tsv
    """
}

process admixture_find_k {
    container FIND_K_IMAGE
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path bed_dir
        path ancestry_map

    output:
        path "admixture_cv_errors.tsv", emit: cv_errors
        path "admixture_k_summary.tsv", emit: k_summary
        path "admixture_k_summary.png", emit: k_summary_plot, optional: true
        path "selected_k.tsv", emit: selected_k
        path "component_anchor_means.tsv", emit: component_anchor_means
        path "component_labels.tsv", emit: component_labels
        path "ancestry_anchor_samples.tsv", emit: ancestry_anchor_samples
        path "reported_ancestry_normalized.tsv", emit: reported_ancestry_normalized
        path "find_k_report.txt", emit: find_k_report
        path "admixture_K*.Q", optional: true
        path "admixture_K*.P", optional: true
        path "admixture_K*_labeled_Q.tsv", optional: true

    script:
    """
    set -euo pipefail

    USER_ID="\$(id -u)"; GROUP_ID="\$(id -g)"
    if ! getent passwd "\${USER_ID}" >/dev/null 2>&1; then
        echo "biovault:x:\${USER_ID}:\${GROUP_ID}:biovault:/tmp:/bin/bash" >> /etc/passwd
    fi
    if ! getent group "\${GROUP_ID}" >/dev/null 2>&1; then
        echo "biovault:x:\${GROUP_ID}:" >> /etc/group
    fi
    export HOME=/tmp

    mkdir -p sex_biased_admixture_find_k/scripts
    cp /opt/biovault/scripts/sex_biased_admixture_find_k/*.py sex_biased_admixture_find_k/scripts/

    source /opt/conda/etc/profile.d/conda.sh
    conda activate biovault_popgen

    python3 sex_biased_admixture_find_k/scripts/find_k.py \\
        --bed-prefix "${bed_dir}/genotypes" \\
        --ancestry-map "${ancestry_map}" \\
        --out-dir "\${PWD}" \\
        --k-min "\${BV_ADMIXTURE_K_MIN:-2}" \\
        --k-max "\${BV_ADMIXTURE_K_MAX:-8}" \\
        --reps "\${BV_ADMIXTURE_REPS:-1}" \\
        --threads "\${BV_THREADS:-8}"
    """
}
