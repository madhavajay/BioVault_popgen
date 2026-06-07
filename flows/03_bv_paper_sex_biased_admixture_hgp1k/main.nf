// BioVault popgen: HGP1K-anchored sex-biased ADMIXTURE.
//
// Joint ADMIXTURE of the study cohort + a baked 1000 Genomes reference
// (AFR/EUR/SAS) at K=3,4,5, run on the combined genome and separately on
// autosomes vs X to detect sex-biased admixture (a component over-represented
// on X). The `sex` participant facet is materialized into sex_mapping.tsv and
// applied (plink2 --update-sex) so male X is encoded haploid — sex is NEVER
// inferred from the data.
//
// Two stages:
//   cohort_bed       (biosynth)          bvs cohort-bed -> study_raw PLINK BED
//   hgp1k_admixture  (biovault-admixture) merge w/ baked reference, QC+LD-prune
//                                         autosomes & X, ADMIXTURE, X-vs-auto.
// The baked reference BED (.docker/reference/hgp1k_admixture) and the analysis
// script live in the biovault-admixture image.

nextflow.enable.dsl=2

if (!params.containsKey('admixture_image')) { params.admixture_image = null }
if (!params.containsKey('biosynth_image'))  { params.biosynth_image  = null }

def BIOSYNTH_IMAGE  = System.getenv('BIOSYNTH_IMAGE') ?: (params.biosynth_image ?: 'ghcr.io/openmined/biosynth:0.1.32')
def ADMIXTURE_IMAGE = params.admixture_image ?: 'ghcr.io/madhavajay/biovault-admixture:latest'

def normalizeSex(String raw) { return (raw ?: '').trim() }

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
            def sex = normalizeSex((record.sex ?: facets.sex)?.toString())
            if (!sex) {
                println "[bv] WARNING: skipping participant ${record.participant_id}: missing required sex facet"
                return []
            }
            return [tuple(record.participant_id.toString(), sex, file(record.genotype_file))]
        }

        def collected = records
            .collect(flat: false)
            .map { items ->
                if (items.isEmpty()) {
                    throw new IllegalArgumentException("No valid participants with genotype files and a sex facet remained")
                }
                tuple(
                    items.collect { it[0] },
                    items.collect { it[1] },
                    items.collect { it[2] }
                )
            }

        def bed = cohort_bed(collected)
        def result = hgp1k_admixture(bed.bed_dir, bed.sex_mapping)

    emit:
        sex_bias_x_vs_auto       = result.sex_bias
        sex_bias_plot            = result.sex_bias_plot
        sex_bias_per_sample      = result.per_sample
        figures                  = result.figures
        figures_pdf              = result.figures_pdf
        component_labels         = result.component_labels
        qc_report                = result.qc_report
        labeled_q_files          = result.labeled_q
        admixture_q_files         = result.q_files
        admixture_p_files         = result.p_files
        admixture_logs           = result.adx_logs
        plink_bed                = bed.bed_dir
        plink_bed_archive        = bed.bed_archive
}

process cohort_bed {
    container BIOSYNTH_IMAGE
    containerOptions '--entrypoint=""'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), val(sexes), path(genotype_files)

    output:
        path "plink_bed",          emit: bed_dir
        path "plink_bed.tar.gz",   emit: bed_archive
        path "sex_mapping.tsv",    emit: sex_mapping

    script:
    def staging = []
    def sexRows = []
    participant_ids.eachWithIndex { pid, idx ->
        def fname = genotype_files[idx].getName()
        staging << "mkdir -p input/${pid} && ln -s \"../../${fname}\" \"input/${pid}/${fname}\""
        sexRows << "${pid}\t${sexes[idx]}"
    }
    def sexText = sexRows.join('\\n')
    """
    set -euo pipefail
    mkdir -p input plink_bed
    ${staging.join('\n    ')}
    { printf 'participant_id\\tsex\\n'; printf '%b\\n' "${sexText}"; } > sex_mapping.tsv
    bvs cohort-bed -i input \\
        --out-prefix plink_bed/genotypes \\
        --snp-info plink_bed/snp_info.tsv
    tar -czf plink_bed.tar.gz plink_bed
    """
}

process hgp1k_admixture {
    container ADMIXTURE_IMAGE
    containerOptions '--platform=linux/amd64 --entrypoint=""'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path bed_dir
        path sex_mapping

    output:
        path "sex_bias_x_vs_auto.tsv",          emit: sex_bias
        path "sex_bias_x_vs_auto.png",          emit: sex_bias_plot, optional: true
        path "sex_bias_x_vs_auto.pdf",          emit: sex_bias_plot_pdf, optional: true
        path "sex_bias_x_vs_auto_reference.tsv", emit: sex_bias_reference, optional: true
        path "sex_bias_per_sample_K*.tsv",      emit: per_sample, optional: true
        path "figure_sex_biased_admixture_K*.png", emit: figures, optional: true
        path "figure_sex_biased_admixture_K*.pdf", emit: figures_pdf, optional: true
        path "component_labels.tsv",            emit: component_labels
        path "qc_report.txt",                   emit: qc_report
        path "admixture_*_labeled_Q.tsv",       emit: labeled_q, optional: true
        path "admixture_*.Q",                   emit: q_files, optional: true
        path "admixture_*.P",                   emit: p_files, optional: true
        path "admixture_*.log",                 emit: adx_logs, optional: true

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

    mkdir -p sex_biased_admixture_hgp1k/scripts
    cp /opt/biovault/scripts/sex_biased_admixture_hgp1k/*.py sex_biased_admixture_hgp1k/scripts/

    python3 sex_biased_admixture_hgp1k/scripts/hgp1k_admixture.py \\
        --study-bed "${bed_dir}/genotypes" \\
        --sex-mapping "${sex_mapping}" \\
        --reference-dir "\${HGP1K_ADMIXTURE_REF:-/opt/biovault/reference/hgp1k_admixture}" \\
        --out-dir "\${PWD}" \\
        --k-values "\${BV_ADMIXTURE_K:-3,4,5}" \\
        --threads "\${BV_THREADS:-8}"
    """
}
