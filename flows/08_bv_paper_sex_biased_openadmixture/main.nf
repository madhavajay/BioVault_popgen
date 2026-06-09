// BioVault popgen: HGP1K-anchored sex-biased OpenADMIXTURE.
//
// Joint OpenADMIXTURE of the study cohort + a baked 1000 Genomes reference
// (AFR/AMR/EAS/EUR/SAS) at K=5. Because OpenADMIXTURE has no ADMIXTURE-style
// --haploid option, the all-sample headline run is autosome-only, and the
// autosome-vs-X sex-bias comparison uses matched female-only study/reference
// samples. The `sex` participant facet is materialized into sex_mapping.tsv;
// sex is NEVER inferred from the data.
//
// Two stages:
//   cohort_bed            (biosynth)             bvs cohort-bed -> study_raw PLINK BED
//   hgp1k_openadmixture   (openadmixture.jl)     merge w/ baked reference, QC+LD-prune
//                                                all-auto and female-only auto/X.
// The BioVault OpenADMIXTURE image supplies the Julia runner, PLINK tools, and
// a baked HGP1K reference directory.

nextflow.enable.dsl=2

if (!params.containsKey('openadmixture_image')) { params.openadmixture_image = null }
if (!params.containsKey('biosynth_image'))  { params.biosynth_image  = null }
if (!params.containsKey('hgp1k_admixture_ref')) { params.hgp1k_admixture_ref = null }

def BIOSYNTH_IMAGE        = System.getenv('BIOSYNTH_IMAGE') ?: (params.biosynth_image ?: 'ghcr.io/openmined/biosynth:0.1.32')
def OPENADMIXTURE_IMAGE   = params.openadmixture_image ?: 'ghcr.io/madhavajay/biovault-openadmixture:latest'

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

        def analysis_script = file("${projectDir}/scripts/hgp1k_openadmixture.py")
        def reference_dir = params.hgp1k_admixture_ref ?: '/opt/biovault/reference/hgp1k_admixture'

        def bed = cohort_bed(collected)
        def result = hgp1k_openadmixture(bed.bed_dir, bed.sex_mapping, analysis_script, reference_dir)

    emit:
        sex_bias_x_vs_auto       = result.sex_bias
        sex_bias_plot            = result.sex_bias_plot
        sex_bias_plot_pdf        = result.sex_bias_plot_pdf
        sex_bias_per_sample      = result.per_sample
        figures                  = result.figures
        figures_pdf              = result.figures_pdf
        component_labels         = result.component_labels
        qc_report                = result.qc_report
        labeled_q_files          = result.labeled_q
        openadmixture_q_files     = result.q_files
        openadmixture_p_files     = result.p_files
        openadmixture_logs        = result.adx_logs
}

process cohort_bed {
    container BIOSYNTH_IMAGE
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), val(sexes), path(genotype_files)

    output:
        path "plink_bed",          emit: bed_dir
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
    echo "[bv] cohort_bed: building PLINK BED for ${participant_ids.size()} participants with sex facets"
    mkdir -p input plink_bed
    ${staging.join('\n    ')}
    { printf 'participant_id\\tsex\\n'; printf '%b\\n' "${sexText}"; } > sex_mapping.tsv
    bvs cohort-bed -i input \\
        --out-prefix plink_bed/genotypes \\
        --snp-info plink_bed/snp_info.tsv
    """
}

process hgp1k_openadmixture {
    container OPENADMIXTURE_IMAGE
    containerOptions '--entrypoint=""'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path bed_dir
        path sex_mapping
        path analysis_script
        val reference_dir

    output:
        path "sex_bias_x_vs_auto.tsv",          emit: sex_bias
        path "sex_bias_x_vs_auto.png",          emit: sex_bias_plot, optional: true
        path "sex_bias_x_vs_auto.pdf",          emit: sex_bias_plot_pdf, optional: true
        path "sex_bias_x_vs_auto_reference.tsv", emit: sex_bias_reference, optional: true
        path "sex_bias_per_sample_K*.tsv",      emit: per_sample, optional: true
        path "figure_sex_biased_openadmixture_K*.png", emit: figures, optional: true
        path "figure_sex_biased_openadmixture_K*.pdf", emit: figures_pdf, optional: true
        path "component_labels.tsv",            emit: component_labels
        path "qc_report.txt",                   emit: qc_report
        path "openadmixture_*_labeled_Q.tsv",   emit: labeled_q, optional: true
        path "openadmixture_*.Q",               emit: q_files, optional: true
        path "openadmixture_*.P",               emit: p_files, optional: true
        path "openadmixture_*.log",             emit: adx_logs, optional: true

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
    export JULIA_DEPOT_PATH="\${JULIA_DEPOT_PATH:-/tmp/.julia:/opt/julia-depot:/root/.julia}"

    mkdir -p sex_biased_openadmixture_hgp1k/scripts
    cp "${analysis_script}" sex_biased_openadmixture_hgp1k/scripts/hgp1k_openadmixture.py
    echo "[bv] hgp1k_openadmixture: study BED=${bed_dir}/genotypes"
    echo "[bv] hgp1k_openadmixture: reference=${reference_dir}"
    echo "[bv] hgp1k_openadmixture: K=\${BV_OPENADMIXTURE_K:-5} threads=\${BV_THREADS:-8}"

    python3 sex_biased_openadmixture_hgp1k/scripts/hgp1k_openadmixture.py \\
        --study-bed "${bed_dir}/genotypes" \\
        --sex-mapping "${sex_mapping}" \\
        --reference-dir "${reference_dir}" \\
        --out-dir "\${PWD}" \\
        --k-values "\${BV_OPENADMIXTURE_K:-5}" \\
        --threads "\${BV_THREADS:-8}"
    """
}
