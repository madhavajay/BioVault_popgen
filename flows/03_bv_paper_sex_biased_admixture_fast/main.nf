// BioVault popgen: sex-biased admixture (X-hemizygosity), `sex`-facet driven.
//
// Byte-identical fast sibling of the original sex_biased_admixture module
// (only the genotype loader is parallelized). The `sex` participant facet is
// materialized into a sex_mapping.tsv (participant_id -> sex) and handed to
// the analysis via BIOVAULT_SEX_MAPPING — exactly how the CLI
// (03_individual_level.sh --sex) feeds it. Sex is NEVER inferred from BAF.
//
// The image bakes the script *contents* flat into
//   /opt/biovault/scripts/sex_biased_admixture_fast/  (fast loader)
//   /opt/biovault/scripts/sex_biased_admixture/        (original analysis)
//   /opt/biovault/tools/                               (shared normalizer)
// fast_sex_biased_admixture.py resolves _BASE = parents[1] and
// ORIG_SCRIPTS = parents[2]/sex_biased_admixture/scripts, so both trees are
// reconstructed under the writable workdir before running.

nextflow.enable.dsl=2

def normalizeSex(String raw) {
    return (raw ?: '').trim()
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
                    throw new IllegalArgumentException("No valid participants with readable genotype files and sex facet remained")
                }
                tuple(
                    items.collect { it[0] },
                    items.collect { it[1] },
                    items.collect { it[2] }
                )
            }

        def result = sex_biased_admixture_fast(collected)

    emit:
        sex_bias_results  = result.sex_bias_results
        nmf_variant_filter_autosomes = result.nmf_variant_filter_autosomes
        nmf_variant_filter_x = result.nmf_variant_filter_x
        sex_bias_plot     = result.sex_bias_plot
        sex_bias_plot_pdf = result.sex_bias_plot_pdf
        pipeline_log      = result.pipeline_log
        errors            = result.errors
        warnings          = result.warnings
}

process sex_biased_admixture_fast {
    container 'ghcr.io/madhavajay/biovault-popgen:0.1.6-fast'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), val(sexes), path(genotype_files)

    output:
        path "sex_bias_results.tsv",                emit: sex_bias_results
        path "nmf_variant_filter_autosomes.tsv",    emit: nmf_variant_filter_autosomes, optional: true
        path "nmf_variant_filter_x.tsv",            emit: nmf_variant_filter_x, optional: true
        path "figure4_sex_biased_admixture.png",    emit: sex_bias_plot, optional: true
        path "figure4_sex_biased_admixture.pdf",    emit: sex_bias_plot_pdf, optional: true
        path "sex_biased_admixture.log",            emit: pipeline_log, optional: true
        path "errors.tsv",                          emit: errors, optional: true
        path "warnings.tsv",                        emit: warnings, optional: true

    script:
    def staging = []
    def sexMap = []
    participant_ids.eachWithIndex { pid, idx ->
        def fname = genotype_files[idx].getName()
        staging << "mkdir -p input/${pid} && ln -s \"../../${fname}\" \"input/${pid}/${fname}\""
        sexMap << "${pid}\t${sexes[idx]}"
    }
    def sexMapText = sexMap.join('\\n')
    """
    set -euo pipefail

    # passwd entry for the host UID (subprocesses use pwd-name lookups)
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

    # Materialize the `sex` facet into the mapping file the analysis reads.
    { printf 'participant_id\\tsex\\n'; printf '%b\\n' "${sexMapText}"; } \\
        > sex_mapping.tsv

    # Reconstruct both baked script trees so fast_sex_biased_admixture.py's
    # _BASE (parents[1]) and ORIG_SCRIPTS (parents[2]/sex_biased_admixture
    # /scripts) resolve inside the writable workdir.
    mkdir -p sex_biased_admixture_fast/scripts \\
             sex_biased_admixture/scripts
    cp /opt/biovault/scripts/sex_biased_admixture_fast/*.py sex_biased_admixture_fast/scripts/
    cp /opt/biovault/scripts/sex_biased_admixture/*.py      sex_biased_admixture/scripts/

    source /opt/conda/etc/profile.d/conda.sh
    conda activate biovault_popgen

    export BIOVAULT_DATA_DIR="\${PWD}/input"
    export BIOVAULT_SEX_MAPPING="\${PWD}/sex_mapping.tsv"
    python3 sex_biased_admixture_fast/scripts/fast_sex_biased_admixture.py

    # Hoist artefacts to the process root so publishDir picks them up.
    cp sex_biased_admixture_fast/results/sex_bias_results.tsv sex_bias_results.tsv
    [ -f sex_biased_admixture_fast/results/nmf_variant_filter_autosomes.tsv ] && \\
        cp sex_biased_admixture_fast/results/nmf_variant_filter_autosomes.tsv nmf_variant_filter_autosomes.tsv || true
    [ -f sex_biased_admixture_fast/results/nmf_variant_filter_x.tsv ] && \\
        cp sex_biased_admixture_fast/results/nmf_variant_filter_x.tsv nmf_variant_filter_x.tsv || true
    [ -f sex_biased_admixture_fast/plots/figure4_sex_biased_admixture.png ] && \\
        cp sex_biased_admixture_fast/plots/figure4_sex_biased_admixture.png figure4_sex_biased_admixture.png || true
    [ -f sex_biased_admixture_fast/plots/figure4_sex_biased_admixture.pdf ] && \\
        cp sex_biased_admixture_fast/plots/figure4_sex_biased_admixture.pdf figure4_sex_biased_admixture.pdf || true
    if [ -f sex_biased_admixture_fast/logs/sex_biased_admixture.log ]; then
        cp sex_biased_admixture_fast/logs/sex_biased_admixture.log sex_biased_admixture.log
    elif [ -f sex_biased_admixture/logs/sex_biased_admixture.log ]; then
        cp sex_biased_admixture/logs/sex_biased_admixture.log sex_biased_admixture.log
    fi
    [ -f sex_biased_admixture_fast/logs/errors.tsv ] && cp sex_biased_admixture_fast/logs/errors.tsv errors.tsv || true
    [ -f sex_biased_admixture_fast/logs/warnings.tsv ] && cp sex_biased_admixture_fast/logs/warnings.tsv warnings.tsv || true
    """
}
