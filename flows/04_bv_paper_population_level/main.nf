// BioVault popgen: population-level FST + AIMs, driven by a `country` facet.
//
// Two containers, orchestrated by Nextflow (the desktop runner pre-pulls every
// per-process `container '...'`):
//
//   split_allele_freq  (BIOSYNTH_IMAGE, default ghcr.io/openmined/biosynth:0.1.27)
//       per-participant `bvs emit-long`, then per-country `bvs aggregate-long`
//       -> allele_freq_<country_norm>.tsv
//
//   population_fst_aims  (ghcr.io/madhavajay/biovault-popgen:0.1.7-fast)
//       FST (load/merge -> WC84 -> visualise) then AIMs (merge w/ bundled
//       gnomAD ref -> differential SNPs -> AIMs panels)
//
// The split step is just `bvs` calls, inlined here (biosynth image's only
// tool) — no script shipped with the flow. The FST/AIMs step consumes the
// scripts BAKED into biovault-popgen at /opt/biovault/scripts/population_level
// (canonical source: 04_population_level/fst_aims_fast/scripts; Dockerfile
// bakes it like the other analysis trees). Editing those scripts needs an
// image rebuild, same as 01_/02_.
//
// Country normalization (must match scripts/popset.py): trim -> lowercase ->
// non-alphanumeric runs to "_" -> strip leading/trailing "_".

nextflow.enable.dsl=2

def BIOSYNTH_IMAGE = System.getenv('BIOSYNTH_IMAGE') ?: 'ghcr.io/openmined/biosynth:0.1.27'

def normalizeCountry(String raw) {
    return (raw ?: '')
        .trim()
        .toLowerCase()
        .replaceAll(/[^a-z0-9]+/, '_')
        .replaceAll(/^_+|_+$/, '')
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
            def country = normalizeCountry((record.country ?: facets.country)?.toString())
            if (!country) {
                println "[bv] WARNING: skipping participant ${record.participant_id}: missing required country facet"
                return []
            }
            return [tuple(record.participant_id.toString(), country, file(record.genotype_file))]
        }

        def collected = records
            .collect(flat: false)
            .map { items ->
                if (items.isEmpty()) {
                    throw new IllegalArgumentException("No valid participants with readable genotype files and country facet remained")
                }
                tuple(
                    items.collect { it[0] },
                    items.collect { it[1] },
                    items.collect { it[2] }
                )
            }

        def split = split_allele_freq(collected)
        def result = population_fst_aims(split.af_dir, split.populations)

    emit:
        allele_freqs       = result.allele_freqs
        fst_matrix         = result.fst_matrix
        merged_allele_freq = result.merged_allele_freq
        master_af_table    = result.master_af_table
        differential_snps  = result.differential_snps
        aims_combined      = result.aims_combined
        summary            = result.summary
        errors             = result.errors
        warnings           = result.warnings
}

process split_allele_freq {
    container BIOSYNTH_IMAGE
    // The biosynth image sets ENTRYPOINT ["bvs"], so Nextflow's
    // `<image> /bin/bash .command.run` becomes `bvs /bin/bash …` and bvs
    // rejects it ("unrecognized subcommand '/bin/bash'"). Clear the
    // entrypoint so Nextflow's own bash wrapper runs (mirrors the CLI's
    // `docker run --entrypoint "" … bvs …` in 04_run_allele_freq.sh).
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), val(countries), path(genotype_files)

    output:
        path "af_out",          emit: af_dir
        env  POPULATIONS,       emit: populations

    script:
    def mapping = []
    def staging = []
    participant_ids.eachWithIndex { pid, idx ->
        def orig = genotype_files[idx].getName()
        def staged = "${pid}__${orig}"
        mapping << "${pid}\t${countries[idx]}\t${staged}"
        staging << "ln -s \"../${orig}\" \"geno/${staged}\""
    }
    def mappingText = mapping.join('\\n')
    """
    set -euo pipefail

    mkdir -p geno af_out
    ${staging.join('\n    ')}

    printf '%b\\n' "${mappingText}" > mapping.tsv

    # Inlined split: just bvs (the biosynth image's only tool) — no baked
    # script needed in this container. Mirrors 04_population_level/fst_aims_fast
    # /scripts/split_allele_freq.sh (kept as the by-hand reference).
    mkdir -p bvlr
    : > successful_mapping.tsv
    : > skipped_participants.tsv
    total_participants=\$(wc -l < mapping.tsv | tr -d ' ')
    processed_participants=0
    participant_start=\$(date +%s)
    while IFS="\$(printf '\\t')" read -r pid country fname; do
        [ -n "\${pid}" ] || continue
        processed_participants=\$((processed_participants + 1))
        src="geno/\${fname}"
        if [ ! -s "\${src}" ]; then
            echo "WARNING: skipping participant \${pid}: missing or empty genotype \${src}" >&2
            printf '%s\\t%s\\t%s\\tmissing_or_empty_genotype\\n' "\${pid}" "\${country}" "\${fname}" >> skipped_participants.tsv
        elif ! bvs emit-long --input "\${src}" --output "bvlr/\${pid}.bvlr" \\
            --participant "\${pid}" >/dev/null; then
            echo "WARNING: skipping participant \${pid}: bvs emit-long failed for \${src}" >&2
            printf '%s\\t%s\\t%s\\temit_long_failed\\n' "\${pid}" "\${country}" "\${fname}" >> skipped_participants.tsv
            rm -f "bvlr/\${pid}.bvlr"
        elif [ ! -s "bvlr/\${pid}.bvlr" ]; then
            echo "WARNING: skipping participant \${pid}: bvs emit-long produced an empty .bvlr" >&2
            printf '%s\\t%s\\t%s\\tempty_bvlr\\n' "\${pid}" "\${country}" "\${fname}" >> skipped_participants.tsv
            rm -f "bvlr/\${pid}.bvlr"
        else
            printf '%s\\t%s\\t%s\\n' "\${pid}" "\${country}" "\${fname}" >> successful_mapping.tsv
        fi
        if [ "\$((processed_participants % 50))" -eq 0 ] || [ "\${processed_participants}" -eq "\${total_participants}" ]; then
            now=\$(date +%s)
            elapsed=\$((now - participant_start))
            [ "\${elapsed}" -gt 0 ] || elapsed=1
            rate=\$(awk -v done="\${processed_participants}" -v sec="\${elapsed}" 'BEGIN { printf "%.1f", done / sec }')
            remaining=\$((total_participants - processed_participants))
            eta=\$(awk -v rem="\${remaining}" -v done="\${processed_participants}" -v sec="\${elapsed}" 'BEGIN { if (done > 0) printf "%.0f", rem * sec / done; else print "0" }')
            echo "Processed \${processed_participants}/\${total_participants} participants (rate \${rate}/s, ETA \${eta}s)"
        fi
    done < mapping.tsv

    [ -s successful_mapping.tsv ] || { echo "ERROR: no participants produced usable .bvlr files" >&2; exit 1; }

    FAIL=0
    : > successful_countries.txt
    total_countries=\$(cut -f2 successful_mapping.tsv | sort -u | wc -l | tr -d ' ')
    processed_countries=0
    for country in \$(cut -f2 successful_mapping.tsv | sort -u); do
        processed_countries=\$((processed_countries + 1))
        cdir="agg_\${country}"; mkdir -p "\${cdir}"; n=0
        while IFS="\$(printf '\\t')" read -r pid c _; do
            [ "\${c}" = "\${country}" ] || continue
            [ -s "bvlr/\${pid}.bvlr" ] && { cp "bvlr/\${pid}.bvlr" "\${cdir}/"; n=\$((n+1)); }
        done < successful_mapping.tsv
        [ "\${n}" -gt 0 ] || { echo "WARNING: country '\${country}' produced 0 .bvlr; skipping" >&2; continue; }
        out="af_out/allele_freq_\${country}.tsv"
        echo "Aggregating country \${processed_countries}/\${total_countries}: \${country} (\${n} participants) -> \${out}"
        if ! bvs aggregate-long --input "\${cdir}" \\
            --matrix-tsv "matrix_\${country}.tsv" \\
            --allele-freq-tsv "\${out}" >/dev/null; then
            echo "WARNING: bvs aggregate-long failed for country '\${country}'; skipping" >&2
            rm -f "\${out}"
            continue
        fi
        if [ ! -s "\${out}" ]; then
            echo "WARNING: empty AF for country '\${country}'; skipping" >&2
            rm -f "\${out}"
            continue
        fi
        printf '%s\\n' "\${country}" >> successful_countries.txt
    done
    [ -s successful_countries.txt ] || { echo "ERROR: no countries produced aggregate allele-frequency files" >&2; exit 1; }

    { printf 'participant_id\\tcountry\\n'; cut -f1,2 successful_mapping.tsv; } \\
        > af_out/country_map.tsv
    if [ -s skipped_participants.tsv ]; then
        { printf 'participant_id\\tfile\\tseverity\\tcode\\tmessage\\n'; \\
          awk -F '\\t' 'BEGIN { OFS="\\t" } { print \$1, \$3, "ERROR", \$4, "country=" \$2 }' skipped_participants.tsv; } \\
            > af_out/errors.tsv
    else
        printf 'participant_id\\tfile\\tseverity\\tcode\\tmessage\\n' > af_out/errors.tsv
    fi
    printf 'file\\tline_no\\tseverity\\tcode\\tmessage\\traw_line\\n' > af_out/warnings.tsv

    POPULATIONS="\$(paste -sd, successful_countries.txt)"
    """
}

process population_fst_aims {
    container 'ghcr.io/madhavajay/biovault-popgen:0.1.7-fast'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path af_dir
        val  populations

    output:
        path "allele_freq_*.tsv",                  emit: allele_freqs
        path "fst_matrix.tsv",                     emit: fst_matrix
        path "merged_allele_freq_annotated.tsv",   emit: merged_allele_freq
        path "master_af_table.tsv",                emit: master_af_table
        path "all_outliers_long.tsv",              emit: differential_snps
        path "aims_combined.tsv",                  emit: aims_combined
        path "population_level_summary.txt",       emit: summary
        path "errors.tsv",                         emit: errors, optional: true
        path "warnings.tsv",                       emit: warnings, optional: true
        path "*.png",                              emit: plots, optional: true
        path "*.pdf",                              emit: plots_pdf, optional: true

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

    bash /opt/biovault/scripts/population_level/run_pipeline.sh \\
        "\${PWD}/${af_dir}" \\
        "\${PWD}/work" \\
        "\${PWD}" \\
        "${populations}"

    for f in "\${PWD}/${af_dir}"/allele_freq_*.tsv; do
        [ -f "\${f}" ] || { echo "ERROR: no aggregate allele_freq_*.tsv files found in ${af_dir}" >&2; exit 1; }
        cp "\${f}" .
    done
    [ -f "\${PWD}/${af_dir}/errors.tsv" ] && cp "\${PWD}/${af_dir}/errors.tsv" errors.tsv || true
    [ -f "\${PWD}/${af_dir}/warnings.tsv" ] && cp "\${PWD}/${af_dir}/warnings.tsv" warnings.tsv || true
    """
}
