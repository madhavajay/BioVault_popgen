// BioVault popgen: population-level FST + AIMs, driven by a `country` facet.
//
// Two containers, orchestrated by Nextflow (the desktop runner pre-pulls every
// per-process `container '...'`):
//
//   split_allele_freq  (ghcr.io/openmined/biosynth:0.1.23)
//       per-participant `bvs emit-long`, then per-country `bvs aggregate-long`
//       -> allele_freq_<country_norm>.tsv
//
//   population_fst_aims  (biovault-popgen:0.1.1)
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
        def records = participants.map { record ->
            def country = normalizeCountry(record.country?.toString())
            if (!country) {
                throw new IllegalArgumentException(
                    "Participant ${record.participant_id} has an empty 'country' " +
                    "facet after normalization. The flow declares " +
                    "required_facets: [country]; the samplesheet should have " +
                    "been rejected upstream."
                )
            }
            tuple(record.participant_id.toString(), country, file(record.genotype_file))
        }

        def collected = records
            .collect(flat: false)
            .map { items ->
                tuple(
                    items.collect { it[0] },
                    items.collect { it[1] },
                    items.collect { it[2] }
                )
            }

        def split = split_allele_freq(collected)
        def result = population_fst_aims(split.af_dir, split.populations)

    emit:
        country_map        = result.country_map
        fst_matrix         = result.fst_matrix
        merged_allele_freq = result.merged_allele_freq
        master_af_table    = result.master_af_table
        differential_snps  = result.differential_snps
        aims_combined      = result.aims_combined
        summary            = result.summary
}

process split_allele_freq {
    container 'ghcr.io/openmined/biosynth:0.1.23'
    // The biosynth image sets ENTRYPOINT ["bvs"], so Nextflow's
    // `<image> /bin/bash .command.run` becomes `bvs /bin/bash …` and bvs
    // rejects it ("unrecognized subcommand '/bin/bash'"). Clear the
    // entrypoint so Nextflow's own bash wrapper runs (mirrors the CLI's
    // `docker run --entrypoint "" … bvs …` in 04_run_allele_freq.sh).
    containerOptions '--entrypoint=""'
    stageInMode 'copy'
    errorStrategy { params.nextflow.error_strategy }
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
        staging << "mv '${orig}' 'geno/${staged}'"
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
    while IFS="\$(printf '\\t')" read -r pid country fname; do
        [ -n "\${pid}" ] || continue
        src="geno/\${fname}"
        [ -s "\${src}" ] || { echo "ERROR: missing genotype \${pid}: \${src}" >&2; exit 1; }
        bvs emit-long --input "\${src}" --output "bvlr/\${pid}.bvlr" \\
            --participant "\${pid}" >/dev/null
    done < mapping.tsv

    FAIL=0
    for country in \$(cut -f2 mapping.tsv | sort -u); do
        cdir="agg_\${country}"; mkdir -p "\${cdir}"; n=0
        while IFS="\$(printf '\\t')" read -r pid c _; do
            [ "\${c}" = "\${country}" ] || continue
            [ -s "bvlr/\${pid}.bvlr" ] && { cp "bvlr/\${pid}.bvlr" "\${cdir}/"; n=\$((n+1)); }
        done < mapping.tsv
        [ "\${n}" -gt 0 ] || { echo "ERROR: country '\${country}' produced 0 .bvlr" >&2; FAIL=1; continue; }
        out="af_out/allele_freq_\${country}.tsv"
        echo "  \${country} (\${n} participants) -> \${out}"
        bvs aggregate-long --input "\${cdir}" \\
            --matrix-tsv "matrix_\${country}.tsv" \\
            --allele-freq-tsv "\${out}" >/dev/null
        [ -s "\${out}" ] || { echo "ERROR: empty AF for \${country}" >&2; FAIL=1; }
    done
    [ "\${FAIL}" -eq 0 ] || { echo "ERROR: per-country split failed" >&2; exit 1; }

    { printf 'participant_id\\tcountry\\n'; cut -f1,2 mapping.tsv; } \\
        > af_out/country_map.tsv

    POPULATIONS="\$(cut -f2 mapping.tsv | sort -u | paste -sd, -)"
    """
}

process population_fst_aims {
    container 'biovault-popgen:0.1.1'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'copy'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        path af_dir
        val  populations

    output:
        path "country_map.tsv",                    emit: country_map
        path "fst_matrix.tsv",                     emit: fst_matrix
        path "merged_allele_freq_annotated.tsv",   emit: merged_allele_freq
        path "master_af_table.tsv",                emit: master_af_table
        path "all_outliers_long.tsv",              emit: differential_snps
        path "aims_combined.tsv",                  emit: aims_combined
        path "population_level_summary.txt",       emit: summary
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
    """
}
