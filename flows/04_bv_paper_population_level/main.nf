// BioVault popgen: population-level FST + AIMs, driven by a `country` facet.
//
// Structure (one Nextflow task PER COUNTRY so each is cached/resumable on its
// own, and disk stays bounded to a single country's .bvlr at a time):
//
//   country_allele_freq  (BIOSYNTH_IMAGE, default ghcr.io/openmined/biosynth:0.1.31)
//       For ONE country: a single fused `bvs fast-allele-freq` parses every
//       participant and aggregates -> allele_freq_<country>.tsv in one pass.
//       No per-participant .bvlr is written (peak disk ~0), parsing runs in
//       parallel across cores. Output is content-identical to the previous
//       emit-long + aggregate-long path (loci/counts/rsid), only sorted by
//       locus_key; downstream FST/AIMs join on locus_key so results are
//       unchanged. Runs one country at a time by default (params.country_forks).
//
//   population_fst_aims  (ghcr.io/madhavajay/biovault-popgen:0.2.2-fast)
//       Gathers every country's outputs, then FST (load/merge -> WC84 ->
//       visualise) and AIMs (merge w/ bundled gnomAD ref -> differential SNPs
//       -> AIMs panels). FST/AIMs scripts are BAKED into biovault-popgen at
//       /opt/biovault/scripts/population_level (editing them needs a rebuild).
//
// Country normalization (must match scripts/popset.py): trim -> lowercase ->
// non-alphanumeric runs to "_" -> strip leading/trailing "_".

nextflow.enable.dsl=2

def BIOSYNTH_IMAGE = System.getenv('BIOSYNTH_IMAGE') ?: 'ghcr.io/openmined/biosynth:0.1.31'

def normalizeCountry(String raw) {
    return (raw ?: '')
        .trim()
        .toLowerCase()
        .replaceAll(/[^a-z0-9]+/, '_')
        .replaceAll(/^_+|_+$/, '')
}

def emitWorkers() {
    // Concurrent `bvs emit-long` processes within one country task.
    // 0/unset = use all cores at runtime (nproc); set params.emit_workers to cap.
    def raw = params.emit_workers ?: 0
    def n = raw as int
    return n > 0 ? n : 0
}

def countryForks() {
    // How many country tasks Nextflow runs at once. Default 1 = one country at a
    // time (peak disk ~one country's .bvlr, ~9 GB for ~180 full-size files).
    // Raise via params.country_forks only if you have the disk + cores.
    def raw = params.country_forks ?: 1
    def n = raw as int
    return n > 0 ? n : 1
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

        // Fan out: one task per country -> independently cached/resumable, and
        // disk bounded to one country at a time.
        def byCountry = records
            .map { pid, country, gfile -> tuple(country, pid, gfile) }
            .groupTuple()
            .map { country, pids, gfiles ->
                println "[bv] country '${country}': ${pids.size()} participants"
                tuple(country, pids, gfiles)
            }

        def perCountry = country_allele_freq(byCountry)
        def result = population_fst_aims(perCountry.country_dir.collect())

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

process country_allele_freq {
    container BIOSYNTH_IMAGE
    // biosynth image sets ENTRYPOINT ["bvs"]; clear it so Nextflow's bash runs.
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    tag { country }
    // Nextflow 25.x strict syntax rejects a user-function call as a directive
    // argument, so inline countryForks()'s logic (null/0 -> 1, else the value).
    maxForks( (params.country_forks ?: 1) as int )
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(country), val(participant_ids), path(genotype_files)

    output:
        path "cc_${country}", emit: country_dir

    script:
    def staging = []
    def mapping = []
    participant_ids.eachWithIndex { pid, idx ->
        def orig = genotype_files[idx].getName()
        def staged = "${pid}__${orig}"
        mapping << "${pid}\t${staged}"
        staging << "ln -s \"../${orig}\" \"geno/${staged}\""
    }
    def mappingText = mapping.join('\\n')
    """
    set -euo pipefail
    cc="cc_${country}"
    mkdir -p geno "\${cc}"
    ${staging.join('\n    ')}
    printf '%b\\n' "${mappingText}" > mapping.tsv   # pid<TAB>staged-filename

    THREADS="${emitWorkers()}"   # 0 = all cores (bvs auto-detects)
    n_in=\$(wc -l < mapping.tsv | tr -d ' ')
    echo "[${country}] fast-allele-freq: \${n_in} participants, threads=\${THREADS}"

    # Pre-scan only to record diagnostics: missing/empty genotype files -> skipped,
    # the rest -> country participant map. bvs also tolerates bad files at runtime
    # (a single unreadable file is skipped, it never aborts the country).
    : > "\${cc}/skipped.tsv"
    : > "\${cc}/country_map.part.tsv"
    while IFS="\$(printf '\\t')" read -r pid fname; do
        [ -n "\${pid}" ] || continue
        if [ ! -s "geno/\${fname}" ]; then
            printf '%s\\t%s\\t%s\\tmissing_or_empty_genotype\\n' "\${pid}" "${country}" "\${fname}" >> "\${cc}/skipped.tsv"
        else
            printf '%s\\t%s\\n' "\${pid}" "${country}" >> "\${cc}/country_map.part.tsv"
        fi
    done < mapping.tsv

    n=\$(wc -l < "\${cc}/country_map.part.tsv" | tr -d ' ')
    : > "\${cc}/warnings.log"
    if [ "\${n}" -gt 0 ]; then
        out="\${cc}/allele_freq_${country}.tsv"
        echo "[${country}] aggregate: \${n} participants -> allele_freq_${country}.tsv"
        # ONE fused pass: parse every participant + aggregate -> allele_freq.
        # No per-participant .bvlr written (peak disk ~0); parallel across cores;
        # per-row warnings are coalesced into warnings.log, not streamed per row.
        # --max-ram-gb 0 = auto (80% of detected/cgroup RAM); set params.max_ram_gb to hard-cap.
        if bvs fast-allele-freq -i geno --threads "\${THREADS}" \\
            --max-ram-gb ${params.max_ram_gb ?: 0} \\
            --allele-freq-tsv "\${out}" >/dev/null 2>"\${cc}/warnings.log"; then
            [ -s "\${out}" ] || { echo "[${country}] WARNING: empty AF" >&2; rm -f "\${out}"; }
        else
            echo "[${country}] WARNING: bvs fast-allele-freq failed" >&2
            rm -f "\${out}"
        fi
    else
        echo "[${country}] WARNING: 0 usable genotype files" >&2
    fi

    # No intermediate .bvlr to clean up; just drop the staged symlinks.
    rm -rf geno
    """
}

process population_fst_aims {
    container 'ghcr.io/madhavajay/biovault-popgen:0.2.2-fast'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path country_dirs

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

    # Gather every country's outputs into one af_out/ for the FST/AIMs step.
    mkdir -p af_out
    : > successful_countries.txt
    printf 'participant_id\\tcountry\\n' > af_out/country_map.tsv
    printf 'participant_id\\tfile\\tseverity\\tcode\\tmessage\\n' > af_out/errors.tsv
    : > all_warnings.log

    for d in cc_*; do
        [ -d "\${d}" ] || continue
        for af in "\${d}"/allele_freq_*.tsv; do
            [ -e "\${af}" ] || continue
            cp "\${af}" af_out/
            base="\$(basename "\${af}" .tsv)"
            printf '%s\\n' "\${base#allele_freq_}" >> successful_countries.txt
        done
        [ -s "\${d}/country_map.part.tsv" ] && cat "\${d}/country_map.part.tsv" >> af_out/country_map.tsv || true
        if [ -s "\${d}/skipped.tsv" ]; then
            awk -F '\\t' 'BEGIN { OFS="\\t" } { print \$1, \$3, "ERROR", \$4, "country=" \$2 }' "\${d}/skipped.tsv" >> af_out/errors.tsv
        fi
        [ -s "\${d}/warnings.log" ] && cat "\${d}/warnings.log" >> all_warnings.log || true
    done

    [ -s successful_countries.txt ] || { echo "ERROR: no countries produced allele-frequency files" >&2; exit 1; }

    if [ -s all_warnings.log ]; then
        awk 'BEGIN { OFS="\\t"; print "file", "line_no", "severity", "code", "message", "raw_line" }
            /^WARNING: / {
                raw=\$0
                msg=\$0
                sub(/^WARNING: /, "", msg)
                n=split(msg, parts, ":")
                file=parts[1]
                line_no=parts[2]
                detail=msg
                if (n >= 3) {
                    prefix=file ":" line_no ": "
                    detail=substr(msg, length(prefix) + 1)
                }
                split(detail, tokens, " ")
                code=tokens[1]
                print file, line_no, "WARNING", code, detail, raw
            }' all_warnings.log > af_out/warnings.tsv
    else
        printf 'file\\tline_no\\tseverity\\tcode\\tmessage\\traw_line\\n' > af_out/warnings.tsv
    fi

    POPULATIONS="\$(sort -u successful_countries.txt | paste -sd, -)"
    echo "Populations: \${POPULATIONS}"

    bash /opt/biovault/scripts/population_level/run_pipeline.sh \\
        "\${PWD}/af_out" \\
        "\${PWD}/work" \\
        "\${PWD}" \\
        "\${POPULATIONS}"

    cp af_out/allele_freq_*.tsv .
    [ -f af_out/errors.tsv ] && cp af_out/errors.tsv errors.tsv || true
    [ -f af_out/warnings.tsv ] && cp af_out/warnings.tsv warnings.tsv || true
    """
}
