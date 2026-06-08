// BioVault popgen: CookHLA HLA imputation from cohort PLINK BED/BIM/FAM.

nextflow.enable.dsl=2

if (!params.containsKey('biosynth_image')) {
    params.biosynth_image = null
}
if (!params.containsKey('cookhla_image')) {
    params.cookhla_image = null
}
if (!params.containsKey('cookhla_input_build')) {
    params.cookhla_input_build = '38'
}
if (!params.containsKey('cookhla_panel')) {
    params.cookhla_panel = 'ALL'
}
if (!params.containsKey('cookhla_reference_panel')) {
    params.cookhla_reference_panel = null
}
if (!params.containsKey('cookhla_genetic_map')) {
    params.cookhla_genetic_map = null
}
if (!params.containsKey('cookhla_average_erate')) {
    params.cookhla_average_erate = null
}
if (!params.containsKey('cookhla_mem')) {
    params.cookhla_mem = '4g'
}
if (!params.containsKey('cookhla_threads')) {
    params.cookhla_threads = '4'
}
if (!params.containsKey('cookhla_multiprocess')) {
    params.cookhla_multiprocess = '4'
}

def BIOSYNTH_IMAGE = System.getenv('BIOSYNTH_IMAGE') ?: (params.biosynth_image ?: 'ghcr.io/openmined/biosynth:0.1.32')
def COOKHLA_IMAGE = params.cookhla_image ?: 'ghcr.io/madhavajay/cookhla-rs:latest'
def EMPTY_HLA_OUTPUTS = '''
            printf 'participant_id\\tcountry\\tfamily_id\\tindividual_id\\thla_gene\\tbroad_2_digit_call\\tspecific_4_digit_call\\tallele_1\\tallele_2\\tgenotype_4digit\\tallele_1_posterior\\tallele_2_posterior\\tcombined_posterior\\treference_panel\\tinput_build\\toutput_prefix\\n' > hla_individual_results.tsv
            printf 'country\\thla_gene\\tgenotype_4digit\\tcount\\tsample_count\\tfrequency\\n' > hla_country_genotype_counts.tsv
            printf 'country\\thla_gene\\tsample_count\\tmean_combined_posterior\\n' > hla_country_gene_summary.tsv
'''

def countryValue(record) {
    def facets = record.facets ?: [:]
    return (
        record.country ?:
        record.Country ?:
        facets.country ?:
        facets.Country ?:
        'Unknown'
    ).toString().trim() ?: 'Unknown'
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
            def country = countryValue(record)
            if (country == 'Unknown') {
                println "[bv] WARNING: participant ${record.participant_id}: missing country facet; using Unknown"
            }
            return [tuple(record.participant_id.toString(), country, file(record.genotype_file))]
        }

        def collected = records
            .collect(flat: false)
            .map { items ->
                if (items.isEmpty()) {
                    throw new IllegalArgumentException("No valid participants with genotype files remained")
                }
                tuple(
                    items.collect { it[0] },
                    items.collect { it[1] },
                    items.collect { it[2] }
                )
            }

        def bed = cohort_bed_for_cookhla(collected)
        def hla = run_cookhla(bed.bed_dir, bed.country_map)

    emit:
        hla_individual_results = hla.individual_results
        hla_country_genotype_counts = hla.country_genotype_counts
        hla_country_gene_summary = hla.country_gene_summary
        errors = hla.errors
        cookhla_log = hla.cookhla_log
        cookhla_raw_archive = hla.raw_archive
}

process cohort_bed_for_cookhla {
    container BIOSYNTH_IMAGE
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), val(countries), path(genotype_files)

    output:
        path "plink_bed", emit: bed_dir
        path "country_map.tsv", emit: country_map

    script:
    def staging = []
    def countryRows = []
    participant_ids.eachWithIndex { pid, idx ->
        def fname = genotype_files[idx].getName()
        staging << "mkdir -p input/${pid} && ln -s \"../../${fname}\" \"input/${pid}/${fname}\""
        countryRows << "${pid}\t${countries[idx].replace('\t', ' ').replace('\n', ' ')}"
    }
    def countryText = countryRows.join('\\n')
    """
    set -euo pipefail
    mkdir -p input plink_bed
    ${staging.join('\n    ')}
    { printf 'participant_id\\tcountry\\n'; printf '%b\\n' "${countryText}"; } \\
        > country_map.tsv
    bvs cohort-bed -i input \\
        --out-prefix plink_bed/genotypes \\
        --snp-info plink_bed/snp_info.tsv
    """
}

process run_cookhla {
    container COOKHLA_IMAGE
    containerOptions '--entrypoint=""'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path bed_dir
        path country_map

    output:
        path "hla_individual_results.tsv", emit: individual_results
        path "hla_country_genotype_counts.tsv", emit: country_genotype_counts
        path "hla_country_gene_summary.tsv", emit: country_gene_summary
        path "errors.tsv", emit: errors
        path "cookhla.log", emit: cookhla_log
        path "cookhla_raw.tar.gz", emit: raw_archive

    script:
    """
    set -euo pipefail

    export PATH=/opt/conda/bin:/opt/conda/condabin:/usr/local/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin:\${PATH}
    WORKDIR="\${PWD}"
    INPUT_PREFIX="\${WORKDIR}/${bed_dir}/genotypes"
    OUTDIR="\${WORKDIR}/cookhla_raw"
    OUT_PREFIX="\${OUTDIR}/cookhla"
    MHC_DIR="\${WORKDIR}/cookhla_mhc_plink"
    MHC_INPUT_PREFIX="\${MHC_DIR}/genotypes_mhc"
    mkdir -p "\${OUTDIR}"
    mkdir -p "\${MHC_DIR}"

    PANEL="${params.cookhla_panel}"
    case "\${PANEL}" in
        ALL|AFR|AMR|EAS|EUR|SAS)
            REFERENCE_PANEL="/opt/cookhla/1000G_REF/1000G_REF.\${PANEL}.chr6.hg18.29mb-34mb.inT1DGC"
            ;;
        CEU|HM_CEU_REF|example/HM_CEU_REF)
            REFERENCE_PANEL="/opt/cookhla/example/HM_CEU_REF"
            ;;
        *)
            {
                printf 'participant_id\\tseverity\\tcode\\tmessage\\n'
                printf 'cohort\\tERROR\\tinvalid_panel\\tInvalid CookHLA panel "%s"; expected one of ALL, AFR, AMR, EAS, EUR, SAS, CEU\\n' "\${PANEL}"
            } > errors.tsv
${EMPTY_HLA_OUTPUTS}
            printf 'Invalid CookHLA panel: %s\\n' "\${PANEL}" > cookhla.log
            tar -czf cookhla_raw.tar.gz cookhla_raw
            exit 0
            ;;
    esac

    PANEL_OVERRIDE="${params.cookhla_reference_panel}"
    if [ -n "\${PANEL_OVERRIDE}" ] && [ "\${PANEL_OVERRIDE}" != "null" ]; then
        REFERENCE_PANEL="${params.cookhla_reference_panel}"
    fi

    GENETIC_MAP="${params.cookhla_genetic_map}"
    AVERAGE_ERATE="${params.cookhla_average_erate}"
    if { [ -n "\${GENETIC_MAP}" ] && [ "\${GENETIC_MAP}" != "null" ]; } || { [ -n "\${AVERAGE_ERATE}" ] && [ "\${AVERAGE_ERATE}" != "null" ]; }; then
        if [ -z "\${GENETIC_MAP}" ] || [ "\${GENETIC_MAP}" = "null" ] || [ -z "\${AVERAGE_ERATE}" ] || [ "\${AVERAGE_ERATE}" = "null" ]; then
            {
                printf 'participant_id\\tseverity\\tcode\\tmessage\\n'
                printf 'cohort\\tERROR\\tpartial_map_override\\tCookHLA requires both cookhla_genetic_map and cookhla_average_erate, or neither so it can generate an adaptive map\\n'
            } > errors.tsv
${EMPTY_HLA_OUTPUTS}
            printf 'Partial CookHLA map override\\n' > cookhla.log
            tar -czf cookhla_raw.tar.gz cookhla_raw
            exit 0
        fi
        MAP_ARGS="--genetic-map \${GENETIC_MAP} --average-erate \${AVERAGE_ERATE}"
    elif [ "\${REFERENCE_PANEL}" = "/opt/cookhla/example/HM_CEU_REF" ]; then
        MAP_ARGS="--genetic-map /opt/cookhla/example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB --average-erate /opt/cookhla/example/AGM.1958BC+HM_CEU_REF.aver.erate"
    else
        {
            printf 'participant_id\\tseverity\\tcode\\tmessage\\n'
            printf 'cohort\\tERROR\\tcookhla_rs_missing_adaptive_map\\tcookhla-rs cannot auto-generate adaptive genetic maps yet; panel %s requires cookhla_genetic_map and cookhla_average_erate, or use panel CEU for the bundled example map\\n' "\${PANEL}"
        } > errors.tsv
${EMPTY_HLA_OUTPUTS}
        printf 'cookhla-rs requires a precomputed adaptive genetic map for panel %s; bundled map exists only for CEU/example/HM_CEU_REF\\n' "\${PANEL}" > cookhla.log
        tar -czf cookhla_raw.tar.gz cookhla_raw
        exit 0
    fi

    printf 'participant_id\\tseverity\\tcode\\tmessage\\n' > errors.tsv

    /opt/conda/bin/plink \\
        --bfile "\${INPUT_PREFIX}" \\
        --chr 6 \\
        --from-bp 29000000 \\
        --to-bp 34000000 \\
        --allow-extra-chr \\
        --make-bed \\
        --out "\${MHC_INPUT_PREFIX}" \\
        >> "\${WORKDIR}/cookhla.log" 2>&1

    if [ ! -s "\${MHC_INPUT_PREFIX}.bim" ]; then
        {
            printf 'participant_id\\tseverity\\tcode\\tmessage\\n'
            printf 'cohort\\tERROR\\tno_mhc_variants\\tNo chr6:29000000-34000000 variants remained after PLINK filtering\\n'
        } > errors.tsv
        printf 'participant_id\\tcountry\\tfamily_id\\tindividual_id\\thla_gene\\tbroad_2_digit_call\\tspecific_4_digit_call\\tallele_1\\tallele_2\\tgenotype_4digit\\tallele_1_posterior\\tallele_2_posterior\\tcombined_posterior\\treference_panel\\tinput_build\\toutput_prefix\\n' > hla_individual_results.tsv
        printf 'country\\thla_gene\\tgenotype_4digit\\tcount\\tsample_count\\tfrequency\\n' > hla_country_genotype_counts.tsv
        printf 'country\\thla_gene\\tsample_count\\tmean_combined_posterior\\n' > hla_country_gene_summary.tsv
        tar -czf cookhla_raw.tar.gz cookhla_raw cookhla_mhc_plink
        exit 0
    fi

    /usr/local/bin/cookhla \\
        --input "\${MHC_INPUT_PREFIX}" \\
        --human-genome "${params.cookhla_input_build}" \\
        --out "\${OUT_PREFIX}" \\
        --reference "\${REFERENCE_PANEL}" \\
        \${MAP_ARGS} \\
        --java-memory "${params.cookhla_mem}" \\
        --nthreads "${params.cookhla_threads}" \\
        --multiprocess "${params.cookhla_multiprocess}" \\
        >> "\${WORKDIR}/cookhla.log" 2>&1

    ALLELES="\${OUT_PREFIX}.MHC.HLA_IMPUTATION_OUT.alleles"
    if [ ! -s "\${ALLELES}" ]; then
        echo "CookHLA did not produce a non-empty alleles file: \${ALLELES}" >&2
        tail -200 cookhla.log >&2 || true
        exit 1
    fi

    awk -F '\\t' 'BEGIN{OFS="\\t"} NR>1{country[\$1]=\$2} END{for (id in country) print id,country[id]}' "${country_map}" > country_map.clean.tsv

    awk -v reference_panel="\${REFERENCE_PANEL}" \\
        -v input_build="hg${params.cookhla_input_build}" \\
        -v output_prefix="\${OUT_PREFIX}" \\
        -v country_file="country_map.clean.tsv" '
        BEGIN {
            FS="[ \\t]+"
            OFS="\\t"
            while ((getline line < country_file) > 0) {
                split(line, c, "\\t")
                country[c[1]] = c[2]
            }
            close(country_file)
            print "participant_id","country","family_id","individual_id","hla_gene","broad_2_digit_call","specific_4_digit_call","allele_1","allele_2","genotype_4digit","allele_1_posterior","allele_2_posterior","combined_posterior","reference_panel","input_build","output_prefix"
        }
        function fmt(g, a) {
            if (a == "" || a == "0") return "0"
            return "HLA-" g "*" substr(a, 1, length(a) - 2) ":" substr(a, length(a) - 1, 2)
        }
        {
            split(\$5, a, ",")
            gene = "HLA-" \$3
            allele1 = fmt(\$3, a[1])
            allele2 = fmt(\$3, a[2])
            genotype = allele1 "/" allele2
            participant = \$2
            c = country[participant]
            if (c == "") c = country[\$1]
            if (c == "") c = "Unknown"
            print participant,c,\$1,\$2,gene,\$4,\$5,allele1,allele2,genotype,\$6,\$7,\$8,reference_panel,input_build,output_prefix
        }
    ' "\${ALLELES}" > hla_individual_results.tsv

    awk -F '\\t' 'BEGIN{OFS="\\t"}
        NR == 1 { next }
        {
            key = \$2 SUBSEP \$5 SUBSEP \$10
            count[key]++
            sample_key = \$2 SUBSEP \$5 SUBSEP \$1
            if (!(sample_key in seen)) {
                seen[sample_key] = 1
                denom[\$2 SUBSEP \$5]++
            }
        }
        END {
            print "country","hla_gene","genotype_4digit","count","sample_count","frequency"
            for (key in count) {
                split(key, parts, SUBSEP)
                d = denom[parts[1] SUBSEP parts[2]]
                freq = (d > 0) ? count[key] / d : 0
                printf "%s\\t%s\\t%s\\t%d\\t%d\\t%.6f\\n", parts[1], parts[2], parts[3], count[key], d, freq
            }
        }
    ' hla_individual_results.tsv > hla_country_genotype_counts.unsorted.tsv
    { head -n 1 hla_country_genotype_counts.unsorted.tsv; tail -n +2 hla_country_genotype_counts.unsorted.tsv | sort -t \$'\\t' -k1,1 -k2,2 -k3,3; } > hla_country_genotype_counts.tsv

    awk -F '\\t' 'BEGIN{OFS="\\t"}
        NR == 1 { next }
        {
            key = \$2 SUBSEP \$5
            sample_key = \$2 SUBSEP \$5 SUBSEP \$1
            if (!(sample_key in seen)) {
                seen[sample_key] = 1
                sample_count[key]++
            }
            posterior_sum[key] += \$13
            posterior_n[key]++
        }
        END {
            print "country","hla_gene","sample_count","mean_combined_posterior"
            for (key in sample_count) {
                split(key, parts, SUBSEP)
                mean = (posterior_n[key] > 0) ? posterior_sum[key] / posterior_n[key] : 0
                printf "%s\\t%s\\t%d\\t%.6f\\n", parts[1], parts[2], sample_count[key], mean
            }
        }
    ' hla_individual_results.tsv > hla_country_gene_summary.unsorted.tsv
    { head -n 1 hla_country_gene_summary.unsorted.tsv; tail -n +2 hla_country_gene_summary.unsorted.tsv | sort -t \$'\\t' -k1,1 -k2,2; } > hla_country_gene_summary.tsv

    tar -czf cookhla_raw.tar.gz cookhla_raw cookhla_mhc_plink
    """
}
