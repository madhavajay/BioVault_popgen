// BioVault popgen: within-cohort QC + PCA sanity check.
//
// Two-container pipeline (like flow 04):
//   cohort_bed   (BIOSYNTH_IMAGE) — fast Rust parse of every genotype file into the
//                cohort PLINK .bed/.bim/.fam via `bvs cohort-bed`. Byte-for-byte
//                identical to fast_pipeline.py's bed, but parallel + bounded-RAM
//                (replaces the slow Python parse/memmap that used to OOM at scale).
//   pca_qc_fast  (biovault-popgen) — runs fast_pipeline.py with BV_PREBUILT_BED set,
//                so it skips parse/memmap/bed and goes straight to the same plink2
//                QC/PCA + plots. Outputs are identical (verified) because QC/PCA
//                runs on the identical bed.

nextflow.enable.dsl=2

def BIOSYNTH_IMAGE = System.getenv('BIOSYNTH_IMAGE') ?: 'ghcr.io/openmined/biosynth:0.1.32'

workflow USER {
    take:
        context
        participants

    main:
        def participantTuples = participants.flatMap { record ->
            def validation = record.validation ?: [status: 'ok', message: '']
            if (validation.status?.toString() != 'ok') {
                println "[bv] WARNING: skipping participant ${record.participant_id}: genotype file ${validation.status} - ${validation.message}"
                return []
            }
            def facets = record.facets ?: [:]
            def country = (record.country ?: facets.country ?: '(unset)').toString().trim() ?: '(unset)'
            def sex = (record.sex ?: facets.sex ?: '(unset)').toString().trim() ?: '(unset)'
            return [tuple(
                record.participant_id.toString(),
                file(record.genotype_file),
                country,
                sex
            )]
        }
        def collected = participantTuples
            .collect(flat: false)
            .map { items ->
                if (items.isEmpty()) {
                    throw new IllegalArgumentException("No valid genotype files remained after BioVault input validation")
                }
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
        def bed = cohort_bed(collected)
        def qc = pca_qc_fast(bed.bed_dir)

    emit:
        eigenvec = qc.eigenvec
        eigenval = qc.eigenval
        pca_plot_points = qc.pca_plot_points
        snp_info = qc.snp_info
        pca_prefiltered_snps = qc.pca_prefiltered_snps
        filtered_snps = qc.filtered_snps
        pca_pc12_plot = qc.pca_pc12_plot
        pca_pc34_plot = qc.pca_pc34_plot
        pipeline_log = qc.pipeline_log
        errors = qc.errors
        warnings = qc.warnings
}

// Fast cohort PLINK bed build (biosynth). Replaces fast_pipeline.py's slow
// parse + uint8 memmap + bed write with one parallel, bounded-RAM Rust pass.
process cohort_bed {
    container BIOSYNTH_IMAGE
    // biosynth image sets ENTRYPOINT ["bvs"]; clear it so Nextflow's bash runs.
    containerOptions '--entrypoint=""'
    stageInMode 'symlink'
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple val(participant_ids), path(genotype_files)

    output:
        path "plink_bed", emit: bed_dir

    script:
    def staging = []
    participant_ids.eachWithIndex { pid, idx ->
        def fname = genotype_files[idx].getName()
        staging << "mkdir -p input/${pid} && ln -s \"../../${fname}\" \"input/${pid}/${fname}\""
    }
    """
    set -euo pipefail
    mkdir -p input plink_bed
    ${staging.join('\n    ')}
    # Parse all samples -> cohort PLINK bed + full SNP-universe snp_info.tsv.
    # Byte-for-byte identical to fast_pipeline.py; no reference DB needed.
    bvs cohort-bed -i input \\
        --out-prefix plink_bed/genotypes \\
        --snp-info plink_bed/snp_info.tsv
    """
}

process pca_qc_fast {
    container 'ghcr.io/madhavajay/biovault-popgen:0.2.4-fast'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'symlink'
    errorStrategy 'terminate'
    maxRetries { params.nextflow.max_retries }

    input:
        path bed_dir

    output:
        path "pca.eigenvec",     emit: eigenvec
        path "pca.eigenval",     emit: eigenval
        path "pca_plot_points.tsv", emit: pca_plot_points
        path "snp_info.tsv",     emit: snp_info
        path "pca_prefiltered_snps.tsv", emit: pca_prefiltered_snps
        path "filtered_snps.tsv", emit: filtered_snps
        path "pca_pc1_pc2.png",  emit: pca_pc12_plot, optional: true
        path "pca_pc3_pc4.png",  emit: pca_pc34_plot, optional: true
        path "fast_pipeline.log", emit: pipeline_log, optional: true
        path "errors.tsv",       emit: errors, optional: true
        path "warnings.tsv",     emit: warnings, optional: true

    script:
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

    # fast_pipeline.py uses BASE_DIR = Path(__file__).parents[1]; reconstruct the
    # layout in the writable workdir so outputs land in <base>/{data,plots,logs}.
    mkdir -p pca_qc_fast/scripts pca_qc_fast/data/plink pca_qc_fast/data/merged pca_qc_fast/data/qc pca_qc_fast/plots pca_qc_fast/logs
    cp /opt/biovault/scripts/pca_qc_fast/*.py pca_qc_fast/scripts/

    source /opt/conda/etc/profile.d/conda.sh
    conda activate biovault_popgen

    # Prebuilt cohort bed from cohort_bed (bvs). fast_pipeline.py skips
    # parse/memmap/bed and runs the same plink2 QC/PCA + plots on it, so the
    # eigenvec is identical (verified). snp_info.tsv comes from bvs cohort-bed.
    export BV_PREBUILT_BED="\${PWD}/${bed_dir}/genotypes"
    cp "${bed_dir}/snp_info.tsv" pca_qc_fast/data/merged/snp_info.tsv 2>/dev/null || true

    set +e
    python3 pca_qc_fast/scripts/fast_pipeline.py
    status=\$?
    set -e
    if [ "\${status}" -ne 0 ]; then
        mkdir -p pca_qc_fast/logs
        reason="pca_qc_fast failed with exit \${status}"
        if [ "\${status}" -eq 137 ]; then
            reason="pca_qc_fast was killed with exit 137 (out of memory during QC/PCA). Increase Docker memory or reduce selected files."
        fi
        {
            printf '%s\\n' "\${reason}"
            printf 'work_dir\\t%s\\n' "\${PWD}"
            printf 'bed_dir\\t%s\\n' "${bed_dir}"
        } > pca_qc_fast/logs/failure_summary.txt
        if [ ! -f pca_qc_fast/logs/errors.tsv ]; then
            printf 'participant_id\\tfile\\tseverity\\tcode\\tmessage\\n' > pca_qc_fast/logs/errors.tsv
        fi
        printf 'COHORT\\t%s\\tERROR\\tPIPELINE_FAILED\\t%s\\n' "${bed_dir}" "\${reason}" >> pca_qc_fast/logs/errors.tsv
        [ -f pca_qc_fast/logs/fast_pipeline.log ] && printf '\\nERROR: %s\\n' "\${reason}" >> pca_qc_fast/logs/fast_pipeline.log
        cp pca_qc_fast/logs/failure_summary.txt failure_summary.txt
        cp pca_qc_fast/logs/errors.tsv errors.tsv
        [ -f pca_qc_fast/logs/fast_pipeline.log ] && cp pca_qc_fast/logs/fast_pipeline.log fast_pipeline.log || true
        echo "ERROR: \${reason}" >&2
        exit "\${status}"
    fi

    # Hoist final artefacts to the process root so publishDir picks them up.
    cp pca_qc_fast/data/pca/pca.eigenvec       pca.eigenvec
    cp pca_qc_fast/data/pca/pca.eigenval       pca.eigenval
    cp pca_qc_fast/data/pca/pca_plot_points.tsv pca_plot_points.tsv
    cp pca_qc_fast/data/merged/snp_info.tsv    snp_info.tsv
    cp pca_qc_fast/data/qc/pca_prefiltered_snps.tsv pca_prefiltered_snps.tsv
    cp pca_qc_fast/data/qc/filtered_snps.tsv   filtered_snps.tsv
    [ -f pca_qc_fast/plots/pca_pc1_pc2.png ] && cp pca_qc_fast/plots/pca_pc1_pc2.png pca_pc1_pc2.png || true
    [ -f pca_qc_fast/plots/pca_pc3_pc4.png ] && cp pca_qc_fast/plots/pca_pc3_pc4.png pca_pc3_pc4.png || true
    [ -f pca_qc_fast/logs/fast_pipeline.log ] && cp pca_qc_fast/logs/fast_pipeline.log fast_pipeline.log || true
    [ -f pca_qc_fast/logs/errors.tsv ] && cp pca_qc_fast/logs/errors.tsv errors.tsv || true
    [ -f pca_qc_fast/logs/warnings.tsv ] && cp pca_qc_fast/logs/warnings.tsv warnings.tsv || true
    """
}
