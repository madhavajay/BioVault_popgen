#!/usr/bin/env python
"""
Dump the variant list from a gnomAD PCA loadings Hail Table to a plain TSV.

Why: with only 10 study samples we can't LD-prune locally (PLINK needs ≥50).
But the loadings HT IS already an LD-pruned, ancestry-informative SNP set
computed by gnomAD on the full HGDP+1kGP panel. Using its variant list as
the study-side filter gives the same effect.

Output (TSV, no header):
    chrom_no_prefix    pos    ref    alt    id_chrpos
    1                  12345  A      G      1:12345

Usage:
    python extract_loadings_variants.py <loadings_ht_path> <out_tsv>
"""

import sys
import os
from pathlib import Path

import hail as hl


def main():
    if len(sys.argv) < 3:
        print("Usage: extract_loadings_variants.py <loadings_ht> <out_tsv>")
        sys.exit(1)

    ht_path = sys.argv[1]
    out_tsv = Path(sys.argv[2])
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    spark_conf = {
        "spark.driver.memory": "4g",
    }
    gcs_jar = os.environ.get("GCS_CONNECTOR_JAR")
    if gcs_jar:
        spark_conf.update({
            "spark.jars": gcs_jar,
            "spark.driver.extraClassPath": gcs_jar,
            "spark.executor.extraClassPath": gcs_jar,
            "spark.hadoop.fs.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem",
            "spark.hadoop.fs.AbstractFileSystem.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS",
            "spark.hadoop.fs.gs.auth.type": "UNAUTHENTICATED",
            "spark.hadoop.fs.gs.auth.service.account.enable": "false",
            "spark.hadoop.fs.gs.auth.null.enable": "true",
        })

    hl.init(default_reference="GRCh38", quiet=True,
            log=str(out_tsv.parent / "hail_extract.log"),
            spark_conf=spark_conf)

    print(f"Reading {ht_path}")
    ht = hl.read_table(ht_path)

    # locus.contig is "chr1" / "chr2" / ... in GRCh38; strip prefix to match
    # plink2 --set-all-var-ids '@:#' output (which uses bare chrom).
    ht = ht.annotate(
        chrom=ht.locus.contig.replace("^chr", ""),
        pos=ht.locus.position,
        ref=ht.alleles[0],
        alt=ht.alleles[1],
    )
    ht = ht.annotate(id_chrpos=ht.chrom + ":" + hl.str(ht.pos))

    # Drop the original keys (locus, alleles) so export gives clean 5-col TSV.
    out_fields = ht.key_by().select("chrom", "pos", "ref", "alt", "id_chrpos")

    # Collapse 3709 source partitions into a single Spark task. The HT is small
    # (~30k rows total) but its partitioning is fine-grained, and per-task
    # scheduling overhead dominates wall time when the JVM runs under
    # Rosetta/qemu emulation (linux/amd64 image on an Apple Silicon host).
    out_fields = out_fields.naive_coalesce(1)
    out_fields.export(str(out_tsv), header=False)

    with open(out_tsv) as f:
        n = sum(1 for _ in f)
    print(f"Wrote {n} loadings variants -> {out_tsv}")


if __name__ == "__main__":
    main()
