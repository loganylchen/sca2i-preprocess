"""REDItools2 de novo invocations + TSV→parquet conversion."""


def _reditools_flags():
    rt = config["reditools"]
    return (
        f"-q {rt['min_read_quality']} "
        f"-bq {rt['min_base_quality']} "
        f"-mrl {rt['min_read_length']} "
        f"-mn {rt['min_coverage']} "
        f"-men {rt['min_edits_per_position']} "
        f"-s {rt['strand']}"
    )


rule reditools_10x:
    """
    Per-group × per-chrom REDItools2 de novo. Slices the BAM with samtools view
    so REDItools sees only one chromosome at a time (parallelism + bounded mem).
    """
    input:
        bam = "results/prep/10x/{sample}/dedup/{group}.bam",
        bai = "results/prep/10x/{sample}/dedup/{group}.bam.bai",
        fa  = "resources/genome/genome.fa",
        fai = "resources/genome/genome.fa.fai",
    output:
        tsv = "results/reditools/10x/{sample}/{group}/{chrom}.tsv.gz",
    params:
        flags = _reditools_flags(),
    threads: config["resources"]["reditools"]["threads"]
    resources:
        mem_mb = config["resources"]["reditools"]["mem_mb"],
    log:
        "logs/reditools/10x/{sample}/{group}_{chrom}.log",
    conda:
        "../envs/reditools.yaml"
    shell:
        r"""
        set -euo pipefail
        TMP=$(mktemp -d)
        trap "rm -rf $TMP" EXIT
        # Slice BAM to one chromosome
        samtools view -b -@ {threads} {input.bam} {wildcards.chrom} > $TMP/in.bam
        samtools index $TMP/in.bam
        reditools.py \
            -f $TMP/in.bam \
            -r {input.fa} \
            -t {threads} \
            -g {wildcards.chrom} \
            {params.flags} \
            -o $TMP/out.tsv 2> {log}
        gzip -c $TMP/out.tsv > {output.tsv}
        """


rule reditools_ss:
    """Per-cell REDItools2 de novo (whole BAM, all chroms in one pass)."""
    input:
        bam = "results/prep/ss/{sample}/dedup/{cell}.bam",
        bai = "results/prep/ss/{sample}/dedup/{cell}.bam.bai",
        fa  = "resources/genome/genome.fa",
        fai = "resources/genome/genome.fa.fai",
    output:
        tsv = "results/reditools/ss/{sample}/{cell}.tsv.gz",
    params:
        flags = _reditools_flags(),
    threads: config["resources"]["reditools"]["threads"]
    resources:
        mem_mb = config["resources"]["reditools"]["mem_mb"],
    log:
        "logs/reditools/ss/{sample}/{cell}.log",
    conda:
        "../envs/reditools.yaml"
    shell:
        r"""
        set -euo pipefail
        TMP=$(mktemp -d)
        trap "rm -rf $TMP" EXIT
        reditools.py \
            -f {input.bam} \
            -r {input.fa} \
            -t {threads} \
            {params.flags} \
            -o $TMP/out.tsv 2> {log}
        gzip -c $TMP/out.tsv > {output.tsv}
        """


rule reditools_to_parquet_10x:
    input:
        tsv = "results/reditools/10x/{sample}/{group}/{chrom}.tsv.gz",
    output:
        parquet = "results/reditools/10x/{sample}/{group}/{chrom}.parquet",
    log:
        "logs/reditools/parquet/10x/{sample}_{group}_{chrom}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/parse_reditools.py"


rule reditools_to_parquet_ss:
    input:
        tsv = "results/reditools/ss/{sample}/{cell}.tsv.gz",
    output:
        parquet = "results/reditools/ss/{sample}/{cell}.parquet",
    log:
        "logs/reditools/parquet/ss/{sample}_{cell}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/parse_reditools.py"
