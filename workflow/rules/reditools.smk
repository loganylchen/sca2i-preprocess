"""REDItools2 de novo invocations + TSV→parquet conversion."""


_STRANDED_CHEMS = {"10x", "10x_v3", "smartseq3", "smart_seq3"}
_UNSTRANDED_CHEMS = {"smartseq2", "smart_seq2"}


def _chem_strand(chem: str) -> int:
    """Map chemistry name to REDItools -s value."""
    c = str(chem).lower()
    if c in _STRANDED_CHEMS:
        return 2  # FR-second-strand (10x v3, Smart-seq3)
    if c in _UNSTRANDED_CHEMS:
        return 0  # Smart-seq2 unstranded
    return int(config["reditools"]["strand"])


def _reditools_flags(wildcards):
    rt = config["reditools"]
    chem = SAMPLES.at[wildcards.sample, "chemistry"]
    s = _chem_strand(chem)
    return (
        f"-q {rt['min_read_quality']} "
        f"-bq {rt['min_base_quality']} "
        f"-mrl {rt['min_read_length']} "
        f"-mn {rt['min_coverage']} "
        f"-men {rt['min_edits_per_position']} "
        f"-s {s}"
    )


def _chem_for(wildcards):
    return SAMPLES.at[wildcards.sample, "chemistry"]


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
        flags = _reditools_flags,
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
        # Slice BAM to one chromosome with MAPQ>=20 (parent NN#3 lock)
        samtools view -b -@ {threads} -q 20 {input.bam:q} {wildcards.chrom} > $TMP/in.bam
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
        flags = _reditools_flags,
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
    params:
        chemistry = _chem_for,
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
    params:
        chemistry = _chem_for,
    log:
        "logs/reditools/parquet/ss/{sample}_{cell}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/parse_reditools.py"
