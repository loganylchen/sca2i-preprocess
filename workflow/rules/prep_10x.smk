"""10x preprocessing — barcode→group split via sinto, then UMI dedup."""


rule build_cell_groups:
    """
    Write a sinto-compatible 2-col TSV mapping cell_barcode -> group_name.
    Group = celltype (within a sample). Donor is implicit (1 sample = 1 donor
    for 10x typically; multi-donor 10x should be pre-split per sample).
    """
    input:
        cells_tsv = config["cells_tsv"],
    output:
        groups_tsv = "results/prep/10x/{sample}/groups.tsv",
    log:
        "logs/prep_10x/{sample}/build_cell_groups.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/prep_groups.py"


checkpoint split_bam_by_group:
    """
    Use sinto to split the cellranger BAM into one BAM per group in a single
    pass. Outputs to a directory; downstream rules trigger the checkpoint
    via `checkpoints.split_bam_by_group.get(sample=...).output[0]` to discover
    the actual per-group BAM filenames at runtime.
    """
    input:
        bam        = lambda wc: SAMPLES.at[wc.sample, "bam"],
        groups_tsv = "results/prep/10x/{sample}/groups.tsv",
    output:
        outdir = directory("results/prep/10x/{sample}/raw"),
    threads: config["resources"]["sinto"]["threads"]
    resources:
        mem_mb = config["resources"]["sinto"]["mem_mb"],
    log:
        "logs/prep_10x/{sample}/split_bam.log",
    conda:
        "../envs/sinto.yaml"
    shell:
        r"""
        mkdir -p {output.outdir}
        cd {output.outdir}
        sinto filterbarcodes \
            -b $OLDPWD/{input.bam} \
            -c $OLDPWD/{input.groups_tsv} \
            -p {threads} \
            > $OLDPWD/{log} 2>&1
        """


def _group_bam(wildcards):
    """Resolve the per-group BAM via the checkpoint output directory."""
    outdir = checkpoints.split_bam_by_group.get(sample=wildcards.sample).output.outdir
    return f"{outdir}/{wildcards.group}.bam"


rule index_group_bam:
    input:
        bam = _group_bam,
    output:
        bai = "results/prep/10x/{sample}/raw/{group}.bam.bai",
    log:
        "logs/prep_10x/{sample}/index_{group}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input.bam} 2> {log}"


rule umi_dedup:
    """UMI-aware deduplication per group (umi_tools dedup)."""
    input:
        bam = _group_bam,
        bai = "results/prep/10x/{sample}/raw/{group}.bam.bai",
    output:
        bam = "results/prep/10x/{sample}/dedup/{group}.bam",
        log = "results/prep/10x/{sample}/dedup/{group}.umi.log",
    threads: config["resources"]["umi_dedup"]["threads"]
    resources:
        mem_mb = config["resources"]["umi_dedup"]["mem_mb"],
    log:
        "logs/prep_10x/{sample}/umi_dedup_{group}.log",
    conda:
        "../envs/umi_tools.yaml"
    shell:
        r"""
        umi_tools dedup \
            --extract-umi-method=tag \
            --umi-tag=UB \
            --cell-tag=CB \
            --per-cell \
            --stdin={input.bam} \
            --stdout={output.bam} \
            --log={output.log} 2> {log}
        """


rule index_dedup_bam:
    input:
        "results/prep/10x/{sample}/dedup/{group}.bam",
    output:
        "results/prep/10x/{sample}/dedup/{group}.bam.bai",
    log:
        "logs/prep_10x/{sample}/index_dedup_{group}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} 2> {log}"
