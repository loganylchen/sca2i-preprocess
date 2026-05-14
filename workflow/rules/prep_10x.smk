"""10x preprocessing — per-barcode split via sinto, then UMI dedup.

Per-cell architecture: each cell barcode produces one BAM. cells.tsv must
contain one row per cell barcode for the sample (obs_id = sample__barcode,
cell_barcode = the 10x barcode).

ChildIOException fix: split_bam_per_cell outputs a sentinel file (not a
``directory()``) so downstream rules can write per-cell BAI/dedup files into
the same ``split/`` directory.
"""


rule build_cell_groups:
    """
    Write a sinto-compatible 2-col TSV mapping cell_barcode -> group_name.
    For per-cell mode, group_name == cell_barcode (one BAM per cell).
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


checkpoint split_bam_per_cell:
    """
    Use sinto to split the cellranger BAM into one BAM per cell barcode in a
    single pass. Sentinel file (not directory) keeps the path tree compatible
    with downstream rules that write into the same dir (BAI, dedup outputs).
    """
    input:
        bam        = lambda wc: SAMPLES.at[wc.sample, "bam"],
        groups_tsv = "results/prep/10x/{sample}/groups.tsv",
    output:
        sentinel = "results/prep/10x/{sample}/split/.sentinel",
    threads: config["resources"]["sinto"]["threads"]
    resources:
        mem_mb = config["resources"]["sinto"]["mem_mb"],
    log:
        "logs/prep_10x/{sample}/split_bam.log",
    conda:
        "../envs/sinto.yaml"
    shell:
        r"""
        set -euo pipefail
        BAM=$(realpath {input.bam:q})
        GROUPS=$(realpath {input.groups_tsv:q})
        LOG=$(realpath {log:q})
        OUTDIR=$(dirname {output.sentinel:q})
        mkdir -p "$OUTDIR"
        cd "$OUTDIR"
        sinto filterbarcodes \
            -b "$BAM" \
            -c "$GROUPS" \
            -p {threads} \
            > "$LOG" 2>&1
        touch .sentinel
        """


def _split_bam(wildcards):
    """Resolve per-cell BAM, gated by the split checkpoint."""
    checkpoints.split_bam_per_cell.get(sample=wildcards.sample)
    return f"results/prep/10x/{wildcards.sample}/split/{wildcards.group}.bam"


rule index_cell_bam:
    input:
        bam = _split_bam,
    output:
        bai = "results/prep/10x/{sample}/split/{group}.bam.bai",
    log:
        "logs/prep_10x/{sample}/index_{group}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input.bam} 2> {log}"


rule umi_dedup:
    """UMI-aware deduplication per cell (umi_tools dedup)."""
    input:
        bam = _split_bam,
        bai = "results/prep/10x/{sample}/split/{group}.bam.bai",
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
