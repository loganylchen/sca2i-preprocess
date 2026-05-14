"""Smart-seq2/3 preprocessing — Picard MarkDuplicates per cell BAM."""


def _ss_input_bam(wildcards):
    cell_dir = SAMPLES.at[wildcards.sample, "cell_bam_dir"].rstrip("/")
    return f"{cell_dir}/{wildcards.cell}.bam"


rule sort_ss:
    input:
        bam = _ss_input_bam,
    output:
        bam = temp("results/prep/ss/{sample}/sorted/{cell}.bam"),
    threads: config["resources"]["markdup"]["threads"]
    resources:
        mem_mb = config["resources"]["markdup"]["mem_mb"],
    log:
        "logs/prep_ss/{sample}/sort_{cell}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} 2> {log}"


rule markdup_ss:
    input:
        bam = "results/prep/ss/{sample}/sorted/{cell}.bam",
    output:
        bam     = "results/prep/ss/{sample}/dedup/{cell}.bam",
        metrics = "results/prep/ss/{sample}/dedup/{cell}.metrics.txt",
    threads: config["resources"]["markdup"]["threads"]
    resources:
        mem_mb = config["resources"]["markdup"]["mem_mb"],
    log:
        "logs/prep_ss/{sample}/markdup_{cell}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        picard MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORT_ORDER=coordinate 2> {log}
        """


rule index_ss_dedup:
    input:
        "results/prep/ss/{sample}/dedup/{cell}.bam",
    output:
        "results/prep/ss/{sample}/dedup/{cell}.bam.bai",
    log:
        "logs/prep_ss/{sample}/index_{cell}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} 2> {log}"
