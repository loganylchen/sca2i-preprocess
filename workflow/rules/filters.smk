"""ADAR motif annotation + Alu overlap annotation. Annotation only — no rows
are dropped here. Downstream sca2i.tl.discover() decides what to use."""


rule annotate_motif_10x:
    input:
        parquet = "results/reditools/10x/{sample}/{group}/{chrom}.parquet",
        fa      = "resources/genome/genome.fa",
        fai     = "resources/genome/genome.fa.fai",
    output:
        parquet = "results/filtered/{sample}/10x__{group}__{chrom}.motif.parquet",
    params:
        window  = config["filters"]["motif_window"],
        do_run  = config["filters"]["annotate_motif"],
    log:
        "logs/filters/motif/10x_{sample}_{group}_{chrom}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/filter_motif.py"


rule annotate_motif_ss:
    input:
        parquet = "results/reditools/ss/{sample}/{cell}.parquet",
        fa      = "resources/genome/genome.fa",
        fai     = "resources/genome/genome.fa.fai",
    output:
        parquet = "results/filtered/{sample}/ss__{cell}.motif.parquet",
    params:
        window  = config["filters"]["motif_window"],
        do_run  = config["filters"]["annotate_motif"],
    log:
        "logs/filters/motif/ss_{sample}_{cell}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/filter_motif.py"


rule annotate_alu_10x:
    input:
        parquet = "results/filtered/{sample}/10x__{group}__{chrom}.motif.parquet",
        bed     = "resources/alu/alu.bed",
    output:
        parquet = "results/filtered/{sample}/10x__{group}__{chrom}.annot.parquet",
    params:
        do_run = config["filters"]["annotate_alu"],
    log:
        "logs/filters/alu/10x_{sample}_{group}_{chrom}.log",
    conda:
        "../envs/pybedtools.yaml"
    script:
        "../scripts/filter_alu.py"


rule annotate_alu_ss:
    input:
        parquet = "results/filtered/{sample}/ss__{cell}.motif.parquet",
        bed     = "resources/alu/alu.bed",
    output:
        parquet = "results/filtered/{sample}/ss__{cell}.annot.parquet",
    params:
        do_run = config["filters"]["annotate_alu"],
    log:
        "logs/filters/alu/ss_{sample}_{cell}.log",
    conda:
        "../envs/pybedtools.yaml"
    script:
        "../scripts/filter_alu.py"
