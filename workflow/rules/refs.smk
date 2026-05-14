"""Reference resources — genome FASTA + Alu BED."""


rule download_genome:
    output:
        fa = "resources/genome/genome.fa",
    params:
        url = config["genome"]["fasta_url"],
    log:
        "logs/refs/download_genome.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        mkdir -p resources/genome
        (curl -fsSL "{params.url}" \
            | (gunzip -c 2>/dev/null || cat) > {output.fa}) 2> {log}
        """


rule index_genome:
    input:
        "resources/genome/genome.fa",
    output:
        "resources/genome/genome.fa.fai",
    log:
        "logs/refs/index_genome.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} 2> {log}"


rule download_alu_raw:
    output:
        "resources/alu/rmsk.out.gz",
    params:
        url = config["alu"]["rmsk_url"],
    log:
        "logs/refs/download_alu.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        mkdir -p resources/alu
        curl -fsSL "{params.url}" -o {output} 2> {log}
        """


rule make_alu_bed:
    """Extract Alu* family entries from RepeatMasker .out into a 3-col BED."""
    input:
        "resources/alu/rmsk.out.gz",
    output:
        "resources/alu/alu.bed",
    log:
        "logs/refs/make_alu_bed.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        (gunzip -c {input} \
            | awk 'NR>3 && $11 ~ /^Alu/ {{ print $5"\t"($6-1)"\t"$7 }}' \
            | sort -k1,1 -k2,2n \
            > {output}) 2> {log}
        """
