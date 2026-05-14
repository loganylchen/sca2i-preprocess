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
        (curl -fsSL {params.url:q} \
            | (gunzip -c 2>/dev/null || cat) > {output.fa}.tmp \
            && samtools faidx {output.fa}.tmp \
            && mv {output.fa}.tmp {output.fa} \
            && rm -f {output.fa}.tmp.fai) 2> {log}
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
        (curl -fsSL {params.url:q} -o {output}.tmp \
            && mv {output}.tmp {output}) 2> {log}
        """


rule make_alu_bed:
    """Extract Alu* family entries from RepeatMasker .out into a 3-col BED.

    Normalises the chrom prefix to match the indexed genome (UCSC ``chr1`` vs
    Ensembl ``1``) so downstream Alu overlap isn't silently empty.
    """
    input:
        rmsk = "resources/alu/rmsk.out.gz",
        fai  = "resources/genome/genome.fa.fai",
    output:
        "resources/alu/alu.bed",
    log:
        "logs/refs/make_alu_bed.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        if [ ! -s {input.fai} ]; then
            echo "ERROR: empty or missing fai: {input.fai}" >&2
            exit 1
        fi
        HAS_CHR=$(head -1 {input.fai} | awk '{{print ($1 ~ /^chr/) ? "1" : "0"}}')
        (gunzip -c {input.rmsk} \
            | awk -v has_chr="$HAS_CHR" 'NR>3 && $11 ~ /^Alu/ {{
                chrom=$5;
                if (has_chr=="0") sub(/^chr/, "", chrom);
                else if (chrom !~ /^chr/) chrom="chr"chrom;
                print chrom"\t"($6-1)"\t"$7
              }}' \
            | sort -k1,1 -k2,2n \
            > {output}) 2> {log}
        """
