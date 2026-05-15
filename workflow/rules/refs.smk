"""Reference resources — genome FASTA + Alu BED + REDItools2 source."""


# Pinned commit SHA for REDItools2 reproducibility.
REDITOOLS2_COMMIT = "1e9d396f8aba058f073d7d27bf1148fe408adaf8"


rule clone_reditools2:
    """REDItools2 ships as bare scripts (no setup.py), so we clone the repo
    at a pinned commit. reditools.smk invokes the resulting reditools.py."""
    output:
        script = "resources/reditools2/src/cineca/reditools.py",
    params:
        commit = REDITOOLS2_COMMIT,
    log:
        "logs/refs/clone_reditools2.log",
    conda:
        "../envs/reditools.yaml"
    shell:
        r"""
        set -euo pipefail
        DEST=resources/reditools2
        if [ -d "$DEST/.git" ]; then
            git -C "$DEST" fetch --quiet origin {params.commit} || true
        else
            rm -rf "$DEST"
            git clone --quiet https://github.com/BioinfoUNIBA/REDItools2.git "$DEST"
        fi
        git -C "$DEST" checkout --quiet {params.commit}
        test -f {output.script}
        """


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
        # RepeatMasker .out columns: $10 = repeat name (e.g. AluSx, AluJb),
        # $11 = repeat class/family (e.g. SINE/Alu). Match column 11 to the
        # SINE/Alu class, which captures every Alu subfamily.
        (gunzip -c {input.rmsk} \
            | awk -v has_chr="$HAS_CHR" 'NR>3 && $11 == "SINE/Alu" {{
                chrom=$5;
                if (has_chr=="0") sub(/^chr/, "", chrom);
                else if (chrom !~ /^chr/) chrom="chr"chrom;
                print chrom"\t"($6-1)"\t"$7
              }}' \
            | sort -k1,1 -k2,2n \
            > {output}) 2> {log}
        """
