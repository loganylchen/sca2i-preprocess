"""Annotate candidate edit sites with Ensembl gene context.

Adds the following columns to ``var`` of ``results/anndata/sca2i_input.h5ad``:

    gene_id        primary overlapping gene (first hit, alphabetical)
    gene_name      Ensembl gene symbol
    gene_biotype   protein_coding / lncRNA / miRNA / ...
    gene_strand    + / -
    region         exon | intron | intergenic
    in_utr         True if site falls in a 5' or 3' UTR

Inputs:
    results/anndata/sca2i_input.h5ad           (read+write — added columns)
    resources/gtf/Homo_sapiens.GRCh38.110.gtf.gz

Outputs:
    results/anndata/sca2i_input.annotated.h5ad
    results/anndata/sites_annotated.tsv
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import pyranges as pr


H5AD_IN = Path("results/anndata/sca2i_input.h5ad")
H5AD_OUT = Path("results/anndata/sca2i_input.annotated.h5ad")
TSV_OUT = Path("results/anndata/sites_annotated.tsv")
GTF = Path("resources/gtf/Homo_sapiens.GRCh38.110.gtf.gz")
CHROMS = {str(c) for c in list(range(1, 23)) + ["X"]}


def _site_pr(var: pd.DataFrame) -> pr.PyRanges:
    df = pd.DataFrame({
        "Chromosome": var["chrom"].astype(str).values,
        "Start": var["pos"].astype(int).values - 1,   # 1-based -> 0-based half-open
        "End":   var["pos"].astype(int).values,
        "site_idx": pd.Series(range(len(var)), dtype="int64").values,
    })
    return pr.PyRanges(df)


def _load_gtf(path: Path, chroms: set[str]) -> pr.PyRanges:
    sys.stderr.write(f"reading {path} ...\n")
    gtf = pr.read_gtf(str(path))
    # PyRanges keeps chromosome names verbatim from the GTF (Ensembl uses bare
    # contig names like "22", "X" — already matches our BAM convention).
    gtf = gtf[gtf.Chromosome.isin(list(chroms))]
    return gtf


def _first_hit(joined: pd.DataFrame, key: str, fields: list[str]) -> pd.DataFrame:
    """Collapse a one-to-many overlap to a single row per site (lex-smallest gene_id)."""
    if joined.empty:
        return pd.DataFrame(columns=[key, *fields])
    # Keep a deterministic pick across runs.
    joined = joined.sort_values([key, "gene_id"])
    return joined.drop_duplicates(subset=[key], keep="first")[[key, *fields]]


def main() -> int:
    if not H5AD_IN.exists():
        sys.stderr.write(f"missing {H5AD_IN}\n"); return 1
    if not GTF.exists():
        sys.stderr.write(f"missing {GTF}\n"); return 1

    a = ad.read_h5ad(H5AD_IN)
    sys.stderr.write(f"loaded {a.shape[0]} cells x {a.shape[1]} sites\n")
    var = a.var.copy()
    sites = _site_pr(var)

    gtf = _load_gtf(GTF, CHROMS)
    genes  = gtf[gtf.Feature == "gene"]
    exons  = gtf[gtf.Feature == "exon"]
    utrs5  = gtf[gtf.Feature == "five_prime_utr"]
    utrs3  = gtf[gtf.Feature == "three_prime_utr"]

    # ---- gene overlap ----
    gene_join = sites.join(genes).df
    gene_cols = ["site_idx", "gene_id", "gene_name", "gene_biotype", "Strand"]
    have = [c for c in gene_cols if c in gene_join.columns]
    gene_join = gene_join[have].rename(columns={"Strand": "gene_strand"})
    gene_hit = _first_hit(
        gene_join, "site_idx", ["gene_id", "gene_name", "gene_biotype", "gene_strand"],
    )

    # ---- exon overlap (boolean) ----
    in_exon = sites.join(exons).df["site_idx"].unique() if not exons.empty else []
    in_utr = (
        set(sites.join(utrs5).df["site_idx"].unique() if not utrs5.empty else [])
        | set(sites.join(utrs3).df["site_idx"].unique() if not utrs3.empty else [])
    )

    # ---- assemble ----
    out = pd.DataFrame({"site_idx": range(len(var))})
    out = out.merge(gene_hit, on="site_idx", how="left")
    out["in_exon"] = out["site_idx"].isin(in_exon)
    out["in_utr"] = out["site_idx"].isin(in_utr)
    out["region"] = "intergenic"
    out.loc[out["gene_id"].notna() & ~out["in_exon"], "region"] = "intron"
    out.loc[out["gene_id"].notna() &  out["in_exon"], "region"] = "exon"

    # Promote to var (cast away any Categorical dtype from PyRanges/GTF parsing
    # before fillna, otherwise "" is not in the category set).
    for col in ("gene_id", "gene_name", "gene_biotype", "gene_strand"):
        var[col] = out[col].astype("object").fillna("").astype(str).values
    var["region"] = out["region"].values
    var["in_utr"] = out["in_utr"].values

    a.var = var
    H5AD_OUT.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(H5AD_OUT, compression="gzip")
    sys.stderr.write(f"wrote {H5AD_OUT}\n")

    # Companion TSV for grep-style inspection
    tsv = var.reset_index().rename(columns={"index": "site_id"})
    tsv.to_csv(TSV_OUT, sep="\t", index=False)
    sys.stderr.write(f"wrote {TSV_OUT} ({len(tsv)} rows)\n")

    # Summary
    sys.stderr.write("---- region breakdown ----\n")
    sys.stderr.write(var["region"].value_counts().to_string() + "\n")
    sys.stderr.write("---- top gene_biotypes ----\n")
    sys.stderr.write(var["gene_biotype"].value_counts().head(10).to_string() + "\n")
    sys.stderr.write("---- top genes ----\n")
    sys.stderr.write(
        var.loc[var["gene_name"] != "", "gene_name"].value_counts().head(15).to_string() + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
