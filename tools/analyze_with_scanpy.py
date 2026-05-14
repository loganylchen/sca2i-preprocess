"""Standard scanpy single-cell pipeline on PBMC 1k v3, overlaid with A-to-I
editing metrics from sca2i-preprocess.

Inputs:
    results/data/pbmc1k/filtered_feature_bc_matrix/        (10x v3 mtx)
    results/anndata/sca2i_input.annotated.h5ad             (or sca2i_input.h5ad)

Outputs (docs/img/smoke_pbmc1k/):
    06_umap_celltype.png             leiden + PBMC-marker cell type assignment
    07_umap_smoke_overlay.png        49 smoke cells highlighted on the full UMAP
    08_umap_edits.png                edit metrics colored on smoke cells
    09_gene_annotation_breakdown.png GTF region + biotype site breakdown
    10_editing_by_celltype.png       per-cluster editing rate violin/box

Side outputs:
    results/anndata/pbmc1k_scanpy.h5ad  (full 10x adata with editing metrics
                                         merged on shared barcodes)
"""

import sys
import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")

MTX_DIR = Path("results/data/pbmc1k/filtered_feature_bc_matrix")
SCA2I_ANNOT = Path("results/anndata/sca2i_input.annotated.h5ad")
SCA2I_PLAIN = Path("results/anndata/sca2i_input.h5ad")
OUTDIR = Path("docs/img/smoke_pbmc1k")
SCANPY_OUT = Path("results/anndata/pbmc1k_scanpy.h5ad")

# Compact PBMC marker set (Seurat tutorial)
PBMC_MARKERS = {
    "T (CD4)":  ["IL7R", "CCR7"],
    "T (CD8)":  ["CD8A", "CD8B"],
    "T_naive":  ["LEF1", "TCF7"],
    "NK":       ["GNLY", "NKG7", "KLRD1"],
    "B":        ["MS4A1", "CD79A"],
    "Mono_CD14":["CD14", "LYZ"],
    "Mono_FCGR3A":["FCGR3A", "MS4A7"],
    "DC":       ["FCER1A", "CST3"],
    "Platelet": ["PPBP"],
}

sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)
sc.settings.verbosity = 1


# ---------------------------------------------------------------------------
def load_expression() -> ad.AnnData:
    if not MTX_DIR.exists():
        sys.exit(f"missing 10x mtx dir: {MTX_DIR}")
    sys.stderr.write(f"reading 10x mtx from {MTX_DIR}\n")
    adata = sc.read_10x_mtx(MTX_DIR, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    sys.stderr.write(f"  -> {adata.shape}\n")
    return adata


def qc_and_cluster(adata: ad.AnnData) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
    sys.stderr.write(f"after QC: {adata.shape}, {adata.obs['leiden'].nunique()} clusters\n")
    return adata


def annotate_celltypes(adata: ad.AnnData) -> ad.AnnData:
    """Score each cluster against PBMC marker panels and assign best label."""
    raw = adata.raw.to_adata() if adata.raw is not None else adata
    for ct, genes in PBMC_MARKERS.items():
        present = [g for g in genes if g in raw.var_names]
        if not present:
            adata.obs[f"score_{ct}"] = 0.0
            continue
        sc.tl.score_genes(raw, gene_list=present, score_name=f"score_{ct}", use_raw=False)
        adata.obs[f"score_{ct}"] = raw.obs[f"score_{ct}"].values

    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    mean_by_cluster = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    best = mean_by_cluster.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["cell_type"] = adata.obs["leiden"].map(best.to_dict()).astype("category")
    return adata


def merge_editing_metrics(adata: ad.AnnData, sca: ad.AnnData) -> ad.AnnData:
    """Match by trailing barcode portion (sca2i obs is 'sample__BARCODE-1')."""
    sca_bc = []
    for oid in sca.obs_names:
        bc = oid.split("__", 1)[1] if "__" in oid else oid
        sca_bc.append(bc)
    sca_bc = np.asarray(sca_bc)

    X = sca.X.tocsr() if sp.issparse(sca.X) else sp.csr_matrix(sca.X)
    E = sca.layers["edits"].tocsr() if sp.issparse(sca.layers["edits"]) else sp.csr_matrix(sca.layers["edits"])
    total_cov = np.asarray(X.sum(axis=1)).ravel()
    total_edit = np.asarray(E.sum(axis=1)).ravel()
    n_sites = np.asarray((X > 0).sum(axis=1)).ravel()
    n_edited = np.asarray((E > 0).sum(axis=1)).ravel()
    global_af = np.where(total_cov > 0, total_edit / np.maximum(total_cov, 1), np.nan)

    metric_df = pd.DataFrame({
        "barcode": sca_bc,
        "edit_total_cov": total_cov,
        "edit_total_k":   total_edit,
        "edit_n_sites":   n_sites,
        "edit_n_edited":  n_edited,
        "edit_global_af": global_af,
    }).set_index("barcode")

    bc_idx = pd.Index(adata.obs_names)
    overlap = metric_df.index.intersection(bc_idx)
    sys.stderr.write(f"smoke cells matched in expression UMAP: {len(overlap)}/{len(sca_bc)}\n")

    for col in metric_df.columns:
        adata.obs[col] = np.nan
        adata.obs.loc[overlap, col] = metric_df.loc[overlap, col].values
    adata.obs["in_smoke_set"] = adata.obs_names.isin(metric_df.index)
    return adata


# ---------------------------------------------------------------------------
def plot_umap_celltype(adata: ad.AnnData, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    sc.pl.umap(adata, color="leiden", ax=axes[0], show=False, frameon=False, title="leiden")
    sc.pl.umap(adata, color="cell_type", ax=axes[1], show=False, frameon=False, title="PBMC cell type")
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_smoke_overlay(adata: ad.AnnData, out: Path) -> None:
    coords = adata.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(coords[:, 0], coords[:, 1], s=6, c="#cccccc", alpha=0.6, label=f"PBMC 1k (n={adata.n_obs})")
    mask = adata.obs["in_smoke_set"].values
    ax.scatter(coords[mask, 0], coords[mask, 1], s=28, c="#E45756", edgecolor="black",
               linewidth=0.4, label=f"smoke cells (n={mask.sum()})")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title("PBMC 1k v3 UMAP — smoke cells overlaid")
    ax.legend(frameon=False, loc="lower left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_edit_metrics(adata: ad.AnnData, out: Path) -> None:
    metrics = [
        ("edit_total_cov", "total coverage (sum over sites)"),
        ("edit_total_k",   "total edits (sum over sites)"),
        ("edit_n_sites",   "# sites covered"),
        ("edit_global_af", "global AF (edits / coverage)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    coords = adata.obsm["X_umap"]
    smoke_mask = adata.obs["in_smoke_set"].values
    for ax, (col, label) in zip(axes.ravel(), metrics):
        ax.scatter(coords[~smoke_mask, 0], coords[~smoke_mask, 1], s=5, c="#e5e5e5", alpha=0.5)
        vals = adata.obs.loc[smoke_mask, col].astype(float).values
        sc_ = ax.scatter(coords[smoke_mask, 0], coords[smoke_mask, 1], c=vals,
                         cmap="viridis", s=36, edgecolor="black", linewidth=0.3)
        fig.colorbar(sc_, ax=ax, fraction=0.04, pad=0.04)
        ax.set_title(label)
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("A-to-I editing metrics on the 49 smoke cells", y=1.005)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_gene_annotation(sca: ad.AnnData, out: Path) -> None:
    var = sca.var.copy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    region_counts = var["region"].value_counts() if "region" in var.columns else pd.Series(dtype=int)
    if region_counts.empty:
        axes[0].text(0.5, 0.5, "no region annotation\n(run annotate_sites_gtf.py)",
                     ha="center", va="center", transform=axes[0].transAxes)
    else:
        order = ["exon", "intron", "intergenic"]
        order = [r for r in order if r in region_counts.index] + \
                [r for r in region_counts.index if r not in order]
        region_counts = region_counts.reindex(order)
        axes[0].bar(region_counts.index, region_counts.values,
                    color=["#54A24B", "#F58518", "#BAB0AC"][:len(region_counts)])
        for i, v in enumerate(region_counts.values):
            axes[0].text(i, v, str(int(v)), ha="center", va="bottom", fontsize=10)
        axes[0].set_ylabel("# candidate sites")
        axes[0].set_title("Sites by genomic region (Ensembl r110)")
        axes[0].spines[["top", "right"]].set_visible(False)

    if "gene_biotype" in var.columns:
        bt = var.loc[var["gene_biotype"] != "", "gene_biotype"].value_counts().head(8)
        axes[1].barh(bt.index[::-1], bt.values[::-1], color="#4C78A8")
        for i, v in enumerate(bt.values[::-1]):
            axes[1].text(v, i, f" {int(v)}", va="center", fontsize=9)
        axes[1].set_xlabel("# sites")
        axes[1].set_title("Top gene biotypes")
        axes[1].spines[["top", "right"]].set_visible(False)
    else:
        axes[1].text(0.5, 0.5, "no biotype annotation", ha="center", va="center",
                     transform=axes[1].transAxes)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_edits_by_celltype(adata: ad.AnnData, out: Path) -> None:
    sub = adata.obs.loc[adata.obs["in_smoke_set"]].copy()
    if sub.empty:
        return
    sub["cell_type"] = sub["cell_type"].astype(str)
    metric = "edit_global_af"
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    order = sub["cell_type"].value_counts().index.tolist()
    sns.boxplot(data=sub, x="cell_type", y=metric, order=order, color="#9ECAE1", ax=axes[0])
    sns.stripplot(data=sub, x="cell_type", y=metric, order=order, color="black", size=3, ax=axes[0])
    axes[0].set_title("Global AF per cell type")
    axes[0].set_xlabel(""); axes[0].set_ylabel("edits / coverage")
    axes[0].tick_params(axis="x", rotation=30)
    axes[0].spines[["top", "right"]].set_visible(False)

    counts = sub["cell_type"].value_counts().reindex(order, fill_value=0)
    axes[1].bar(counts.index, counts.values, color="#B279A2")
    for i, v in enumerate(counts.values):
        axes[1].text(i, v, str(int(v)), ha="center", va="bottom", fontsize=9)
    axes[1].set_title("# smoke cells per cell type")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("# cells")
    axes[1].tick_params(axis="x", rotation=30)
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sca_path = SCA2I_ANNOT if SCA2I_ANNOT.exists() else SCA2I_PLAIN
    if not sca_path.exists():
        sys.exit(f"missing sca2i AnnData: {sca_path}")
    sys.stderr.write(f"reading {sca_path}\n")
    sca = ad.read_h5ad(sca_path)

    adata = load_expression()
    adata = qc_and_cluster(adata)
    adata = annotate_celltypes(adata)
    adata = merge_editing_metrics(adata, sca)

    SCANPY_OUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(SCANPY_OUT, compression="gzip")
    sys.stderr.write(f"wrote {SCANPY_OUT}\n")

    plot_umap_celltype(adata,  OUTDIR / "06_umap_celltype.png")
    plot_smoke_overlay(adata,  OUTDIR / "07_umap_smoke_overlay.png")
    plot_edit_metrics(adata,   OUTDIR / "08_umap_edits.png")
    plot_gene_annotation(sca,  OUTDIR / "09_gene_annotation_breakdown.png")
    plot_edits_by_celltype(adata, OUTDIR / "10_editing_by_celltype.png")

    # Summary printout
    print(adata.obs["cell_type"].value_counts().to_string())
    if "region" in sca.var.columns:
        print("--- site regions ---")
        print(sca.var["region"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
