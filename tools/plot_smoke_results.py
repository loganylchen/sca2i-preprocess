"""Render smoke-test diagnostic plots for docs/SMOKE_PBMC1K.md.

Inputs:
    results/anndata/sca2i_input.h5ad

Outputs (PNG, written to docs/img/smoke_pbmc1k/):
    01_sites_per_cell.png         — bar chart, cells sorted by detection rate
    02_edit_fraction_hist.png     — histogram of per-(cell, site) AF = edits / coverage
    03_motif_context_heatmap.png  — 4x4 ctx5 x ctx3 site count, with ADAR-pref overlay
    04_chrom_breakdown.png        — sites per chromosome stacked by motif_ok
    05_sparsity_pattern.png       — sparse spy plot, cells x sites
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp


H5AD = Path("results/anndata/sca2i_input.h5ad")
OUTDIR = Path("docs/img/smoke_pbmc1k")
BASES = ["A", "C", "G", "T"]


def main() -> int:
    if not H5AD.exists():
        sys.stderr.write(f"missing {H5AD}\n")
        return 1
    OUTDIR.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(H5AD)
    n_cells, n_sites = a.shape
    sys.stderr.write(f"loaded {n_cells} cells x {n_sites} sites\n")

    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    E = a.layers["edits"].tocsr() if sp.issparse(a.layers["edits"]) else sp.csr_matrix(a.layers["edits"])

    # ---- (1) sites per cell ----
    sites_per_cell = np.asarray((X > 0).sum(axis=1)).ravel()
    order = np.argsort(-sites_per_cell)
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.bar(np.arange(n_cells), sites_per_cell[order], color="#4C78A8")
    ax.set_xlabel("cell rank")
    ax.set_ylabel("# sites covered")
    ax.set_title(f"Per-cell detection rate (n={n_cells} cells, {n_sites} candidate sites)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "01_sites_per_cell.png", dpi=150)
    plt.close(fig)

    # ---- (2) edit fraction histogram ----
    cov = X.toarray()
    edits = E.toarray()
    mask = cov > 0
    af = np.divide(edits, cov, where=mask, out=np.zeros_like(cov, dtype=float))
    af_vals = af[mask]
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.hist(af_vals, bins=40, color="#F58518", edgecolor="white")
    ax.set_xlabel("edit fraction = edits / coverage")
    ax.set_ylabel("# (cell, site) observations")
    ax.set_title(f"AF distribution across {mask.sum()} non-empty entries")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "02_edit_fraction_hist.png", dpi=150)
    plt.close(fig)

    # ---- (3) motif context heatmap (ctx5 x ctx3) ----
    var = a.var.copy()
    ctx_table = (
        var.groupby(["ctx5", "ctx3"]).size().unstack(fill_value=0).reindex(index=BASES, columns=BASES, fill_value=0)
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    im0 = axes[0].imshow(ctx_table.values, cmap="Blues")
    axes[0].set_xticks(range(4), BASES)
    axes[0].set_yticks(range(4), BASES)
    axes[0].set_xlabel("3' context")
    axes[0].set_ylabel("5' context")
    axes[0].set_title("Site counts by ctx5 x ctx3")
    for i in range(4):
        for j in range(4):
            v = int(ctx_table.values[i, j])
            axes[0].text(j, i, str(v), ha="center", va="center",
                         color="white" if v > ctx_table.values.max() * 0.5 else "black", fontsize=9)
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    # ADAR1-preference overlay: 5' depleted G; 3' enriched G.
    pref = np.zeros((4, 4))
    for i, b5 in enumerate(BASES):
        for j, b3 in enumerate(BASES):
            pref[i, j] = (b5 != "G") and (b3 == "G")
    axes[1].imshow(pref, cmap="Greens", vmin=0, vmax=1)
    axes[1].set_xticks(range(4), BASES)
    axes[1].set_yticks(range(4), BASES)
    axes[1].set_xlabel("3' context")
    axes[1].set_ylabel("5' context")
    n_ok = int(var["adar_motif_ok"].sum())
    axes[1].set_title(f"ADAR1 preference (5'!=G & 3'=G): {n_ok}/{len(var)} sites")
    for i in range(4):
        for j in range(4):
            axes[1].text(j, i, "OK" if pref[i, j] else "", ha="center", va="center",
                         color="black", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUTDIR / "03_motif_context_heatmap.png", dpi=150)
    plt.close(fig)

    # ---- (4) chrom breakdown ----
    chrom_split = (
        var.groupby(["chrom", "adar_motif_ok"]).size().unstack(fill_value=0)
        .rename(columns={False: "motif_fail", True: "motif_ok"})
    )
    chrom_split = chrom_split.reindex(sorted(chrom_split.index, key=lambda c: (c != "X", c)))
    fig, ax = plt.subplots(figsize=(6, 3.5))
    bottom = np.zeros(len(chrom_split))
    for col, color in [("motif_fail", "#BAB0AC"), ("motif_ok", "#54A24B")]:
        if col in chrom_split.columns:
            vals = chrom_split[col].values
            ax.bar(chrom_split.index, vals, bottom=bottom, label=col, color=color)
            bottom = bottom + vals
    ax.set_ylabel("# candidate sites")
    ax.set_title("Sites per chromosome (smoke restricted to chr22, X)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "04_chrom_breakdown.png", dpi=150)
    plt.close(fig)

    # ---- (5) sparsity pattern ----
    fig, ax = plt.subplots(figsize=(9, 3.5))
    coo = X.tocoo()
    ax.scatter(coo.col, coo.row, s=2, c="#4C78A8", marker=",", linewidths=0, alpha=0.7)
    ax.set_xlim(-0.5, n_sites - 0.5)
    ax.set_ylim(n_cells - 0.5, -0.5)
    ax.set_xlabel("site index (chr22 -> chrX)")
    ax.set_ylabel("cell")
    ax.set_title(f"Coverage sparsity pattern  (nnz={X.nnz}/{n_cells * n_sites}, density={X.nnz / (n_cells * n_sites):.2%})")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "05_sparsity_pattern.png", dpi=150)
    plt.close(fig)

    # ---- summary stats ----
    print(pd.DataFrame({
        "metric": [
            "cells", "sites", "median sites/cell", "max sites/cell",
            "nonzero entries", "mean AF", "median AF", "motif_ok sites",
            "chr22 sites", "chrX sites",
        ],
        "value": [
            n_cells, n_sites, int(np.median(sites_per_cell)), int(sites_per_cell.max()),
            int(mask.sum()), float(af_vals.mean()), float(np.median(af_vals)),
            int(var["adar_motif_ok"].sum()),
            int((var["chrom"] == "22").sum()), int((var["chrom"] == "X").sum()),
        ],
    }).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
