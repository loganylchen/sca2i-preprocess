"""Aggregate every per-observation annotated parquet into a single sparse
AnnData consumable by ``sca2i.tl.discover()``.

Layout:
- ``adata.X``                = coverage ``n``           (CSR int32, sparse)
- ``adata.layers['edits']``  = edit count ``k``         (CSR int32, sparse)
- ``adata.var``              = unique sites (chrom, pos, strand, ref, alt,
                                 ctx5, ctx3, adar_motif_ok, in_alu)
- ``adata.obs``              = per-observation metadata (sample, donor,
                                 chemistry, celltype, cell_barcode if 10x)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

snakemake = snakemake  # type: ignore[name-defined]

parquets = [Path(p) for p in snakemake.input.parquets]
cells_tsv = Path(snakemake.input.cells_tsv)
samples_tsv = Path(snakemake.input.samples_tsv)
out_path = Path(snakemake.output.h5ad)
out_path.parent.mkdir(parents=True, exist_ok=True)
drop_empty = bool(snakemake.params.drop_empty)


# ---------------------------------------------------------------------------
# Parse obs_id from filename
# ---------------------------------------------------------------------------
# Path layouts written by filters.smk:
#   results/filtered/{sample}/10x__{group}__{chrom}.annot.parquet
#   results/filtered/{sample}/ss__{cell}.annot.parquet
RE_10X = re.compile(r"^10x__(?P<group>[^_]+(?:_[^_]+)*?)__(?P<chrom>[^.]+)\.annot\.parquet$")
RE_SS  = re.compile(r"^ss__(?P<cell>.+)\.annot\.parquet$")


def _parse(path: Path) -> tuple[str, str, str | None]:
    """Return (sample, obs_id, chrom_or_none)."""
    sample = path.parent.name
    name = path.name
    m = RE_10X.match(name)
    if m:
        group = m.group("group")
        chrom = m.group("chrom")
        return sample, f"{sample}__{group}", chrom
    m = RE_SS.match(name)
    if m:
        cell = m.group("cell")
        return sample, f"{sample}__{cell}", None
    raise ValueError(f"Cannot parse obs from filename: {path}")


# ---------------------------------------------------------------------------
# First pass: collect per-obs row groups
# ---------------------------------------------------------------------------
per_obs: dict[str, list[pd.DataFrame]] = {}
for p in parquets:
    sample, obs_id, _chrom = _parse(p)
    df = pd.read_parquet(p)
    if df.empty:
        continue
    df["sample_id"] = sample
    df["obs_id"] = obs_id
    per_obs.setdefault(obs_id, []).append(df)

if not per_obs:
    raise RuntimeError("No non-empty input parquets — nothing to aggregate.")

# Concatenate per-obs and union sites
obs_ids = sorted(per_obs.keys())
site_keys: dict[tuple[str, int, int, str, str], int] = {}
site_meta: list[tuple[str, int, int, str, str, str, str, bool, bool]] = []


def _site_key(row) -> tuple[str, int, int, str, str]:
    return (str(row["chrom"]), int(row["pos"]), int(row["strand"]), str(row["ref"]), str(row["alt"]))


for obs_id in obs_ids:
    big = pd.concat(per_obs[obs_id], ignore_index=True)
    per_obs[obs_id] = [big]  # collapsed
    for _, r in big.iterrows():
        k = _site_key(r)
        if k not in site_keys:
            site_keys[k] = len(site_keys)
            site_meta.append(
                (
                    *k,
                    str(r.get("ctx5", "") or ""),
                    str(r.get("ctx3", "") or ""),
                    bool(r.get("adar_motif_ok", False)),
                    bool(r.get("in_alu", False)),
                )
            )

n_obs = len(obs_ids)
n_var = len(site_keys)
sys.stderr.write(f"Aggregating {n_obs} observations x {n_var} unique sites\n")

# ---------------------------------------------------------------------------
# Build sparse matrices (COO -> CSR)
# ---------------------------------------------------------------------------
rows: list[int] = []
cols: list[int] = []
n_data: list[int] = []
k_data: list[int] = []

for i, obs_id in enumerate(obs_ids):
    big = per_obs[obs_id][0]
    for _, r in big.iterrows():
        j = site_keys[_site_key(r)]
        rows.append(i)
        cols.append(j)
        n_data.append(int(r["n"]))
        k_data.append(int(r["k"]))

X_n = sp.coo_matrix(
    (np.asarray(n_data, dtype=np.int32), (np.asarray(rows), np.asarray(cols))),
    shape=(n_obs, n_var),
).tocsr()
X_k = sp.coo_matrix(
    (np.asarray(k_data, dtype=np.int32), (np.asarray(rows), np.asarray(cols))),
    shape=(n_obs, n_var),
).tocsr()

# ---------------------------------------------------------------------------
# Build obs / var DataFrames
# ---------------------------------------------------------------------------
cells = pd.read_csv(cells_tsv, sep="\t", dtype=str, comment="#").set_index("obs_id")
samples = pd.read_csv(samples_tsv, sep="\t", dtype=str, comment="#").set_index("sample_id")

obs_df = pd.DataFrame(index=pd.Index(obs_ids, name="obs_id"))
for col in ["sample_id", "chemistry", "donor", "celltype", "cell_barcode"]:
    if col in cells.columns:
        obs_df[col] = [cells.at[oid, col] if oid in cells.index else "" for oid in obs_ids]
    else:
        obs_df[col] = ""

var_df = pd.DataFrame(
    site_meta,
    columns=[
        "chrom", "pos", "strand", "ref", "alt",
        "ctx5", "ctx3", "adar_motif_ok", "in_alu",
    ],
)
var_df.index = pd.Index(
    [f"{c}:{p}:{s}:{r}>{a}" for c, p, s, r, a in zip(
        var_df["chrom"], var_df["pos"], var_df["strand"], var_df["ref"], var_df["alt"],
    )],
    name="site_id",
)

adata = ad.AnnData(X=X_n, obs=obs_df, var=var_df)
adata.layers["edits"] = X_k

# ---------------------------------------------------------------------------
# Optional empty-site drop
# ---------------------------------------------------------------------------
if drop_empty:
    nz_per_site = (X_n != 0).sum(axis=0)
    nz_per_site = np.asarray(nz_per_site).ravel()
    keep = nz_per_site > 0
    if not keep.all():
        sys.stderr.write(
            f"Dropping {int((~keep).sum())} all-zero sites; keeping {int(keep.sum())}\n"
        )
        adata = adata[:, keep].copy()

# Tag provenance
adata.uns["sca2i_preprocess"] = {
    "version": "0.1.0",
    "source": "sca2i-preprocess Snakemake workflow",
    "n_input_parquets": len(parquets),
}

adata.write_h5ad(out_path, compression="gzip")
sys.stderr.write(
    f"Wrote AnnData {adata.shape} to {out_path} "
    f"(nnz_X={adata.X.nnz}, nnz_edits={adata.layers['edits'].nnz})\n"
)
