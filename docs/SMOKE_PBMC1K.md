# PBMC 1k v3 Per-Cell Smoke Test

End-to-end validation of `sca2i-preprocess` on the 10x Genomics PBMC 1k v3
public demo with per-barcode 10x architecture (one obs row per cell, not
pseudobulk by celltype).

## Inputs

| Source | URL / Path |
|---|---|
| BAM | `https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam` (4.1 GB) |
| BAI | `https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam.bai` |
| Filtered barcodes | `pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz` → `barcodes.tsv.gz` (1222 cells) |
| Genome FASTA | Ensembl `release-110` GRCh38 primary assembly (contigs without `chr` prefix — matches BAM) |
| RepeatMasker | UCSC `hg38.fa.out.gz` |
| REDItools2 | `git@1e9d396f8aba058f073d7d27bf1148fe408adaf8` (pinned) |

## Configuration deviations from production

These are smoke-only relaxations; revert before running production batches.

| Key | Smoke | Production |
|---|---|---|
| `reditools.min_coverage` | 3 | 20 (parent NN#3 lock) |
| `genome.chromosomes` | `["22", "X"]` | full autosomes + X (Y/MT excluded) |
| Cells per sample | 50 (seed=0 subsample via `tools/make_cells_tsv.py`) | full barcode list |
| Cores | 1 | per-cluster (see `profiles/slurm`) |

## How to reproduce

```bash
git checkout validation/pbmc1k-percell

# 1) Stage BAM + barcodes manually into results/data/pbmc1k/
#    (the BAM is too large to host in-repo)
mkdir -p results/data/pbmc1k
# ... curl / scp the 4 files listed above ...
tar -xzf results/data/pbmc1k/pbmc_1k_v3_matrix.tar.gz \
    -C results/data/pbmc1k/

# 2) Generate samples.tsv + cells.tsv for the smoke run
python tools/make_cells_tsv.py \
    --barcodes results/data/pbmc1k/filtered_feature_bc_matrix/barcodes.tsv.gz \
    --bam      results/data/pbmc1k/pbmc_1k_v3_bam.bam \
    --sample   pbmc1k --n-cells 50 --seed 0

# 3) Pre-build conda envs
snakemake --use-conda --conda-create-envs-only --cores 1

# 4) Run. PYTHONNOUSERSITE avoids shadowing by ~/.local pandas.
PYTHONNOUSERSITE=1 snakemake --use-conda --cores 1 -p --keep-going
```

## Acceptance result

`results/anndata/sca2i_input.h5ad`:

| Metric | Value |
|---|---|
| shape (cells × sites) | **49 × 2147** |
| X dtype, nnz | int32, 2353 |
| `layers["edits"]` dtype, nnz | int32, 2353 |
| `obs.columns` | `sample_id, chemistry, donor, cell_type, cell_barcode` |
| `var.columns` | `chrom, pos, strand, ref, alt, ctx5, ctx3, adar_motif_ok, in_alu` |
| Chromosomes seen | `22, X` |
| `adar_motif_ok=True` | 356 / 2147 (16.6%) |
| Wall time | ~32 min (single core) |
| Snakemake jobs | 502 / 502 |

One cell of 50 dropped: it had zero candidate edits passing `min_coverage=3`
on chr22 + chrX.

## Fixes captured on this branch

- **`from __future__ import annotations`** removed from all snakemake
  `script:` files — Snakemake 8/9 injects its preamble before the user
  script body, which makes `from __future__` no longer be the first
  statement and raises `SyntaxError`.
- **`workflow/envs/sinto.yaml`** pins `setuptools <80` — sinto imports
  `pkg_resources` at startup, which setuptools >=80 dropped.
- **REDItools2 install** switched from broken `pip install
  git+...` (the repo has no `setup.py`) to a new
  `clone_reditools2` rule that clones the pinned commit into
  `resources/reditools2/`; `reditools.smk` invokes the resulting
  `reditools.py` via `python`.
- **`split_bam_per_cell`** outputs a sentinel file
  (`results/prep/10x/{sample}/split/.sentinel`) instead of
  `directory(...)`, avoiding Snakemake 9's `ChildIOException` when
  downstream rules write `.bai` siblings into the same directory.
- **Per-cell 10x architecture**: `common.groups_for_sample()` now returns
  unique `cell_barcode` values; `prep_groups.py` writes bc→bc rows so
  sinto emits one BAM per cell. `cells.tsv` carries one obs row per cell.

## Known pre-existing defects (not addressed here)

- `make_alu_bed` filters `$11 ~ /^Alu/` but in RepeatMasker `.out` the
  Alu family name is in column 10 (column 11 is the class
  `SINE/Alu`). Effect: `resources/alu/alu.bed` is empty; `in_alu`
  annotation is all-False. Annotation-only — no impact on counts. One-char
  fix when production-readying.

## Runtime env requirement

Set `PYTHONNOUSERSITE=1` when invoking snakemake on this machine — there
is a partially-broken `~/.local/lib/python3.12/site-packages/pandas` that
will shadow conda env pandas (missing `pytz`) and break every Python rule.
Consider adding to the user's shell profile, or document in `profiles/`.
