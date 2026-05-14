"""One-shot helper: build samples.tsv + cells.tsv for the PBMC 1k v3 smoke test.

Reads filtered barcodes from the cellranger demo bundle (barcodes.tsv.gz) and
emits a per-cell sheet (one row per subsampled barcode).

Usage (run from repo root):
    python tools/make_cells_tsv.py \
        --barcodes results/data/pbmc1k/filtered_feature_bc_matrix/barcodes.tsv.gz \
        --bam      results/data/pbmc1k/pbmc_1k_v3_bam.bam \
        --sample   pbmc1k \
        --n-cells  50 \
        --seed     0

This rewrites config/samples.tsv and config/cells.tsv. The previous example_*
stub rows are intentionally replaced so snakemake only runs the smoke-test
sample. Original stubs live in git history.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from random import Random


def _read_barcodes(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--barcodes", required=True, type=Path)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--sample", default="pbmc1k")
    ap.add_argument("--chemistry", default="10x")
    ap.add_argument("--donor", default="pbmc_donor")
    ap.add_argument("--cell-type", default="PBMC")
    ap.add_argument("--n-cells", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--samples-tsv", default="config/samples.tsv", type=Path)
    ap.add_argument("--cells-tsv", default="config/cells.tsv", type=Path)
    args = ap.parse_args()

    if not args.barcodes.exists():
        sys.stderr.write(f"barcodes file not found: {args.barcodes}\n")
        return 1
    if not args.bam.exists():
        sys.stderr.write(f"BAM file not found: {args.bam}\n")
        return 1

    bam_abs = args.bam.resolve()
    barcodes_abs = args.barcodes.resolve()

    bcs = _read_barcodes(args.barcodes)
    if not bcs:
        sys.stderr.write(f"no barcodes parsed from {args.barcodes}\n")
        return 1
    rng = Random(args.seed)
    rng.shuffle(bcs)
    subset = bcs[: args.n_cells]
    sys.stderr.write(
        f"Picked {len(subset)} of {len(bcs)} barcodes (seed={args.seed})\n"
    )

    # samples.tsv
    samples_lines = [
        "sample_id\tchemistry\tdonor\tbam\tcell_bam_dir\tbarcodes_tsv",
        f"{args.sample}\t{args.chemistry}\t{args.donor}\t{bam_abs}\t\t{barcodes_abs}",
    ]
    args.samples_tsv.write_text("\n".join(samples_lines) + "\n")
    sys.stderr.write(f"Wrote {args.samples_tsv}\n")

    # cells.tsv (one row per subsampled barcode)
    header = "obs_id\tsample_id\tchemistry\tdonor\tcell_type\tcell_barcode\tbam_basename"
    rows = [header]
    for bc in subset:
        # Make obs_id filesystem-safe by replacing '-' nothing-needed (allowed
        # by wildcard_constraints). Keep barcode verbatim.
        obs_id = f"{args.sample}__{bc}"
        rows.append(
            f"{obs_id}\t{args.sample}\t{args.chemistry}\t{args.donor}\t{args.cell_type}\t{bc}\t"
        )
    args.cells_tsv.write_text("\n".join(rows) + "\n")
    sys.stderr.write(f"Wrote {args.cells_tsv} ({len(subset)} cells)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
