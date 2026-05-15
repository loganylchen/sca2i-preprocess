"""Build a sinto-compatible barcode -> group TSV from the cells.tsv sheet.

Per-cell 10x mode: one row per cell, with group_name == cell_barcode so sinto
emits one BAM per cell.

Output: 2-column TSV (barcode<TAB>group_name), as required by ``sinto
filterbarcodes -c``.
"""

from pathlib import Path
import sys

import pandas as pd

snakemake = snakemake  # type: ignore[name-defined]  # injected by Snakemake

cells_tsv = Path(snakemake.input.cells_tsv)
sample = snakemake.wildcards.sample
out_path = Path(snakemake.output.groups_tsv)
out_path.parent.mkdir(parents=True, exist_ok=True)

cells = pd.read_csv(cells_tsv, sep="\t", dtype=str, comment="#")
sub = cells[(cells["sample_id"] == sample) & (cells["chemistry"] == "10x")]
if sub.empty:
    sys.stderr.write(f"No 10x rows in {cells_tsv} for sample {sample}\n")
    sys.exit(1)

rows: list[tuple[str, str]] = []
for _, row in sub.iterrows():
    bc = (row.get("cell_barcode") or "").strip()
    if not bc:
        sys.stderr.write(
            f"obs_id {row['obs_id']}: missing cell_barcode in per-cell mode\n"
        )
        sys.exit(1)
    rows.append((bc, bc))  # group_name = barcode -> one BAM per cell

out_df = pd.DataFrame(rows, columns=["barcode", "group"]).drop_duplicates()
out_df.to_csv(out_path, sep="\t", index=False, header=False)
sys.stderr.write(
    f"Wrote {len(out_df)} per-cell barcode->group rows to {out_path}\n"
)
