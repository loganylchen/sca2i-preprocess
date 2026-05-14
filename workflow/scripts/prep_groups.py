"""Build a sinto-compatible barcode -> group TSV from the cells.tsv sheet.

For 10x samples, cell_barcode is expected to be either:
- a single comma-separated string of barcodes for that group, OR
- empty, in which case the script looks for a sibling file
  ``<cell_bam_dir>/<obs_id>.barcodes.txt`` (one barcode per line).

Output: 2-column TSV (barcode<TAB>group_name), as required by ``sinto filterbarcodes -c``.
"""
from __future__ import annotations

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
    group = row["cell_type"]
    bcs_raw = (row.get("cell_barcode") or "").strip()
    if bcs_raw:
        for bc in bcs_raw.split(","):
            bc = bc.strip()
            if bc:
                rows.append((bc, group))
    else:
        # Look for sidecar file
        sidecar = cells_tsv.parent / f"{row['obs_id']}.barcodes.txt"
        if not sidecar.exists():
            sys.stderr.write(
                f"obs_id {row['obs_id']}: no barcodes inline and no sidecar at {sidecar}\n"
            )
            sys.exit(1)
        for bc in sidecar.read_text().splitlines():
            bc = bc.strip()
            if bc:
                rows.append((bc, group))

out_df = pd.DataFrame(rows, columns=["barcode", "group"]).drop_duplicates()
out_df.to_csv(out_path, sep="\t", index=False, header=False)
sys.stderr.write(
    f"Wrote {len(out_df)} barcode->group rows ({out_df['group'].nunique()} groups) to {out_path}\n"
)
