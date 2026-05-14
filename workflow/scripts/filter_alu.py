"""Annotate Alu repeat overlap (Critic NN#5).

Adds a single boolean column ``in_alu``. Annotation only — no rows are dropped.
"""
from __future__ import annotations

from pathlib import Path
import sys
import tempfile

import pandas as pd

snakemake = snakemake  # type: ignore[name-defined]

in_parquet = Path(snakemake.input.parquet)
bed_path = Path(snakemake.input.bed)
out_parquet = Path(snakemake.output.parquet)
out_parquet.parent.mkdir(parents=True, exist_ok=True)

do_run = bool(snakemake.params.do_run)

df = pd.read_parquet(in_parquet)

if not do_run or df.empty:
    df["in_alu"] = pd.Series([False] * len(df), dtype="boolean")
    df.to_parquet(out_parquet, compression="zstd", index=False)
    sys.stderr.write(f"Alu annotation skipped (do_run={do_run}, rows={len(df)})\n")
    sys.exit(0)

import pybedtools  # noqa: E402

with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
    sites_bed = Path(tmp.name)
    for _, row in df.iterrows():
        # 0-based half-open BED
        start = int(row["pos"]) - 1
        end = int(row["pos"])
        tmp.write(f"{row['chrom']}\t{start}\t{end}\n")

try:
    sites = pybedtools.BedTool(str(sites_bed)).sort()
    alu = pybedtools.BedTool(str(bed_path))
    overlaps = sites.intersect(alu, c=True, sorted=False)
    in_alu_lookup: dict[tuple[str, int], bool] = {}
    for feat in overlaps:
        chrom = feat.chrom
        end = int(feat.end)  # 1-based pos
        cnt = int(feat.fields[-1])
        in_alu_lookup[(chrom, end)] = cnt > 0
finally:
    sites_bed.unlink(missing_ok=True)

in_alu = [
    in_alu_lookup.get((c, int(p)), False)
    for c, p in zip(df["chrom"], df["pos"])
]
df["in_alu"] = pd.Series(in_alu, dtype="boolean")
df.to_parquet(out_parquet, compression="zstd", index=False)
sys.stderr.write(
    f"Alu annotated {len(df)} rows; {int(df['in_alu'].fillna(False).sum())} fall in Alu\n"
)
