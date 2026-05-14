"""Convert a REDItools2 TSV(.gz) to parquet, restricted to A>G (and T>C on
reverse strand) candidate edits.

REDItools2 columns (canonical, v1.x):
    Region  Position  Reference  Strand  Coverage-q30  MeanQ
    BaseCount[A,C,G,T]  AllSubs  Frequency  gCoverage-q30  gMeanQ
    gBaseCount[A,C,G,T]  gAllSubs  gFrequency

We keep:
    chrom, pos (1-based), strand, ref, alt, k (alt count), n (coverage), AF
"""
from __future__ import annotations

from pathlib import Path
import gzip
import sys

import pandas as pd

snakemake = snakemake  # type: ignore[name-defined]

in_path = Path(snakemake.input.tsv)
out_path = Path(snakemake.output.parquet)
out_path.parent.mkdir(parents=True, exist_ok=True)


def _parse_basecount(s: str) -> tuple[int, int, int, int]:
    """REDItools writes BaseCount as e.g. '[12, 0, 3, 0]'."""
    s = s.strip().lstrip("[").rstrip("]")
    parts = [int(x.strip()) for x in s.split(",")]
    if len(parts) != 4:
        raise ValueError(f"Bad BaseCount field: {s!r}")
    a, c, g, t = parts
    return a, c, g, t


def _iter_rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_region   = header.index("Region")
            i_pos      = header.index("Position")
            i_ref      = header.index("Reference")
            i_strand   = header.index("Strand")
            i_cov      = header.index("Coverage-q30")
            i_basect   = header.index("BaseCount[A,C,G,T]")
        except ValueError as exc:
            raise RuntimeError(f"Unexpected REDItools header in {path}: {header}") from exc
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ref = f[i_ref]
            try:
                strand = int(f[i_strand])
            except ValueError:
                strand = 0
            n = int(f[i_cov])
            a, c, g, t = _parse_basecount(f[i_basect])
            # A>G on + strand, T>C on - strand (=A>G on the transcribed strand)
            if ref == "A":
                k = g
                alt = "G"
            elif ref == "T":
                k = c
                alt = "C"
            else:
                continue
            if k <= 0 or n <= 0:
                continue
            yield {
                "chrom":  f[i_region],
                "pos":    int(f[i_pos]),
                "strand": strand,
                "ref":    ref,
                "alt":    alt,
                "k":      k,
                "n":      n,
                "AF":     k / n,
            }


rows = list(_iter_rows(in_path))
df = pd.DataFrame.from_records(
    rows,
    columns=["chrom", "pos", "strand", "ref", "alt", "k", "n", "AF"],
)
# Down-cast for storage
df = df.astype(
    {"pos": "int64", "strand": "int8", "k": "int32", "n": "int32", "AF": "float32"}
)
df.to_parquet(out_path, compression="zstd", index=False)
sys.stderr.write(f"Wrote {len(df)} A>I candidate rows to {out_path}\n")
