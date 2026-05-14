"""Annotate ADAR sequence-context motif. Adds two columns to the parquet:

- ``ctx5``: ``window``-bp upstream base(s) on the transcribed strand.
- ``ctx3``: ``window``-bp downstream base(s) on the transcribed strand.
- ``adar_motif_ok``: bool — True if 5' base is U/A (depleted G) AND 3' base is G
  (enriched). Standard ADAR1 sequence preference (Eggington et al. 2011).

This is annotation only — no rows are dropped. The downstream scA2I discoverer
decides whether to use the column as a hard filter.
"""
from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
from pyfaidx import Fasta

snakemake = snakemake  # type: ignore[name-defined]

in_parquet = Path(snakemake.input.parquet)
fa_path = Path(snakemake.input.fa)
out_parquet = Path(snakemake.output.parquet)
out_parquet.parent.mkdir(parents=True, exist_ok=True)

window = int(snakemake.params.window)
do_run = bool(snakemake.params.do_run)

df = pd.read_parquet(in_parquet)

if not do_run or df.empty:
    df["ctx5"] = pd.Series([""] * len(df), dtype="string")
    df["ctx3"] = pd.Series([""] * len(df), dtype="string")
    df["adar_motif_ok"] = pd.Series([False] * len(df), dtype="boolean")
    df.to_parquet(out_parquet, compression="zstd", index=False)
    sys.stderr.write(f"Motif annotation skipped (do_run={do_run}, rows={len(df)})\n")
    sys.exit(0)

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


fa = Fasta(str(fa_path), as_raw=True, sequence_always_upper=True)

ctx5: list[str] = []
ctx3: list[str] = []
ok: list[bool] = []

for chrom, pos, strand in zip(df["chrom"], df["pos"], df["strand"]):
    pos = int(pos)  # 1-based
    if chrom not in fa:
        ctx5.append("")
        ctx3.append("")
        ok.append(False)
        continue
    seq = fa[chrom]
    seq_len = len(seq)
    lo = max(0, pos - 1 - window)
    hi = min(seq_len, pos + window)
    region = str(seq[lo:hi])  # 0-based slice
    if len(region) < 1 + 2 * window:
        ctx5.append("")
        ctx3.append("")
        ok.append(False)
        continue
    # Revcomp based on ref base (DNA-strand of the edited site), not REDItools
    # strand code, so chemistry-aware strand handling upstream stays decoupled
    # from motif context extraction.
    seq_ref = region[window:window + 1]
    if seq_ref == "T":
        region = _revcomp(region)
    five = region[:window]
    three = region[-window:]
    ctx5.append(five)
    ctx3.append(three)
    # ADAR1 preference: 5' depleted G (so T/A at -1 is "ok"); 3' enriched G.
    # Note: DNA-only reference — no "U" base appears post-revcomp.
    five_ok = five[-1] in {"T", "A"}
    three_ok = three[0] == "G"
    ok.append(bool(five_ok and three_ok))

df["ctx5"] = pd.Series(ctx5, dtype="string")
df["ctx3"] = pd.Series(ctx3, dtype="string")
df["adar_motif_ok"] = pd.Series(ok, dtype="boolean")
df.to_parquet(out_parquet, compression="zstd", index=False)
sys.stderr.write(
    f"Motif annotated {len(df)} rows; {int(df['adar_motif_ok'].fillna(False).sum())} pass ADAR motif\n"
)
