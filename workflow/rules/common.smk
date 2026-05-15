"""Shared helpers — chemistry dispatch, observation-unit fan-out.

Per-cell 10x mode: each row in cells.tsv with chemistry='10x' represents a
single cell barcode (one obs row). ``groups_for_sample`` returns the unique
cell_barcode values for the sample.
"""
from __future__ import annotations


def sample_chemistry(sample: str) -> str:
    """Return chemistry string for a sample id."""
    return SAMPLES.at[sample, "chemistry"]


def is_10x(sample: str) -> bool:
    return sample_chemistry(sample) == "10x"


def is_ss(sample: str) -> bool:
    return sample_chemistry(sample) in {"smartseq2", "smartseq3"}


def cells_for_sample(sample: str):
    """All obs_id rows belonging to a sample."""
    return CELLS[CELLS["sample_id"] == sample]


def groups_for_sample(sample: str):
    """10x per-cell: list of unique cell_barcode values for the sample."""
    sub = cells_for_sample(sample)
    return sorted(sub["cell_barcode"].dropna().unique().tolist())


def get_obs_units(sample: str):
    """
    Per-chemistry list of observation units (the unit at which REDItools is
    invoked).
      - 10x: one unit per cell barcode (per-cell mode).
      - SS:  one unit per cell BAM.
    Returns a list of (obs_id, unit_name) tuples where unit_name is the
    wildcard used in result paths (``{group}`` for 10x, ``{cell}`` for SS).
    """
    sub = cells_for_sample(sample)
    if is_10x(sample):
        return [(f"{sample}__{bc}", bc) for bc in groups_for_sample(sample)]
    return [(row["obs_id"], row["bam_basename"]) for _, row in sub.iterrows()]


def reditools_obs_inputs(wildcards):
    """
    Aggregate-rule input function: returns the list of per-obs annotated parquet
    files for a given sample, choosing the right per-chemistry path layout.
    """
    sample = wildcards.sample
    out = []
    if is_10x(sample):
        for obs_id, group in get_obs_units(sample):
            for chrom in CHROMS:
                out.append(
                    f"results/filtered/{sample}/10x__{group}__{chrom}.annot.parquet"
                )
    else:
        for obs_id, cell in get_obs_units(sample):
            out.append(f"results/filtered/{sample}/ss__{cell}.annot.parquet")
    return out


def all_aggregate_inputs():
    """Top-level: every sample's annotated parquets, flat."""
    out = []
    for sample in SAMPLES.index:
        out.extend(
            reditools_obs_inputs(type("W", (), {"sample": sample})())
        )
    return out
