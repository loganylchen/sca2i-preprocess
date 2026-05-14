"""Shared helpers — chemistry dispatch, observation-unit fan-out."""
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
    """10x: list of celltype groups (one BAM per group after sinto split)."""
    sub = cells_for_sample(sample)
    return sorted(sub["celltype"].dropna().unique().tolist())


def get_obs_units(sample: str):
    """
    Per-chemistry list of observation units (the unit at which REDItools is
    invoked). For 10x this is celltype groups; for SS it is per-cell BAMs.
    Returns a list of (obs_id, unit_name) tuples.
    """
    sub = cells_for_sample(sample)
    if is_10x(sample):
        return [(f"{sample}__{g}", g) for g in groups_for_sample(sample)]
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
