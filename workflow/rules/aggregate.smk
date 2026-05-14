"""Aggregate per-obs annotated parquets into a single sparse AnnData."""


def _aggregate_inputs(wildcards):
    """Collect every annotated parquet across all samples."""
    return all_aggregate_inputs()


rule aggregate_anndata:
    input:
        parquets   = _aggregate_inputs,
        cells_tsv  = config["cells_tsv"],
        samples_tsv = config["samples_tsv"],
    output:
        h5ad = config["aggregate"]["out_h5ad"],
    params:
        drop_empty = config["aggregate"]["drop_empty_sites"],
    log:
        "logs/aggregate/aggregate_anndata.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/aggregate_anndata.py"
