# sca2i-preprocess

Snakemake 8.x workflow that turns aligned RNA-seq BAMs (10x Chromium **or**
Smart-seq2/3) into a single AnnData object suitable for the
[scA2I](https://github.com/loganaq/scA2I) site-level A-to-I editing detector.

The workflow is intentionally **anchor-free** by default: REDItools2 is run in
de novo mode so REDIportal v3 can later be used as an *independent benchmark*
of the algorithm's filtering effect, not as an input filter.

---

## What it does

1. Prepares per-observation BAMs:
   - **10x**: `sinto filterbarcodes` splits the cellranger BAM into
     per-celltype/per-donor pseudobulk groups, then `umi_tools dedup`.
   - **Smart-seq2/3**: `picard MarkDuplicates` per per-cell BAM (cells already
     come as one BAM each).
2. Runs REDItools2 de novo per observation × chromosome (10x) or per cell
   (SS2/3). No `--bed_file` anchor.
3. Annotates each candidate site with:
   - ADAR sequence motif context (5' depleted G/A, 3' enriched G).
   - Alu repeat overlap (RepeatMasker).
4. Aggregates everything into a sparse AnnData:
   - `.X` = coverage `n` (CSR sparse).
   - `.layers["edits"]` = edit count `k`.
   - `.var` = site coordinates + motif/Alu annotations.
   - `.obs` = per-cell (SS2/3) or per-pseudobulk-group (10x) metadata.

The resulting `.h5ad` is the input contract for `sca2i.tl.discover()`.

## Quickstart

```bash
# 1. install snakemake + conda
mamba create -n snake -c bioconda -c conda-forge "snakemake>=8" "snakemake-wrapper-utils"
mamba activate snake

# 2. point config/samples.tsv and config/cells.tsv at your data
$EDITOR config/config.yaml config/samples.tsv config/cells.tsv

# 3. dry-run
make dryrun

# 4. run
make test JOBS=8
```

## Layout

```
sca2i-preprocess/
├── config/
│   ├── config.yaml            # workflow params (genome, REDItools, filters)
│   ├── samples.tsv            # one row per sample (chemistry, donor, BAM path)
│   └── cells.tsv              # one row per cell (SS2/3) or per group (10x)
├── workflow/
│   ├── Snakefile              # rule all + includes
│   ├── rules/
│   │   ├── common.smk         # helpers (chemistry dispatch, obs unit fan-out)
│   │   ├── refs.smk           # genome + Alu BED download/index
│   │   ├── prep_10x.smk       # sinto split + umi_tools dedup
│   │   ├── prep_ss.smk        # picard MarkDuplicates per cell
│   │   ├── reditools.smk      # REDItools2 invocations + parquet conversion
│   │   ├── filters.smk        # motif annotation, Alu intersect
│   │   └── aggregate.smk      # build final AnnData
│   ├── scripts/               # python helpers (one per rule that needs one)
│   ├── envs/                  # conda env YAMLs
│   └── schemas/               # JSON Schema for config validation
├── profiles/
│   ├── default/               # local execution
│   └── slurm/                 # cluster execution
├── Makefile                   # convenience targets
└── .github/workflows/ci.yml   # snakemake -n + schema lint on PR
```

## Outputs

```
results/
├── refs/
│   ├── genome.fa{,.fai}
│   └── alu.bed
├── prep/
│   ├── 10x/{sample}/groups/{group}.dedup.bam{,.bai}
│   └── ss/{sample}/{cell}.dedup.bam{,.bai}
├── reditools/
│   ├── 10x/{sample}/{group}/{chrom}.parquet
│   └── ss/{sample}/{cell}.parquet
├── filtered/
│   └── {sample}/{obs_unit}.annot.parquet     # adds motif + Alu cols
└── anndata/
    └── sca2i_input.h5ad                      # final deliverable
```

## Key design choices

| Choice | Why |
| --- | --- |
| De novo REDItools2 (no `--bed_file`) | Lets REDIportal v3 serve as an independent benchmark, not an input filter. Costs ~100x more compute but enables novel-site discovery + filter-effect quantification. |
| 10x → pseudobulk groups, SS2/3 → per cell | Matches scA2I's locked **A2** chemistry-aware decision: per-cell is too sparse at 10x coverage. |
| Per-chrom chunking for 10x | 10x BAMs can be huge; per-chrom REDItools fan-out parallelises trivially and bounds memory. |
| Sparse CSR AnnData | Site × cell matrices are >99% zero at scA2I LoD; dense would OOM. |
| Annotate, don't filter | Motif/Alu columns ride along; downstream `sca2i.tl.discover()` decides. |

## Compute estimates (rough)

| Stage | Per-sample wall (10x, 5k cells, 30k UMI/cell) | Per-cell wall (SS2/3, 1M reads/cell) |
| --- | --- | --- |
| sinto split | ~10 min | n/a |
| umi_tools dedup (per group) | 5–15 min | n/a |
| picard MarkDuplicates | n/a | ~30 s |
| REDItools2 de novo (per chrom × group, mid-coverage) | 20–60 min | 5–20 min |
| Aggregate AnnData | <5 min | <5 min |

Numbers assume `--threads 4` per REDItools call and ~24 chrom × ~6 groups for a typical PBMC 10x sample.

## Configuration

See `config/config.yaml` for all options. The schema at
`workflow/schemas/config.schema.yaml` is enforced by Snakemake at start-up.

Per-sample inputs live in `config/samples.tsv` (one row per sample/donor), and
per-observation metadata lives in `config/cells.tsv` (one row per Smart-seq cell
**or** per 10x pseudobulk group).

## License

MIT — see [LICENSE](LICENSE).
