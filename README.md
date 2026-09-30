# OMIX GSEA Visualization Legacy

Generate publication-ready Gene Set Enrichment Analysis (GSEA) visualizations:

- enrichment score (ES) curves;
- ranked-gene (RNK) panels; and
- leading-edge (LE) expression heatmaps.

## Canonical OMIX module

| Item | Location |
| --- | --- |
| Canonical module | [OMIX GSEA Visualization Legacy](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-GSEA-Visualization-Legacy) |
| Interface contract | [schemas/interface.yml](https://github.com/NIDAP-Community/OMIX/blob/main/modules/OMIX-GSEA-Visualization-Legacy/schemas/interface.yml) |
| Development contract | [OMIX module contract](https://github.com/NIDAP-Community/OMIX/blob/main/docs/module-contract.md) |
| Released source reference | [OMIX_MODULE_SOURCE.md](OMIX_MODULE_SOURCE.md) |

The canonical module owns scientific behavior, portable CLI operation, tests,
and the reusable input/output contract. This repository is its Code Ocean
deployment adapter.

## What this deployment adds

- The Code Ocean App Panel and capsule entry point.
- Recursive workflow-result discovery and explicit DEG/metadata overrides.
- The input bundle and result layout used by the OMIX DEG-to-GSEA workflow.

## Workflow inputs

The capsule expects three inputs:

1. **MSigDB database** — the same MSigDB release used for upstream GSEA.
2. **Filtered GSEA results** — `filtered_gsea_results.csv` from OMIX GSEA
   Filters.
3. **DEG Analysis result bundle** — one OMIX DEG Analysis Result containing
   both `DEG_Analysis.csv` and `Sample_Metadata.csv` in the same folder.

The DEG table supplies the ranking statistics and sample-level expression. If
batch adjustment was included in DEG Analysis, its appended expression columns
are batch-adjusted voom-scale log-CPM values and are used directly for the LE
heatmap. The matching metadata controls sample annotation and ordering.

Keep `DEG_Analysis.csv` and `Sample_Metadata.csv` together. The adapter stops
on ambiguous or mismatched bundles instead of silently using files from
different analyses.

For an ad-hoc override, users may optionally provide both a DEG table and its
matching sample-metadata table in the app panel. The two explicit files take
precedence over the attached bundle and cannot be supplied independently.

## Use in Code Ocean

Attach the three inputs above, then choose the pathway selection and display
settings in the app panel. Workflow Results may mount in generated
subdirectories below each named data input. The adapter searches recursively
within `/data/msigdb`, `/data/gsea_filter_results`, and `/data/deg-training`
independently, so a file from one input cannot be selected for another.

The input bundle may also contain the DEG run summary and diagnostic images.
Those provenance files are preserved but ignored by the visualization step.

## Run locally

For an explicit local run, provide the four canonical files directly:

```bash
Rscript code/main.R \
  --msigdb_database /path/to/MSigDB_v2023_2.rds \
  --gsea_filter_results /path/to/filtered_gsea_results.csv \
  --deg_table /path/to/DEG_Analysis.csv \
  --sample_metadata /path/to/Sample_Metadata.csv \
  --plots_to_include ES+RNK+LE \
  --output_dir results
```

The Code Ocean adapter resolves the paired files from its combined attached
DEG result when the two explicit selectors are blank. Local runs should pass
both files explicitly; one without the other is rejected.

## Outputs

- `GSEA-Vis-Enrichment-Plots.pdf` — selected pathway plots in a multi-page
  PDF.
- `GSEA-Vis-RunningES.csv` — ranked-gene and running-enrichment-score values.
- `GSEA-Vis-Pathway-Bubble-*.png` plus a manifest — standardized bubble
  plots across all collections, within collections, and for each collection.
  These use adjusted p-values/FDR by default and show 20 pathways by default;
  the App Panel can switch to nominal p-values, another top-N limit, or shared
  color limits across collections.

Run logs also report input-consistency checks, including whether the filtered
GSEA pathways and DEG ranking statistics are compatible.

## Environment and reproducibility

Use the pinned capsule environment defined in `environment/`. Retain the
MSigDB release, filtered-GSEA result identity, DEG-result bundle, and selected
plot parameters with every released figure.

## Repository layout

```text
code/
  main.R                         Code Ocean and command-line entry point
  functions/gsea_enrichment_plot.R  Preserved visualization implementation
  README.md                      Detailed adapter documentation
.codeocean/                      App-panel and default-data configuration
environment/                     Capsule Docker environment
```

## Related OMIX repositories

- [OMIX DEG Analysis](https://github.com/NIDAP-Community/OMIX-DEG-Analysis)
- [OMIX GSEA Filters Legacy](https://github.com/NIDAP-Community/OMIX-GSEA-Filters-Legacy)
- [OMIX GSEA Preranked Legacy](https://github.com/NIDAP-Community/OMIX-GSEA-Preranked-Legacy)

“Legacy” identifies the established visualization implementation this adapter
preserves; it does not change the required input validation or provenance
checks.

## For developers

Read [AGENTS.md](AGENTS.md) and [OMIX_MODULE_SOURCE.md](OMIX_MODULE_SOURCE.md)
before editing. Reusable scientific changes belong in the canonical OMIX module
and are exported here only after validation.
