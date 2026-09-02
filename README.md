# OMIX GSEA Visualization Legacy

Generate publication-ready Gene Set Enrichment Analysis (GSEA) visualizations:

- enrichment score (ES) curves;
- ranked-gene (RNK) panels; and
- leading-edge (LE) expression heatmaps.

This repository is the **Code Ocean adapter** for the legacy OMIX GSEA
visualization implementation. It packages the run entry point, Code Ocean app
panel, and reproducible input contract. The plotting implementation is kept
intact while its inputs follow the current OMIX workflow handoff.

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
subdirectories below `/data`; the adapter discovers the required files
recursively.

The input bundle may also contain the DEG run summary and diagnostic images.
Those provenance files are preserved but ignored by the visualization step.

## Run locally

The equivalent command-line interface is:

```bash
Rscript code/main.R \
  --msigdb_database /path/to/MSigDB_v2023_2.rds \
  --gsea_filter_results /path/to/filtered_gsea_results.csv \
  --deg_analysis_results /path/to/deg-analysis-result \
  --plots_to_include ES+RNK+LE \
  --output_dir results
```

`--deg_analysis_results` must be the directory containing both portable DEG
output tables, not just `DEG_Analysis.csv` alone.

Alternatively, provide both `--deg_table` and `--sample_metadata` to override
the bundle with an explicit matched pair.

## Outputs

- `GSEA-Vis-Enrichment-Plots.pdf` — selected pathway plots in a multi-page
  PDF.
- `GSEA-Vis-RunningES.csv` — ranked-gene and running-enrichment-score values.

Run logs also report input-consistency checks, including whether the filtered
GSEA pathways and DEG ranking statistics are compatible.

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
