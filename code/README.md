# OMIX GSEA Visualization Legacy

Generates enrichment plots from OMIX GSEA Filters Legacy, including:
- **ES (Enrichment Score)** plot with running score
- **RNK (Rank)** plot showing gene rankings
- **LE (Leading Edge)** heatmap with expression data

## Input Requirements

### Required
- `gsea_filter_result.rds` or `gsea_filter_result.csv`- GSEA results table filtered recomended
- `deg_table.rds` or `deg_table.csv` - DEG results (ranking statistics)
- `sample_metadata.rds` or `sample_metadata.csv` - sample id and grouping columns

### Optional
- `batch_result.rds` - NormalBatch-corrected expression data (required for LE heatmaps)

## Usage

### Command Line

```bash
# Basic usage with defaults (reads from /data)
./run

# Specify custom paths
./run /data/my_gsea_filter.rds /data/my_gsea_preranked.rds

# With batch correction for LE heatmaps
./run /data/gsea_filter.rds /data/gsea_preranked.rds /data/batch_result.rds

# Filter by specific contrasts
./run /data/gsea_filter.rds /data/gsea_preranked.rds "" "Treatment_vs_Control,Drug_vs_Vehicle"

# Control number of top pathways per contrast
./run /data/gsea_filter.rds /data/gsea_preranked.rds "" "" 3

# Choose plot types (ES, ES+RNK, ES+RNK+LE, LE)
./run /data/gsea_filter.rds /data/gsea_preranked.rds "" "" 1 "ES+RNK"
```

### Direct R Script Usage

```bash
Rscript /code/run_gsea_visualization.R \
  --gsea_filter_result /data/gsea_filter_result.rds \
  --sample_metadata.rds /data/sample_metadata.rds \
  --batch_result /data/batch_result.rds \
  --plot_contrasts "Treatment_vs_Control,Drug_vs_Vehicle" \
  --top_n_pathways 3 \
  --plots_to_include "ES+RNK+LE" \
  --max_plots_in_pdf 50 \
  --heatmap_transform "z-score" \
  --output_dir /results
```

### Available Options

```
--gsea_filter_result PATH        GSEA Filter result RDS file
--gsea_preranked_result PATH     GSEA Preranked result RDS file
--batch_result PATH              Batch Correction result RDS (optional)
--gsdb_result PATH               GSDB result RDS (optional)
--plot_contrasts "c1,c2,..."     Comma-separated contrasts to plot (default: all)
--plot_all_pathways              Plot all pathways (default: FALSE)
--top_n_pathways N               Top N pathways per contrast/collection (default: 1)
--preview_contrasts "c1,c2"      Generate preview PNGs for these contrasts
--max_plots_in_pdf N             Maximum plots in PDF (default: 50)
--plots_to_include TYPE          ES, ES+RNK, ES+RNK+LE, or LE (default: ES+RNK+LE)
--running_score_line_color       "ES sign" or "green" (default: ES sign)
--heatmap_transform TYPE         z-score, center by row mean, or none
--max_le_genes_heatmap N         Max genes in LE heatmap (default: 50)
--cluster_le_heatmap_rows        Cluster heatmap rows (default: TRUE)
--show_le_heatmap_gene_names     Show gene names in heatmap (default: TRUE)
--pdf_width W                    PDF width in inches (default: 8.5)
--pdf_height H                   PDF height in inches (default: 6.5)
--output_dir PATH                Output directory (default: /results)
```

## Outputs

All outputs are saved to `/results`:

- `GSEA-Vis-Enrichment-Plots-<timestamp>.pdf` - Multi-page PDF with all plots
- `GSEA-Vis-Preview-<timestamp>.png` - Preview image (first plot or specified contrast)
- `GSEA-Vis-Manifest-<timestamp>.csv` - Manifest of all generated plots
- `GSEA-Vis-RunningES-<timestamp>.csv` - Running enrichment score data
- `GSEA-Vis-Skipped-<timestamp>.csv` - List of skipped plots (if any)

## Examples

### Plot top 5 pathways for each contrast

```bash
./run /data/gsea_filter.rds /data/gsea_preranked.rds "" "" 5
```

### ES+RNK plots only (no heatmap)

```bash
./run /data/gsea_filter.rds /data/gsea_preranked.rds "" "" 1 "ES+RNK"
```

### Filter by specific contrasts

```bash
./run /data/gsea_filter.rds /data/gsea_preranked.rds "" "WT_vs_KO,Treated_vs_Control"
```

## Dependencies

- R packages: `ggplot2`, `patchwork`, `ComplexHeatmap`, `optparse`, `arrow`, `dplyr`
- The visualization functions are sourced from `GSEA_Visualization_Local_v1.R`

## Notes

- Default input paths assume GSEA Filter results are mounted at `/data`
- LE heatmaps require `batch_result` with expression matrix and sample metadata
- Set `--max_plots_in_pdf` to control memory usage for large result sets
- Use `--preview_contrasts` to generate standalone PNG files for specific contrasts
