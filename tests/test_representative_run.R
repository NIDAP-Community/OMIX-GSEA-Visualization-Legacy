#!/usr/bin/env Rscript

# End-to-end local adapter fixture. This is intentionally separate from the
# lightweight CI contract job because it requires the complete r-pathway
# plotting environment, including OmixPathwayPlots and ComplexHeatmap.

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."))

required_packages <- c(
  "optparse", "ggplot2", "patchwork", "ComplexHeatmap", "circlize",
  "OmixPathwayPlots"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Representative run requires the r-pathway environment; missing: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

source(file.path(root, "code", "functions", "gsea_enrichment_plot.R"))

requested_fixture_root <- Sys.getenv("OMIX_REPRESENTATIVE_FIXTURE_ROOT", unset = "")
fixture_root <- if (nzchar(requested_fixture_root)) {
  normalizePath(requested_fixture_root, mustWork = FALSE)
} else {
  tempfile("gsea-vis-representative-")
}
input_root <- file.path(fixture_root, "inputs")
output_root <- file.path(fixture_root, "results")
dir.create(input_root, recursive = TRUE)
if (!nzchar(requested_fixture_root)) {
  on.exit(unlink(fixture_root, recursive = TRUE, force = TRUE), add = TRUE)
}

genes <- sprintf("Gene%03d", seq_len(100L))
statistics <- setNames(seq(5, -5, length.out = length(genes)), genes)
pathway_genes <- genes[c(1:10, 31:35, 80:84)]
running <- gsea_vis_running_es_data(statistics, pathway_genes)
leading_edge <- gsea_vis_reconstructed_leading_edge_genes(running)
peak_es <- running$running_es[running$max_deviation][[1L]]

gsea <- data.frame(
  contrast = "B-A",
  pathways_database = "MSigDB",
  collection = "H",
  pathway = "HALLMARK_TEST_PATHWAY",
  pval = 0.001,
  padj = 0.01,
  ES = peak_es,
  NES = 1.8,
  size = length(pathway_genes),
  leadingEdge = paste(leading_edge, collapse = ","),
  size_leadingEdge = length(leading_edge),
  inPathway = paste(pathway_genes, collapse = ","),
  species = "Mouse",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

deg <- data.frame(
  GeneName = genes,
  `B-A_tstat` = unname(statistics),
  A1 = seq(4, 8, length.out = length(genes)),
  A2 = seq(4.2, 8.2, length.out = length(genes)),
  B1 = seq(5, 9, length.out = length(genes)),
  B2 = seq(5.2, 9.2, length.out = length(genes)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- data.frame(
  Sample = c("A1", "A2", "B1", "B2"),
  Group = c("A", "A", "B", "B"),
  stringsAsFactors = FALSE
)
msigdb <- data.frame(
  pathways_database = "MSigDB",
  species = "Mouse",
  collection = "H",
  gene_set_name = "HALLMARK_TEST_PATHWAY",
  gene_symbol = pathway_genes,
  stringsAsFactors = FALSE
)

paths <- list(
  gsea = file.path(input_root, "filtered_gsea_results.csv"),
  deg = file.path(input_root, "DEG_Analysis.csv"),
  metadata = file.path(input_root, "Sample_Metadata.csv"),
  msigdb = file.path(input_root, "MSigDB_test.csv")
)
write.csv(gsea, paths$gsea, row.names = FALSE)
write.csv(deg, paths$deg, row.names = FALSE)
write.csv(metadata, paths$metadata, row.names = FALSE)
write.csv(msigdb, paths$msigdb, row.names = FALSE)

command <- c(
  file.path(root, "code", "main.R"),
  "--msigdb_database", paths$msigdb,
  "--gsea_filter_results", paths$gsea,
  "--deg_table", paths$deg,
  "--sample_metadata", paths$metadata,
  "--output_dir", output_root,
  "--plots_to_include", "ES+RNK+LE",
  "--top_n_pathways", "20",
  "--pathway_bubble_plots", "TRUE",
  "--pathway_bubble_top_n", "20",
  "--pathway_bubble_significance_statistic", "padj"
)
status <- system2("Rscript", command, stdout = TRUE, stderr = TRUE)
exit_status <- attr(status, "status")
if (is.null(exit_status)) exit_status <- 0L
if (!identical(as.integer(exit_status), 0L)) {
  stop(paste(c("Representative adapter run failed:", status), collapse = "\n"), call. = FALSE)
}

expected <- c(
  "GSEA-Vis-Enrichment-Plots.pdf",
  "GSEA-Vis-RunningES.csv",
  "GSEA-Vis-Pathway-Bubble_manifest.csv"
)
missing <- expected[!file.exists(file.path(output_root, expected))]
if (length(missing) > 0L) {
  stop("Representative run did not create: ", paste(missing, collapse = ", "), call. = FALSE)
}

png_files <- list.files(
  output_root,
  pattern = "^GSEA-Vis-Pathway-Bubble.*\\.png$",
  full.names = TRUE
)
stopifnot(
  length(png_files) >= 3L,
  all(file.info(file.path(output_root, expected))$size > 100L),
  all(file.info(png_files)$size > 100L),
  nrow(read.csv(file.path(output_root, "GSEA-Vis-RunningES.csv"))) > 0L
)

# Verify the files really are PDF/PNG images rather than empty placeholders.
pdf_magic <- readBin(file.path(output_root, expected[[1L]]), "raw", n = 4L)
stopifnot(identical(rawToChar(pdf_magic), "%PDF"))
for (png_file in png_files) {
  png_magic <- readBin(png_file, "raw", n = 8L)
  stopifnot(identical(png_magic, as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))))
}

message("GSEA Visualization representative run and output checks passed: ", output_root)
