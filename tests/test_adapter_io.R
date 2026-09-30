#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."))
Sys.setenv(OMIX_ADAPTER_TEST_MODE = "1")
Sys.setenv(OMIX_ADAPTER_CODE_ROOT = file.path(root, "code"))
source(file.path(root, "code", "main.R"))

fixture_root <- tempfile("gsea-vis-adapter-")
dir.create(file.path(fixture_root, "msigdb", "release"), recursive = TRUE)
dir.create(file.path(fixture_root, "gsea_filter_results", "run"), recursive = TRUE)
dir.create(file.path(fixture_root, "deg-training", "run"), recursive = TRUE)
on.exit(unlink(fixture_root, recursive = TRUE, force = TRUE), add = TRUE)

msigdb <- file.path(fixture_root, "msigdb", "release", "MSigDB_test.csv")
gsea <- file.path(fixture_root, "gsea_filter_results", "run", "filtered_gsea_results.csv")
deg <- file.path(fixture_root, "deg-training", "run", "DEG_Analysis.csv")
metadata <- file.path(fixture_root, "deg-training", "run", "Sample_Metadata.csv")

write.csv(data.frame(collection = "H", gene_set_name = "HALLMARK_TEST", gene_symbol = "A"), msigdb, row.names = FALSE)
write.csv(data.frame(contrast = "B-A", collection = "H", pathway = "HALLMARK_TEST", inPathway = "A,B"), gsea, row.names = FALSE)
write.csv(data.frame(GeneName = c("A", "B"), `B-A_tstat` = c(2, -1), check.names = FALSE), deg, row.names = FALSE)
write.csv(data.frame(Sample = c("A1", "B1"), Group = c("A", "B")), metadata, row.names = FALSE)

stopifnot(
  identical(find_msigdb_file(file.path(fixture_root, "msigdb")), msigdb),
  identical(find_gsea_filter_file(file.path(fixture_root, "gsea_filter_results")), gsea)
)
bundle <- find_deg_analysis_bundle(file.path(fixture_root, "deg-training"))
stopifnot(
  identical(bundle$deg_table, deg),
  identical(bundle$sample_metadata, metadata)
)

partial_error <- tryCatch(
  { resolve_explicit_deg_inputs(deg, NULL); NULL },
  error = identity
)
stopifnot(
  inherits(partial_error, "error"),
  grepl("requires both", conditionMessage(partial_error), fixed = TRUE)
)

explicit <- resolve_explicit_deg_inputs(deg, metadata)
stopifnot(
  identical(explicit$deg_table, deg),
  identical(explicit$sample_metadata, metadata)
)

duplicate_dir <- file.path(fixture_root, "msigdb", "second")
dir.create(duplicate_dir)
file.copy(msigdb, file.path(duplicate_dir, "MSigDB_copy.csv"))
ambiguous <- tryCatch(
  { find_msigdb_file(file.path(fixture_root, "msigdb")); NULL },
  error = identity
)
stopifnot(
  inherits(ambiguous, "error"),
  grepl(msigdb, conditionMessage(ambiguous), fixed = TRUE),
  grepl("MSigDB_copy.csv", conditionMessage(ambiguous), fixed = TRUE)
)

message("GSEA Visualization adapter I/O checks passed")
