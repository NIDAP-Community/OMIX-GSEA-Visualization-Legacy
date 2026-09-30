#!/usr/bin/env Rscript

# GSEA Visualization CLI - reads GSEA Filter results and creates enrichment plots
# Usage: Rscript main.R [options]

suppressPackageStartupMessages({
  library(optparse)
})

# Source the core visualization function from the functions directory.
# In Code Ocean, main.R is under /code and the core file is under /code/functions.
script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
test_code_root <- Sys.getenv("OMIX_ADAPTER_CODE_ROOT", unset = "")
script_dir <- if (nzchar(test_code_root)) {
  normalizePath(test_code_root, mustWork = TRUE)
} else if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE))
} else {
  getwd()
}
core_candidates <- c(
  file.path(script_dir, "functions", "gsea_enrichment_plot.R"),
  "/code/functions/gsea_enrichment_plot.R"
)
core_file <- core_candidates[file.exists(core_candidates)][1]
if (is.na(core_file)) {
  stop(
    "ERROR: gsea_enrichment_plot.R was not found under functions/ or /code/functions/.",
    call. = FALSE
  )
}
source(core_file)

# Validate that the CLI and core file are from the same interface generation.
required_core_functions <- c(
  "GSEA_Visualization_Local",
  "gsea_vis_bool",
  "gsea_vis_select_rows",
  "gsea_vis_authoritative_membership",
  "gsea_vis_validate_pathway_consistency"
)
missing_core_functions <- required_core_functions[!vapply(
  required_core_functions,
  exists,
  logical(1),
  mode = "function",
  inherits = TRUE
)]
required_core_parameters <- c(
  "selected_rows",
  "heatmap_gene_order",
  "heatmap_sample_order",
  "heatmap_gene_clustering_distance",
  "heatmap_gene_clustering_method",
  "heatmap_sample_clustering_distance",
  "heatmap_sample_clustering_method"
)
missing_core_parameters <- if (exists("GSEA_Visualization_Local", mode = "function")) {
  setdiff(required_core_parameters, names(formals(GSEA_Visualization_Local)))
} else {
  required_core_parameters
}
if (length(missing_core_functions) > 0 || length(missing_core_parameters) > 0) {
  stop(
    paste0(
      "ERROR: main.R is not compatible with the sourced core file: ",
      normalizePath(core_file, winslash = "/", mustWork = FALSE),
      ". Deploy the matching functions/gsea_enrichment_plot.R from the same package.",
      if (length(missing_core_functions) > 0) {
        paste0(" Missing functions: ", paste(missing_core_functions, collapse = ", "), ".")
      } else {
        ""
      },
      if (length(missing_core_parameters) > 0) {
        paste0(" Missing parameters: ", paste(missing_core_parameters, collapse = ", "), ".")
      } else {
        ""
      }
    ),
    call. = FALSE
  )
}

# CLI-local choice validator. Keeping this in the entry point avoids coupling
# command-line normalization to a helper defined by the plotting implementation.
gsea_cli_choice <- function(value, choices, default, parameter) {
  if (is.null(value) || length(value) == 0) {
    value <- default
  } else {
    value <- trimws(as.character(value[[1]]))
    if (!nzchar(value)) {
      value <- default
    }
  }
  if (!value %in% choices) {
    stop(
      sprintf(
        "Invalid %s: '%s'. Allowed values: %s.",
        parameter,
        value,
        paste(choices, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  value
}

#' Parse command-line arguments
get_args <- function() {
  option_list <- list(
    make_option(c("--msigdb_database"), type = "character", default = NULL,
                help = "Path to an MSigDB database (.rds or .csv). Used to restore pathway membership when it is absent from filtered GSEA results."),
    make_option(c("--gsea_filter_results"), type = "character", default = NULL,
                help = "Path to GSEA Filter results file (CSV or RDS flat table)"),
    make_option(c("--deg_table"), type = "character", default = NULL,
                help = "Optional explicit DEG Analysis table (CSV or RDS). Must be supplied together with --sample_metadata and overrides the attached DEG Analysis result bundle."),
    make_option(c("--sample_metadata"), type = "character", default = NULL,
                help = "Optional explicit sample metadata table (CSV or RDS). Must be supplied together with --deg_table and overrides the attached DEG Analysis result bundle."),
    make_option(c("--output_dir"), type = "character", default = "/results",
                help = "Output directory [default: /results]"),
    make_option(c("--contrast_filter"), type = "character", default = "none",
                help = "Contrast filter mode: none / keep / remove [default: none]"),
    make_option(c("--contrasts"), type = "character", default = NULL,
                help = "Comma-separated contrast names for keep/remove mode"),
    make_option(c("--top_n_pathways"), type = "integer", default = 20L,
                help = "Top pathways per contrast/collection. 0 = all pathways. [default: 20]"),
    make_option(c("--top_n_by_sign"), type = "logical", default = FALSE,
                help = "When TRUE, select top_n_pathways separately for positively (ES>0) and negatively (ES<0) enriched pathways. [default: FALSE]"),
    make_option(c("--pathway_bubble_plots"), type = "logical", default = TRUE,
                help = "Write shared OMIX pathway bubble plots in addition to legacy enrichment panels [default: TRUE]"),
    make_option(c("--pathway_bubble_top_n"), type = "integer", default = 20L,
                help = "Top pathways in each shared bubble-plot selection. 0 = all pathways. [default: 20]"),
    make_option(c("--pathway_bubble_significance_statistic"), type = "character", default = "padj",
                help = "padj (FDR) or pval used for shared bubble-plot selection [default: padj]"),
    make_option(c("--collection_color_scale"), type = "character", default = "independent",
                help = "independent or shared collection-specific bubble-plot colour scales [default: independent]"),
    make_option(c("--max_plots_in_pdf"), type = "integer", default = 0L,
                help = "Total plots in PDF (global ceiling). 0 = no limit. [default: 0]"),
    make_option(c("--plots_to_include"), type = "character", default = "ES+RNK+LE",
                help = "Plot types: ES, ES+RNK, ES+LE, ES+RNK+LE, LE [default: ES+RNK+LE]"),
    make_option(c("--running_score_line_color"), type = "character", default = "ES sign",
                help = "Line color: 'ES sign' or 'green' [default: ES sign]"),
    make_option(c("--add_max_deviation_line"), type = "character", default = "both",
                help = "Max deviation guide line: 'coordinate', 'horizontal', 'both', or 'none' [default: both]"),
    make_option(c("--rank_area_color"), type = "character", default = "red/blue by Gene score",
                help = "RNK area fill color: 'grey' or 'red/blue by Gene score' [default: red/blue by Gene score]"),
    make_option(c("--show_es_rank_bar"), type = "logical", default = FALSE,
                help = "Show colored rank bar strip and +/- signs in ES panel [default: FALSE]"),
    make_option(c("--show_rnk_peak_line"), type = "logical", default = TRUE,
                help = "Show dashed vertical line at peak ES rank in RNK panel [default: TRUE]"),
    make_option(c("--show_rnk_le_highlight"), type = "logical", default = TRUE,
                help = "Show leading-edge highlight rectangle in RNK panel [default: TRUE]"),
    make_option(c("--show_es_le_highlight"), type = "logical", default = TRUE,
                help = "Show leading-edge highlight rectangle in ES panel [default: TRUE]"),
    make_option(c("--heatmap_transform"), type = "character", default = "z-score",
                help = "Heatmap transform: 'z-score', 'center by row mean', 'center by row median', or 'none' [default: z-score]"),
    make_option(c("--max_le_genes_heatmap"), type = "integer", default = 50,
                help = "Max genes in LE heatmap [default: 50]"),
    make_option(c("--heatmap_gene_order"), type = "character", default = "rank",
                help = "Gene ordering: rank / cluster / input [default: rank]"),
    make_option(c("--heatmap_sample_order"), type = "character", default = "group",
                help = "Sample ordering: group / cluster / input [default: group]"),
    make_option(c("--heatmap_gene_clustering_distance"), type = "character", default = "euclidean",
                help = "Gene clustering distance when gene order is cluster [default: euclidean]"),
    make_option(c("--heatmap_gene_clustering_method"), type = "character", default = "complete",
                help = "Gene hierarchical clustering method [default: complete]"),
    make_option(c("--heatmap_sample_clustering_distance"), type = "character", default = "euclidean",
                help = "Sample clustering distance when sample order is cluster [default: euclidean]"),
    make_option(c("--heatmap_sample_clustering_method"), type = "character", default = "complete",
                help = "Sample hierarchical clustering method [default: complete]"),
    make_option(c("--show_le_heatmap_gene_names"), type = "logical", default = TRUE,
                help = "Show gene names in heatmap [default: TRUE]"),
    make_option(c("--show_le_heatmap_sample_names"), type = "logical", default = FALSE,
                help = "Show sample names in heatmap [default: FALSE]"),
    make_option(c("--show_le_heatmap_rank_labels"), type = "logical", default = TRUE,
                help = "Prefix heatmap gene labels with compact leading-edge ranks, e.g. 1. MYC [default: TRUE]"),
    make_option(c("--heatmap_gene_names_column"), type = "character", default = "GeneName",
                help = "Gene column name in expression table [default: GeneName]"),
    make_option(c("--heatmap_sample_names_column"), type = "character", default = "Sample",
                help = "Sample column name in metadata [default: Sample]"),
    make_option(c("--heatmap_group_column"), type = "character", default = "Group",
                help = "Group column name in metadata [default: Group]"),
    make_option(c("--pdf_width"), type = "numeric", default = 8.5,
                help = "PDF width in inches [default: 8.5]"),
    make_option(c("--pdf_height"), type = "numeric", default = 6.5,
                help = "PDF height in inches [default: 6.5]")
  )
  
  parser <- OptionParser(
    usage = "%prog [options]",
    description = "Generate GSEA enrichment plots from GSEA Filter results",
    option_list = option_list
  )
  
  optparse::parse_args(parser)
}

#' Parse comma-separated string to vector
parse_comma_list <- function(value) {
  if (is.null(value) || length(value) == 0 || !nzchar(trimws(value[[1]]))) {
    return(character(0))
  }

  values <- trimws(unlist(strsplit(as.character(value[[1]]), ",", fixed = TRUE)))
  unique(values[!is.na(values) & nzchar(values)])
}

#' Treat empty or whitespace-only strings as NULL
non_empty <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NULL)
  }

  value <- trimws(as.character(value[[1]]))
  if (nzchar(value)) value else NULL
}

#' Read a CSV or RDS tabular input with a caller-specific CSV reader
read_tabular_input <- function(path, description, csv_reader) {
  if (is.null(path) || !file.exists(path)) {
    stop(sprintf("ERROR: %s not found at: %s", description, path %||% "<not provided>"), call. = FALSE)
  }

  cat(sprintf("Loading %s from: %s\n", description, path))
  extension <- tolower(tools::file_ext(path))
  value <- tryCatch(
    {
      if (identical(extension, "rds")) {
        readRDS(path)
      } else if (identical(extension, "csv")) {
        csv_reader(path)
      } else {
        stop(sprintf("Unsupported file format: %s (expected .rds or .csv)", extension))
      }
    },
    error = function(error) {
      stop(sprintf("ERROR: Failed to load %s: %s", description, conditionMessage(error)), call. = FALSE)
    }
  )

  if (is.matrix(value)) {
    value <- as.data.frame(value, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(value)) {
    stop(sprintf("ERROR: %s must contain a data frame or matrix.", description), call. = FALSE)
  }

  value
}

#' Normalize an unnamed first identifier column from CSV exports
normalize_first_identifier_column <- function(data, default_name) {
  if (!is.data.frame(data) || ncol(data) == 0) {
    return(data)
  }

  first_name <- names(data)[[1]]
  first_values <- data[[1]]
  unnamed <- is.na(first_name) || !nzchar(first_name) || first_name %in% c("...1", "row.names")
  identifier_like <- is.character(first_values) || is.factor(first_values)
  if (unnamed && identifier_like) {
    names(data)[[1]] <- default_name
  }

  data
}

#' Resolve an explicit DEG table plus metadata override as one logical bundle.
#'
#' Both files are required together. This deliberately avoids mixing a
#' user-uploaded DEG table with metadata discovered from another workflow
#' Result/Data Asset.
resolve_explicit_deg_inputs <- function(deg_table_path, sample_metadata_path) {
  supplied <- c(
    deg_table = !is.null(deg_table_path),
    sample_metadata = !is.null(sample_metadata_path)
  )
  if (!any(supplied)) {
    return(NULL)
  }
  if (!all(supplied)) {
    missing <- names(supplied)[!supplied]
    stop(
      paste0(
        "ERROR: Explicit DEG input requires both --deg_table and --sample_metadata. Missing: ",
        paste(missing, collapse = ", "), ". Leave both blank to use the attached DEG Analysis result bundle."
      ),
      call. = FALSE
    )
  }

  for (item in c(deg_table = deg_table_path, sample_metadata = sample_metadata_path)) {
    if (!file.exists(item)) {
      stop(sprintf("ERROR: Explicit input file not found: %s", item), call. = FALSE)
    }
  }

  list(
    directory = "<explicit file override>",
    deg_table = deg_table_path,
    sample_metadata = sample_metadata_path
  )
}

#' Read the flat filtered GSEA result without dropping its first column
read_gsea_filter <- function(path) {
  data <- read_tabular_input(
    path,
    "GSEA Filter results",
    function(csv_path) {
      utils::read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  )

  required <- c("contrast", "collection", "pathway")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      sprintf("ERROR: GSEA Filter results are missing required column(s): %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  data
}

#' Read DEG Analysis sample-level expression without assuming the first column is row names
read_expression_table <- function(path) {
  data <- read_tabular_input(
    path,
    "DEG Analysis sample-level expression",
    function(csv_path) {
      utils::read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  )
  data <- normalize_first_identifier_column(data, "gene")

  known_gene_columns <- c("gene", "Gene", "gene_name", "GeneName", "symbol", "Symbol")
  if (!any(known_gene_columns %in% names(data))) {
    row_ids <- rownames(data)
    meaningful_row_ids <-
      !is.null(row_ids) && length(row_ids) == nrow(data) && nrow(data) > 0 &&
      all(!is.na(row_ids) & nzchar(row_ids)) &&
      !all(row_ids == as.character(seq_len(nrow(data))))
    if (meaningful_row_ids) {
      data <- data.frame(gene = row_ids, data, check.names = FALSE, stringsAsFactors = FALSE)
      rownames(data) <- NULL
      cat("  Added 'gene' column from row names for heatmap compatibility.\n")
    }
  }

  data
}

#' Read sample metadata without dropping its sample identifier column
read_sample_metadata <- function(path) {
  data <- read_tabular_input(
    path,
    "sample metadata",
    function(csv_path) {
      utils::read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  )
  data <- normalize_first_identifier_column(data, "Sample")

  known_sample_columns <- c("Sample", "sample", "sample_id", "SampleID")
  if (!any(known_sample_columns %in% names(data))) {
    row_ids <- rownames(data)
    meaningful_row_ids <-
      !is.null(row_ids) && length(row_ids) == nrow(data) && nrow(data) > 0 &&
      all(!is.na(row_ids) & nzchar(row_ids)) &&
      !all(row_ids == as.character(seq_len(nrow(data))))
    if (meaningful_row_ids) {
      data <- data.frame(Sample = row_ids, data, check.names = FALSE, stringsAsFactors = FALSE)
      rownames(data) <- NULL
      cat("  Added 'Sample' column from row names for metadata compatibility.\n")
    }
  }

  data
}

#' Construct the compact batch-result structure used by the heatmap code
construct_batch_result <- function(expression_path, metadata_path) {
  if (is.null(expression_path) || is.null(metadata_path)) {
    return(NULL)
  }

  expression_data <- read_expression_table(expression_path)
  sample_metadata <- read_sample_metadata(metadata_path)

  cat(sprintf("  Expression table: %d genes x %d data columns\n", nrow(expression_data), ncol(expression_data)))
  cat(sprintf("  Sample metadata: %d samples x %d columns\n", nrow(sample_metadata), ncol(sample_metadata)))

  list(
    final_expression = expression_data,
    metadata = sample_metadata
  )
}

#' List matching files below a data root.
#'
#' Workflow-connected Results are mounted in a generated subdirectory of
#' /data, so recursive discovery is required rather than hard-coded mounts.
find_data_files <- function(data_root = "/data", pattern) {
  if (!dir.exists(data_root)) {
    return(character(0))
  }
  sort(list.files(
    data_root,
    pattern = pattern,
    ignore.case = TRUE,
    full.names = TRUE,
    recursive = TRUE
  ))
}

#' Require exactly one matching input file.
find_unique_data_file <- function(label, pattern, data_root = "/data") {
  candidates <- find_data_files(data_root = data_root, pattern = pattern)
  if (length(candidates) == 0L) {
    return(NULL)
  }
  if (length(candidates) > 1L) {
    stop(
      sprintf(
        "ERROR: Multiple %s files were found below %s. Attach exactly one matching input: %s",
        label,
        data_root,
        paste(candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  candidates[[1]]
}

#' Resolve the filtered GSEA result.  CSV is preferred for the portable
#' workflow contract; RDS remains accepted for historical local runs.
find_gsea_filter_file <- function(data_root = "/data") {
  csv_file <- find_unique_data_file(
    label = "filtered GSEA results CSV",
    pattern = "^filtered_gsea_results\\.csv$",
    data_root = data_root
  )
  if (!is.null(csv_file)) {
    return(csv_file)
  }
  find_unique_data_file(
    label = "filtered GSEA results RDS",
    pattern = "^filtered_gsea_results\\.rds$",
    data_root = data_root
  )
}

#' Resolve one MSigDB database supplied as its own workflow input.
find_msigdb_file <- function(data_root = "/data") {
  find_unique_data_file(
    label = "MSigDB database",
    pattern = "msigdb.*\\.(rds|csv)$",
    data_root = data_root
  )
}

#' Locate the two interoperable files in one OMIX DEG Analysis result bundle.
#'
#' The bundle can also include diagnostics and run_summary.txt; those are
#' intentionally ignored.  Only the two portable handoff tables are used.
find_deg_analysis_bundle <- function(data_root = "/data") {
  deg_tables <- find_data_files(
    data_root = data_root,
    pattern = "^DEG_Analysis\\.csv$"
  )
  metadata_tables <- find_data_files(
    data_root = data_root,
    pattern = "^Sample_Metadata\\.csv$"
  )

  if (length(deg_tables) == 0L || length(metadata_tables) == 0L) {
    missing <- c(
      if (length(deg_tables) == 0L) "DEG_Analysis.csv",
      if (length(metadata_tables) == 0L) "Sample_Metadata.csv"
    )
    stop(
      paste0(
        "ERROR: The DEG Analysis input bundle must contain ",
        paste(missing, collapse = " and "),
        ". Attach one OMIX DEG Analysis Result with both files together."
      ),
      call. = FALSE
    )
  }

  deg_dirs <- dirname(deg_tables)
  metadata_dirs <- dirname(metadata_tables)
  bundle_dirs <- intersect(deg_dirs, metadata_dirs)
  if (length(bundle_dirs) == 0L) {
    stop(
      paste0(
        "ERROR: DEG_Analysis.csv and Sample_Metadata.csv were found, but not in the same result bundle. ",
        "Attach one combined OMIX DEG Analysis Result."
      ),
      call. = FALSE
    )
  }
  if (length(bundle_dirs) > 1L) {
    stop(
      sprintf(
        "ERROR: Multiple OMIX DEG Analysis result bundles were found below %s: %s",
        data_root,
        paste(bundle_dirs, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  bundle_dir <- bundle_dirs[[1]]
  list(
    directory = bundle_dir,
    deg_table = deg_tables[match(bundle_dir, deg_dirs)],
    sample_metadata = metadata_tables[match(bundle_dir, metadata_dirs)]
  )
}

#' Read and validate the portable five-column MSigDB database table.
read_msigdb_database <- function(path) {
  database <- read_tabular_input(
    path,
    "MSigDB database",
    function(csv_path) {
      utils::read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  )
  required <- c("collection", "gene_set_name", "gene_symbol")
  missing <- setdiff(required, names(database))
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "ERROR: MSigDB database is missing required column(s): %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  database
}

#' Return TRUE when a GSEA row already carries a usable pathway-membership list.
has_pathway_membership <- function(row) {
  fields <- intersect(c("inPathway_orthologs", "inPathway"), names(row))
  any(vapply(fields, function(field) {
    length(gsea_vis_gene_vector(row[[field]])) > 0L
  }, logical(1)))
}

#' Restore absent memberships from the attached MSigDB database.
#'
#' GSEA Filter normally carries inPathway already.  Keeping it intact retains
#' analysis provenance and avoids loading the 5M-row bundled database unless
#' the original legacy fallback is actually needed.
restore_msigdb_membership <- function(gsea_rows, msigdb_path) {
  missing_membership <- !vapply(
    seq_len(nrow(gsea_rows)),
    function(index) has_pathway_membership(gsea_rows[index, , drop = FALSE]),
    logical(1)
  )
  if (!any(missing_membership)) {
    cat("  MSigDB database: attached (embedded pathway membership already present)\n")
    return(gsea_rows)
  }

  database <- read_msigdb_database(msigdb_path)
  selected_rows <- which(missing_membership)
  selected_collections <- unique(as.character(gsea_rows$collection[selected_rows]))
  selected_pathways <- unique(as.character(gsea_rows$pathway[selected_rows]))
  database <- database[
    as.character(database$collection) %in% selected_collections &
      as.character(database$gene_set_name) %in% selected_pathways,
    , drop = FALSE
  ]
  if (nrow(database) == 0L) {
    stop(
      "ERROR: The supplied MSigDB database has no records for the selected filtered GSEA pathways.",
      call. = FALSE
    )
  }

  matches_value <- function(values, target) {
    tolower(trimws(as.character(values))) == tolower(trimws(as.character(target)))
  }
  if (!"inPathway" %in% names(gsea_rows)) {
    gsea_rows$inPathway <- NA_character_
  }
  for (index in selected_rows) {
    row <- gsea_rows[index, , drop = FALSE]
    candidates <- database[
      matches_value(database$collection, row$collection[[1]]) &
        matches_value(database$gene_set_name, row$pathway[[1]]),
      , drop = FALSE
    ]
    if ("pathways_database" %in% names(row) && "pathways_database" %in% names(candidates) &&
        !is.na(row$pathways_database[[1]]) && nzchar(as.character(row$pathways_database[[1]]))) {
      candidates <- candidates[matches_value(candidates$pathways_database, row$pathways_database[[1]]), , drop = FALSE]
    }
    if ("species" %in% names(row) && "species" %in% names(candidates) &&
        !is.na(row$species[[1]]) && nzchar(as.character(row$species[[1]]))) {
      candidates <- candidates[matches_value(candidates$species, row$species[[1]]), , drop = FALSE]
    }
    genes <- unique(as.character(candidates$gene_symbol))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0L) {
      stop(
        sprintf(
          "ERROR: MSigDB did not provide genes for %s / %s. Check that database version, species, collection, and pathway name match the filtered GSEA result.",
          row$collection[[1]],
          row$pathway[[1]]
        ),
        call. = FALSE
      )
    }
    gsea_rows$inPathway[[index]] <- paste(genes, collapse = ",")
  }
  cat(sprintf("  MSigDB database: restored pathway membership for %d selected pathway(s)\n", length(selected_rows)))
  gsea_rows
}

#' Detect the DEG gene identifier column from column names
find_deg_gene_column <- function(column_names) {
  candidates <- intersect(
    c("GeneName", "Gene", "gene", "gene_name", "gene_id", "gene_symbol", "symbol", "SYMBOL", "GENE"),
    column_names
  )
  if (length(candidates) > 0) candidates[[1]] else column_names[[1]]
}

#' Detect preferred wide-format DEG ranking columns
find_wide_deg_columns <- function(column_names) {
  preferred_suffixes <- c("_tstat", "_logFC", "_log2FC", "_logFoldChange", "_FC")
  for (suffix in preferred_suffixes) {
    hits <- grep(paste0(suffix, "$"), column_names, value = TRUE, ignore.case = TRUE)
    if (length(hits) > 0) {
      return(list(suffix = suffix, columns = hits))
    }
  }
  NULL
}

#' Read only DEG columns needed for the selected contrasts when the input is CSV
read_deg_columns <- function(path, requested_contrasts) {
  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "rds")) {
    value <- readRDS(path)
    if (!is.data.frame(value)) {
      stop("ERROR: DEG RDS must contain a data frame.", call. = FALSE)
    }
    return(value)
  }
  if (!identical(extension, "csv")) {
    stop(sprintf("ERROR: Unsupported DEG format: %s", extension), call. = FALSE)
  }

  header <- utils::read.csv(path, nrows = 0, check.names = FALSE, stringsAsFactors = FALSE)
  column_names <- names(header)
  if (length(column_names) == 0) {
    stop("ERROR: DEG CSV has no columns.", call. = FALSE)
  }

  gene_column <- find_deg_gene_column(column_names)
  wide <- find_wide_deg_columns(column_names)
  if (!is.null(wide)) {
    contrast_names <- sub(paste0(wide$suffix, "$"), "", wide$columns, ignore.case = TRUE)
    keep_wide <- if (length(requested_contrasts) > 0) {
      wide$columns[contrast_names %in% requested_contrasts]
    } else {
      wide$columns
    }
    if (length(keep_wide) == 0) {
      stop(
        sprintf(
          "ERROR: None of the selected contrast(s) were found in the wide DEG table: %s",
          paste(requested_contrasts, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    keep_columns <- unique(c(gene_column, keep_wide))
  } else {
    rank_column <- intersect(
      c("stat", "statistic", "log2FoldChange", "log2FC", "logFC", "LFC", "t", "score"),
      column_names
    )
    if (length(rank_column) == 0) {
      stop(
        "ERROR: DEG table has no recognized ranking column and no wide-format ranking columns.",
        call. = FALSE
      )
    }
    contrast_column <- intersect(
      c("contrast", "condition", "comparison", "group", "Contrast"),
      column_names
    )
    contrast_keep <- if (length(contrast_column) > 0) contrast_column[[1]] else character(0)
    keep_columns <- unique(c(gene_column, rank_column[[1]], contrast_keep))
  }

  column_classes <- rep("NULL", length(column_names))
  column_classes[match(keep_columns, column_names)] <- NA
  data <- utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = column_classes
  )
  normalize_first_identifier_column(data, "gene")
}

#' Build ranked statistics only for contrasts that can be plotted
load_deg_table_as_ranked_stats <- function(path, requested_contrasts = character(0)) {
  if (is.null(path) || !file.exists(path)) {
    stop(sprintf("ERROR: DEG table not found at: %s", path %||% "<not provided>"), call. = FALSE)
  }

  requested_contrasts <- unique(as.character(requested_contrasts))
  requested_contrasts <- requested_contrasts[!is.na(requested_contrasts) & nzchar(requested_contrasts)]
  cat(sprintf("Loading DEG table from: %s\n", path))

  deg <- tryCatch(
    read_deg_columns(path, requested_contrasts),
    error = function(error) {
      stop(sprintf("ERROR: Failed to load DEG table: %s", conditionMessage(error)), call. = FALSE)
    }
  )
  if (!is.data.frame(deg) || ncol(deg) == 0) {
    stop("ERROR: DEG table must be a non-empty data frame.", call. = FALSE)
  }

  known_gene_columns <- intersect(
    c("Gene", "gene", "gene_name", "gene_id", "gene_symbol", "symbol", "SYMBOL", "GENE"),
    names(deg)
  )
  if (length(known_gene_columns) > 0) {
    gene_column <- known_gene_columns[[1]]
  } else {
    row_ids <- rownames(deg)
    meaningful_row_ids <-
      !is.null(row_ids) && length(row_ids) == nrow(deg) && nrow(deg) > 0 &&
      all(!is.na(row_ids) & nzchar(row_ids)) &&
      !all(row_ids == as.character(seq_len(nrow(deg))))
    if (meaningful_row_ids) {
      deg <- data.frame(gene = row_ids, deg, check.names = FALSE, stringsAsFactors = FALSE)
      rownames(deg) <- NULL
      gene_column <- "gene"
      cat("  Added 'gene' column from DEG row names.\n")
    } else {
      gene_column <- names(deg)[[1]]
    }
  }
  cat(sprintf("  Gene column     : %s\n", gene_column))

  wide <- find_wide_deg_columns(names(deg))
  if (!is.null(wide)) {
    contrast_names <- sub(paste0(wide$suffix, "$"), "", wide$columns, ignore.case = TRUE)
    if (length(requested_contrasts) > 0) {
      keep <- contrast_names %in% requested_contrasts
      wide$columns <- wide$columns[keep]
      contrast_names <- contrast_names[keep]
    }
    if (length(wide$columns) == 0) {
      stop("ERROR: DEG table contains no ranking columns for the selected contrasts.", call. = FALSE)
    }

    genes <- as.character(deg[[gene_column]])
    ranked_stats <- setNames(
      lapply(wide$columns, function(column_name) {
        values <- suppressWarnings(as.numeric(deg[[column_name]]))
        statistics <- setNames(values, genes)
        statistics <- statistics[
          !is.na(statistics) & !duplicated(names(statistics)) &
            !is.na(names(statistics)) & nzchar(names(statistics))
        ]
        sort(statistics, decreasing = TRUE)
      }),
      contrast_names
    )
    cat(sprintf("  Layout          : wide (ranking = *%s)\n", wide$suffix))
  } else {
    rank_column <- intersect(
      c("stat", "statistic", "log2FoldChange", "log2FC", "logFC", "LFC", "t", "score"),
      names(deg)
    )
    if (length(rank_column) == 0) {
      stop(
        "ERROR: DEG table has no recognized ranking column (stat, log2FoldChange, logFC) and no wide-format ranking columns.",
        call. = FALSE
      )
    }
    rank_column <- rank_column[[1]]
    contrast_column <- intersect(
      c("contrast", "condition", "comparison", "group", "Contrast"),
      names(deg)
    )

    genes <- as.character(deg[[gene_column]])
    statistics <- suppressWarnings(as.numeric(deg[[rank_column]]))
    contrast_values <- if (length(contrast_column) > 0) {
      cat(sprintf("  Contrast column : %s\n", contrast_column[[1]]))
      as.character(deg[[contrast_column[[1]]]])
    } else {
      rep("default", nrow(deg))
    }

    valid <-
      !is.na(genes) & nzchar(genes) & !is.na(statistics) &
      !is.na(contrast_values) & nzchar(contrast_values)
    if (length(requested_contrasts) > 0) {
      valid <- valid & contrast_values %in% requested_contrasts
    }
    genes <- genes[valid]
    statistics <- statistics[valid]
    contrast_values <- contrast_values[valid]
    if (length(genes) == 0) {
      stop("ERROR: DEG table produced no valid rows for the selected contrasts.", call. = FALSE)
    }

    ranked_stats <- lapply(
      split(seq_along(genes), contrast_values),
      function(index) {
        values <- setNames(statistics[index], genes[index])
        values <- values[!duplicated(names(values))]
        sort(values, decreasing = TRUE)
      }
    )
    cat(sprintf("  Layout          : long (ranking = %s)\n", rank_column))
  }

  ranked_stats <- ranked_stats[vapply(ranked_stats, length, integer(1)) > 0]
  if (length(ranked_stats) == 0) {
    stop("ERROR: No ranked statistics were generated for the selected contrasts.", call. = FALSE)
  }
  cat(sprintf("  Contrasts loaded: %s\n", paste(names(ranked_stats), collapse = ", ")))
  cat(sprintf("  Genes (first)   : %d\n", length(ranked_stats[[1]])))
  ranked_stats
}

#' Normalize and validate CLI parameters once
normalize_args <- function(args) {
  file_parameters <- c(
    "msigdb_database", "gsea_filter_results", "deg_table", "sample_metadata"
  )
  for (parameter in file_parameters) {
    args[[parameter]] <- non_empty(args[[parameter]])
  }

  args$contrast_filter <- gsea_cli_choice(
    non_empty(args$contrast_filter) %||% "none",
    c("none", "keep", "remove"),
    "none",
    "contrast_filter"
  )
  args$contrast_names <- parse_comma_list(args$contrasts)
  if (identical(args$contrast_filter, "keep") && length(args$contrast_names) == 0) {
    stop("ERROR: Contrast filter mode 'keep' requires at least one contrast name.", call. = FALSE)
  }

  args$plots_to_include <- gsea_cli_choice(
    non_empty(args$plots_to_include) %||% "ES+RNK+LE",
    c("ES", "ES+RNK", "ES+LE", "ES+RNK+LE", "LE"),
    "ES+RNK+LE",
    "plots_to_include"
  )
  args$heatmap_transform <- gsea_cli_choice(
    non_empty(args$heatmap_transform) %||% "z-score",
    c("z-score", "center by row mean", "center by row median", "none"),
    "z-score",
    "heatmap_transform"
  )
  args$heatmap_gene_order <- gsea_cli_choice(
    non_empty(args$heatmap_gene_order) %||% "rank",
    c("rank", "cluster", "input"),
    "rank",
    "heatmap_gene_order"
  )
  args$heatmap_sample_order <- gsea_cli_choice(
    non_empty(args$heatmap_sample_order) %||% "group",
    c("group", "cluster", "input"),
    "group",
    "heatmap_sample_order"
  )

  distance_choices <- c(
    "euclidean", "maximum", "manhattan", "canberra", "binary",
    "minkowski", "pearson", "spearman", "kendall"
  )
  method_choices <- c(
    "ward.D", "ward.D2", "single", "complete", "average",
    "mcquitty", "median", "centroid"
  )
  args$heatmap_gene_clustering_distance <- gsea_cli_choice(
    non_empty(args$heatmap_gene_clustering_distance) %||% "euclidean",
    distance_choices,
    "euclidean",
    "heatmap_gene_clustering_distance"
  )
  args$heatmap_gene_clustering_method <- gsea_cli_choice(
    non_empty(args$heatmap_gene_clustering_method) %||% "complete",
    method_choices,
    "complete",
    "heatmap_gene_clustering_method"
  )
  args$heatmap_sample_clustering_distance <- gsea_cli_choice(
    non_empty(args$heatmap_sample_clustering_distance) %||% "euclidean",
    distance_choices,
    "euclidean",
    "heatmap_sample_clustering_distance"
  )
  args$heatmap_sample_clustering_method <- gsea_cli_choice(
    non_empty(args$heatmap_sample_clustering_method) %||% "complete",
    method_choices,
    "complete",
    "heatmap_sample_clustering_method"
  )

  args$heatmap_gene_names_column <- non_empty(args$heatmap_gene_names_column) %||% "GeneName"
  args$heatmap_sample_names_column <- non_empty(args$heatmap_sample_names_column) %||% "Sample"
  args$heatmap_group_column <- non_empty(args$heatmap_group_column) %||% "Group"
  args$rank_area_color <- non_empty(args$rank_area_color) %||% "red/blue by Gene score"
  args$output_dir <- non_empty(args$output_dir) %||% "/results"

  args$pdf_width <- suppressWarnings(as.numeric(args$pdf_width %||% 8.5))
  args$pdf_height <- suppressWarnings(as.numeric(args$pdf_height %||% 6.5))
  if (!is.finite(args$pdf_width) || args$pdf_width <= 0) {
    stop("ERROR: pdf_width must be a positive number.", call. = FALSE)
  }
  if (!is.finite(args$pdf_height) || args$pdf_height <= 0) {
    stop("ERROR: pdf_height must be a positive number.", call. = FALSE)
  }

  top_n_raw <- suppressWarnings(as.integer(args$top_n_pathways %||% 20L))
  args$plot_all_pathways <- is.na(top_n_raw) || top_n_raw <= 0L
  args$top_n_pathways <- if (args$plot_all_pathways) 1L else top_n_raw
  args$top_n_by_sign <- gsea_vis_bool(args$top_n_by_sign, FALSE)
  args$pathway_bubble_plots <- gsea_vis_bool(args$pathway_bubble_plots, TRUE)
  args$pathway_bubble_significance_statistic <- gsea_cli_choice(
    non_empty(args$pathway_bubble_significance_statistic) %||% "padj",
    c("padj", "pval"),
    "padj",
    "pathway_bubble_significance_statistic"
  )
  args$collection_color_scale <- gsea_cli_choice(
    non_empty(args$collection_color_scale) %||% "independent",
    c("independent", "shared"),
    "independent",
    "collection_color_scale"
  )
  bubble_top_n_raw <- suppressWarnings(as.integer(args$pathway_bubble_top_n %||% 20L))
  if (is.na(bubble_top_n_raw) || bubble_top_n_raw < 0L) {
    stop("ERROR: pathway_bubble_top_n must be zero or a positive integer.", call. = FALSE)
  }
  args$pathway_bubble_top_n <- bubble_top_n_raw

  max_plots_raw <- suppressWarnings(as.integer(args$max_plots_in_pdf %||% 0L))
  args$max_plots_in_pdf <- if (is.na(max_plots_raw) || max_plots_raw <= 0L) 0L else max_plots_raw
  args$max_plots_unlimited <- args$max_plots_in_pdf == 0L
  max_le_genes_raw <- suppressWarnings(as.integer(args$max_le_genes_heatmap %||% 50L))
  args$max_le_genes_heatmap <- if (is.na(max_le_genes_raw) || max_le_genes_raw < 1L) 50L else max_le_genes_raw

  args$running_score_line_color_core <- switch(
    non_empty(args$running_score_line_color) %||% "ES sign",
    "red/blue by ES" = "ES sign",
    "red/blue by NES" = "ES sign",
    "ES sign" = "ES sign",
    "green" = "green",
    stop("ERROR: running_score_line_color must be 'ES sign', 'red/blue by ES', or 'green'.", call. = FALSE)
  )
  args$add_max_deviation_line_core <- switch(
    non_empty(args$add_max_deviation_line) %||% "both",
    "x-coordinate" = "coordinate",
    "y-coordinate" = "horizontal",
    "xy-coordinate" = "both",
    "coordinate" = "coordinate",
    "horizontal" = "horizontal",
    "both" = "both",
    "none" = FALSE,
    stop("ERROR: Invalid add_max_deviation_line value.", call. = FALSE)
  )

  args
}

#' Resolve exact selected rows once and reuse them for preflight and plotting
resolve_selection <- function(gsea_filter, config) {
  available_contrasts <- unique(as.character(gsea_filter$contrast))
  available_contrasts <- available_contrasts[!is.na(available_contrasts) & nzchar(available_contrasts)]

  plot_contrasts <- switch(
    config$contrast_filter,
    keep = config$contrast_names,
    remove = if (length(config$contrast_names) > 0) {
      setdiff(available_contrasts, config$contrast_names)
    } else {
      character(0)
    },
    none = character(0)
  )

  selected <- gsea_vis_select_rows(
    gsea_filter,
    plot_contrasts = plot_contrasts,
    plot_all_pathways = config$plot_all_pathways,
    top_n_pathways = config$top_n_pathways,
    top_n_by_sign = config$top_n_by_sign
  )
  if (!is.data.frame(selected) || nrow(selected) == 0) {
    stop("ERROR: No GSEA result rows matched the requested contrast and pathway selection.", call. = FALSE)
  }

  selected_before_cap <- nrow(selected)
  if (config$max_plots_in_pdf > 0L && nrow(selected) > config$max_plots_in_pdf) {
    selected <- utils::head(selected, config$max_plots_in_pdf)
  }

  list(
    plot_contrasts = plot_contrasts,
    selected_rows = selected,
    selected_before_cap = selected_before_cap,
    selected_after_cap = nrow(selected)
  )
}

#' Write standardized OMIX pathway-bubble figures alongside legacy panels
write_pathway_bubble_outputs <- function(
    gsea_results,
    output_dir,
    top_n_pathways,
    collection_color_scale,
    significance_statistic,
    contrasts = NULL) {
  if (!requireNamespace("OmixPathwayPlots", quietly = TRUE)) {
    stop(
      "ERROR: OmixPathwayPlots is required for shared pathway-bubble plots. ",
      "Use the current OMIX r-pathway runtime or disable shared bubble plots."
    )
  }
  plots <- OmixPathwayPlots::plot_pathway_bubble_set(
    gsea_results,
    input_format = "gsea",
    p_value_column = significance_statistic,
    top_n_pathways = top_n_pathways,
    selection_scopes = c(
      "combined_single_panel",
      "across_all_collections",
      "within_each_collection"
    ),
    selection_contrasts = contrasts,
    plot_contrasts = contrasts,
    collection_color_scale = collection_color_scale
  )
  OmixPathwayPlots::save_pathway_bubble_set(
    plots,
    output_dir = output_dir,
    file_prefix = "GSEA-Vis-Pathway-Bubble"
  )
}

#' Load the DEG Analysis expression/metadata bundle only when an LE panel remains requested.
resolve_heatmap_inputs <- function(config, deg_bundle) {
  plots_to_include <- config$plots_to_include
  if (!grepl("LE", plots_to_include, fixed = TRUE)) {
    return(list(plots_to_include = plots_to_include, batch_result = NULL))
  }

  # LE is an explicit user selection. Invalid paired inputs must fail rather
  # than silently changing the requested scientific output to ES+RNK.
  list(
    plots_to_include = plots_to_include,
    batch_result = construct_batch_result(
      deg_bundle$deg_table,
      deg_bundle$sample_metadata
    )
  )
}

#' Load all inputs in dependency order, minimizing unnecessary DEG processing
load_inputs <- function(config, data_root = "/data") {
  # Each hidden Code Ocean input has a fixed mount. Restricting discovery to
  # its own mount prevents unrelated workflow artifacts elsewhere under /data
  # from becoming accidental candidates.
  dataset_roots <- list(
    msigdb = file.path(data_root, "msigdb"),
    gsea_filter = file.path(data_root, "gsea_filter_results"),
    deg_bundle = file.path(data_root, "deg-training")
  )

  msigdb_path <- config$msigdb_database %||%
    find_msigdb_file(dataset_roots$msigdb)
  if (is.null(msigdb_path)) {
    stop(
      "ERROR: MSigDB database not found. Attach one MSigDB database data asset or provide --msigdb_database.",
      call. = FALSE
    )
  }
  if (!file.exists(msigdb_path)) {
    stop(
      sprintf("ERROR: MSigDB database not found at: %s", msigdb_path),
      call. = FALSE
    )
  }
  cat(sprintf("  MSigDB database: %s\n", msigdb_path))

  gsea_filter_path <- config$gsea_filter_results %||%
    find_gsea_filter_file(dataset_roots$gsea_filter)
  if (is.null(gsea_filter_path)) {
    stop(
      "ERROR: filtered_gsea_results.csv not found. Attach one GSEA Filter Result or provide --gsea_filter_results.",
      call. = FALSE
    )
  }
  gsea_filter <- read_gsea_filter(gsea_filter_path)
  cat(sprintf("  GSEA Filter rows: %d, columns: %d\n", nrow(gsea_filter), ncol(gsea_filter)))

  selection <- resolve_selection(gsea_filter, config)
  selection$selected_rows <- restore_msigdb_membership(selection$selected_rows, msigdb_path)
  requested_deg_contrasts <- unique(as.character(selection$selected_rows$contrast))

  explicit_deg_inputs <- resolve_explicit_deg_inputs(
    config$deg_table,
    config$sample_metadata
  )
  deg_bundle <- if (!is.null(explicit_deg_inputs)) {
    cat("  DEG Analysis input: explicit DEG table + sample metadata override\n")
    explicit_deg_inputs
  } else {
    find_deg_analysis_bundle(dataset_roots$deg_bundle)
  }
  cat(sprintf("  DEG Analysis bundle: %s\n", deg_bundle$directory))
  ranked_stats <- load_deg_table_as_ranked_stats(
    deg_bundle$deg_table,
    requested_deg_contrasts
  )

  heatmap <- resolve_heatmap_inputs(config, deg_bundle)
  list(
    gsea_filter = gsea_filter,
    selected_rows = selection$selected_rows,
    selection = selection,
    gsea_preranked = list(ranked_stats = ranked_stats),
    batch_result = heatmap$batch_result,
    plots_to_include = heatmap$plots_to_include
  )
}

#' Print exact preflight configuration and selected plot count
report_configuration <- function(config, inputs) {
  selection <- inputs$selection
  cat("\nConfiguration:\n")
  cat(sprintf("  Contrast filter: %s\n", config$contrast_filter))
  if (!identical(config$contrast_filter, "none") && length(config$contrast_names) > 0) {
    cat(sprintf("  Contrasts (%s): %s\n", config$contrast_filter, paste(config$contrast_names, collapse = ", ")))
  }
  cat(sprintf(
    "  Resolved contrasts: %s\n",
    paste(unique(as.character(inputs$selected_rows$contrast)), collapse = ", ")
  ))
  cat(sprintf("  Top N pathways: %s\n", if (config$plot_all_pathways) "all" else config$top_n_pathways))
  cat(sprintf("  Top N by ES sign: %s\n", if (config$top_n_by_sign) "yes" else "no"))
  cat(sprintf("  Max plots in PDF: %s\n", if (config$max_plots_unlimited) "no limit" else config$max_plots_in_pdf))
  cat(sprintf("  Plots selected: %d\n", selection$selected_after_cap))
  cat(sprintf("  Panels: %s\n", inputs$plots_to_include))

  if (grepl("LE", inputs$plots_to_include, fixed = TRUE)) {
    cat(sprintf("  Heatmap gene order: %s\n", config$heatmap_gene_order))
    if (identical(config$heatmap_gene_order, "cluster")) {
      cat(sprintf(
        "  Gene clustering: %s distance / %s linkage\n",
        config$heatmap_gene_clustering_distance,
        config$heatmap_gene_clustering_method
      ))
    }
    cat(sprintf("  Heatmap sample order: %s\n", config$heatmap_sample_order))
    if (identical(config$heatmap_sample_order, "cluster")) {
      cat(sprintf(
        "  Sample clustering: %s distance / %s linkage\n",
        config$heatmap_sample_clustering_distance,
        config$heatmap_sample_clustering_method
      ))
    }
    cat(sprintf("  Show heatmap gene names: %s\n", if (gsea_vis_bool(config$show_le_heatmap_gene_names, TRUE)) "yes" else "no"))
    cat(sprintf("  Show leading-edge rank numbers: %s\n", if (gsea_vis_bool(config$show_le_heatmap_rank_labels, TRUE)) "yes" else "no"))
  }
  cat(sprintf("  Output directory: %s\n\n", config$output_dir))

  if (selection$selected_before_cap > 50L) {
    cat(sprintf(
      "NOTE: %d plots matched the selection. Large plot sets may produce a large PDF.\n",
      selection$selected_before_cap
    ))
  }
  if (selection$selected_before_cap > selection$selected_after_cap) {
    cat(sprintf(
      "NOTE: Output is capped at %d plots; %d selected rows will not be plotted.\n\n",
      selection$selected_after_cap,
      selection$selected_before_cap - selection$selected_after_cap
    ))
  }
}

#' Execute the core visualization with already selected rows
run_visualization <- function(config, inputs) {
  cat("Generating plots...\n")
  result <- GSEA_Visualization_Local(
    gsea_filter_result = inputs$gsea_filter,
    gsea_preranked_result = inputs$gsea_preranked,
    selected_rows = inputs$selected_rows,
    batch_result = inputs$batch_result,
    max_plots_in_pdf = config$max_plots_in_pdf,
    stop_if_too_many_plots = FALSE,
    plots_to_include = inputs$plots_to_include,
    running_score_line_color = config$running_score_line_color_core,
    add_max_deviation_line = config$add_max_deviation_line_core,
    rank_area_color = config$rank_area_color,
    show_es_rank_bar = gsea_vis_bool(config$show_es_rank_bar, FALSE),
    show_rnk_peak_line = gsea_vis_bool(config$show_rnk_peak_line, TRUE),
    show_rnk_le_highlight = gsea_vis_bool(config$show_rnk_le_highlight, TRUE),
    show_es_le_highlight = gsea_vis_bool(config$show_es_le_highlight, TRUE),
    heatmap_gene_names_column = config$heatmap_gene_names_column,
    heatmap_sample_names_column = config$heatmap_sample_names_column,
    heatmap_group_column = config$heatmap_group_column,
    heatmap_transform = config$heatmap_transform,
    max_le_genes_heatmap = config$max_le_genes_heatmap,
    heatmap_gene_order = config$heatmap_gene_order,
    heatmap_sample_order = config$heatmap_sample_order,
    heatmap_gene_clustering_distance = config$heatmap_gene_clustering_distance,
    heatmap_gene_clustering_method = config$heatmap_gene_clustering_method,
    heatmap_sample_clustering_distance = config$heatmap_sample_clustering_distance,
    heatmap_sample_clustering_method = config$heatmap_sample_clustering_method,
    show_le_heatmap_gene_names = config$show_le_heatmap_gene_names,
    show_le_heatmap_sample_names = config$show_le_heatmap_sample_names,
    show_le_heatmap_rank_labels = config$show_le_heatmap_rank_labels,
    pdf_width = config$pdf_width,
    pdf_height = config$pdf_height,
    output_dir = config$output_dir
  )

  if (isTRUE(config$pathway_bubble_plots)) {
    bubble_contrasts <- if (length(inputs$selection$plot_contrasts) == 0L) {
      NULL
    } else {
      inputs$selection$plot_contrasts
    }
    bubble_outputs <- write_pathway_bubble_outputs(
      gsea_results = inputs$gsea_filter,
      output_dir = config$output_dir,
      top_n_pathways = config$pathway_bubble_top_n,
      collection_color_scale = config$collection_color_scale,
      significance_statistic = config$pathway_bubble_significance_statistic,
      contrasts = bubble_contrasts
    )
    result$files$pathway_bubble_manifest <- bubble_outputs$manifest
    cat(sprintf("Shared pathway-bubble manifest: %s\n", bubble_outputs$manifest))
  }

  result
}

#' Report generated files and concise summary
report_result <- function(result) {
  cat(sprintf("\n%s\n\n", result$message))
  if (is.data.frame(result$manifest) && nrow(result$manifest) > 0) {
    cat("Generated files:\n")
    if (!is.null(result$files$pdf) && file.exists(result$files$pdf)) {
      cat(sprintf("  PDF:                    %s\n", result$files$pdf))
    }
    if (!is.null(result$files$running_es) && file.exists(result$files$running_es)) {
      cat(sprintf("  Running ES (CSV):       %s\n", result$files$running_es))
    }
    if (!is.null(result$files$consistency) && file.exists(result$files$consistency)) {
      cat(sprintf("  Input consistency:      %s\n", result$files$consistency))
    }
    if (!is.null(result$files$skipped) && file.exists(result$files$skipped)) {
      cat(sprintf("  Skipped pathways:       %s\n", result$files$skipped))
    }

    cat("\nSummary:\n")
    cat(sprintf("  Total plots: %d\n", nrow(result$manifest)))
    cat(sprintf("  Contrasts: %d\n", length(unique(result$manifest$contrast))))
    cat(sprintf("  Pathways: %d\n", length(unique(result$manifest$pathway))))
  }

  if (is.list(result$consistency)) {
    cat("\nInput consistency checks:\n")
    cat(sprintf(
      "  Complete pathway matches: %d\n",
      result$consistency$all_genes_strict %||% 0L
    ))
    cat(sprintf(
      "  Partial pathway matches with complete leading edge: %d\n",
      result$consistency$leading_edge_strict %||% 0L
    ))
    cat(sprintf(
      "  Reconstructed ES boundary differences: %d\n",
      result$consistency$boundary_mismatch %||% 0L
    ))
    cat(sprintf(
      "  Skipped consistency failures: %d\n",
      result$consistency$failed %||% 0L
    ))
  }

  if (length(result$skipped) > 0) {
    cat(sprintf("\nSkipped %d plot(s):\n", length(result$skipped)))
    for (reason in utils::head(result$skipped, 5)) {
      cat(sprintf("  - %s\n", reason))
    }
    if (length(result$skipped) > 5) {
      cat(sprintf("  ... and %d more (see skipped file)\n", length(result$skipped) - 5))
    }
  }
}

#' Main execution
main <- function() {
  cat("=== GSEA Enrichment Plot Generator ===\n\n")
  config <- normalize_args(get_args())
  inputs <- load_inputs(config)
  report_configuration(config, inputs)
  result <- run_visualization(config, inputs)
  report_result(result)
  cat("\n=== GSEA Enrichment Plot generation complete ===\n")
  invisible(result)
}

if (!interactive() && !identical(Sys.getenv("OMIX_ADAPTER_TEST_MODE"), "1")) {
  result <- tryCatch(
    main(),
    error = function(error) {
      cat(sprintf("\nERROR: %s\n", conditionMessage(error)))
      quit(status = 1)
    }
  )
  quit(status = 0)
}
