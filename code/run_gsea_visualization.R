#!/usr/bin/env Rscript

# GSEA Visualization CLI - reads GSEA Filter results and creates enrichment plots
# Usage: Rscript run_gsea_visualization.R [options]

suppressPackageStartupMessages({
  library(optparse)
})

# Source the core visualization function
source("/code/GSEA_Visualization_Local_v1.R")

#' Parse command-line arguments
get_args <- function() {
  option_list <- list(
    make_option(c("--gsea_filter_results"), type = "character", default = NULL,
                help = "Path to GSEA Filter results file (CSV or RDS flat table)"),
    make_option(c("--deg_table"), type = "character", default = NULL,
                help = "Path to DEG table (CSV or RDS with gene, ranking stat, and optionally contrast columns)"),
    make_option(c("--normalized_expression"), type = "character", default = NULL,
                help = "Path to normalized expression matrix (CSV or RDS, genes x samples) [required for LE heatmaps]"),
    make_option(c("--sample_metadata"), type = "character", default = NULL,
                help = "Path to sample metadata file (CSV or RDS with sample grouping) [required for LE heatmaps]"),
    make_option(c("--gsdb_result"), type = "character", default = NULL,
                help = "Path to GSDB result RDS file [optional]"),
    make_option(c("--contrast_filter"), type = "character", default = "none",
                help = "Contrast filter mode: none / keep / remove [default: none]"),
    make_option(c("--contrasts"), type = "character", default = "",
                help = "Comma-separated contrast names for keep/remove mode [default: '']"),
    make_option(c("--top_n_pathways"), type = "integer", default = 1L,
                help = "Top pathways per contrast/collection. 0 = all pathways. [default: 1]"),
    make_option(c("--max_plots_in_pdf"), type = "integer", default = 0L,
                help = "Total plots in PDF (global ceiling). 0 = no limit. [default: 0]"),
    make_option(c("--plots_to_include"), type = "character", default = "ES+RNK+LE",
                help = "Plot types: ES, ES+RNK, ES+RNK+LE, LE [default: ES+RNK+LE]"),
    make_option(c("--running_score_line_color"), type = "character", default = "ES sign",
                help = "Line color: 'ES sign' or 'green' [default: ES sign]"),
    make_option(c("--add_max_deviation_line"), type = "character", default = "coordinate",
                help = "Max deviation guide line: 'coordinate', 'full', or 'none' [default: coordinate]"),
    make_option(c("--rank_area_color"), type = "character", default = "grey",
                help = "RNK area fill color: 'grey' or 'red_blue' [default: grey]"),
    make_option(c("--heatmap_transform"), type = "character", default = "z-score",
                help = "Heatmap transform: 'z-score', 'center by row mean', 'none' [default: z-score]"),
    make_option(c("--max_le_genes_heatmap"), type = "integer", default = 50,
                help = "Max genes in LE heatmap [default: 50]"),
    make_option(c("--cluster_le_heatmap_rows"), type = "logical", default = TRUE,
                help = "Cluster heatmap rows [default: TRUE]"),
    make_option(c("--cluster_le_heatmap_columns"), type = "logical", default = FALSE,
                help = "Cluster heatmap columns [default: FALSE]"),
    make_option(c("--show_le_heatmap_gene_names"), type = "logical", default = TRUE,
                help = "Show gene names in heatmap [default: TRUE]"),
    make_option(c("--show_le_heatmap_sample_names"), type = "logical", default = FALSE,
                help = "Show sample names in heatmap [default: FALSE]"),
    make_option(c("--show_le_heatmap_rank_labels"), type = "logical", default = TRUE,
                help = "Show rank position labels in heatmap [default: TRUE]"),
    make_option(c("--order_le_heatmap_rows_by_rank"), type = "logical", default = TRUE,
                help = "Order heatmap rows by gene rank [default: TRUE]"),
    make_option(c("--heatmap_gene_names_column"), type = "character", default = "gene",
                help = "Gene column name in expression table [default: gene]"),
    make_option(c("--heatmap_sample_names_column"), type = "character", default = "Sample",
                help = "Sample column name in metadata [default: Sample]"),
    make_option(c("--heatmap_group_column"), type = "character", default = "Group",
                help = "Group column name in metadata [default: Group]"),
    make_option(c("--pdf_width"), type = "numeric", default = 8.5,
                help = "PDF width in inches [default: 8.5]"),
    make_option(c("--pdf_height"), type = "numeric", default = 6.5,
                help = "PDF height in inches [default: 6.5]"),
    make_option(c("--output_dir"), type = "character", default = "/results",
                help = "Output directory [default: /results]")
  )
  
  parser <- OptionParser(
    usage = "%prog [options]",
    description = "Generate GSEA enrichment plots from GSEA Filter results",
    option_list = option_list
  )
  
  optparse::parse_args(parser)
}

#' Load RDS file with error handling
load_rds <- function(path, description) {
  if (is.null(path) || !file.exists(path)) {
    stop(sprintf("ERROR: %s not found at: %s", description, path), call. = FALSE)
  }
  
  cat(sprintf("Loading %s from: %s\n", description, path))
  tryCatch(
    readRDS(path),
    error = function(e) {
      stop(sprintf("ERROR: Failed to load %s: %s", description, e$message), call. = FALSE)
    }
  )
}

#' Parse comma-separated string to vector
parse_comma_list <- function(value) {
  if (is.null(value) || nchar(trimws(value)) == 0) {
    return(character(0))
  }
  trimws(unlist(strsplit(value, ",", fixed = TRUE)))
}

#' Load file with automatic format detection
load_file <- function(path, description) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  
  cat(sprintf("Loading %s from: %s\n", description, path))
  ext <- tolower(tools::file_ext(path))
  
  tryCatch({
    if (ext == "rds") {
      readRDS(path)
    } else if (ext == "csv") {
      utils::read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      stop(sprintf("Unsupported file format: %s (expected .rds or .csv)", ext))
    }
  }, error = function(e) {
    stop(sprintf("Failed to load %s: %s", description, e$message), call. = FALSE)
  })
}

#' Construct batch result from separate expression and metadata files (for LE heatmaps)
construct_batch_result <- function(expression_path, metadata_path) {
  if (is.null(expression_path) && is.null(metadata_path)) {
    return(NULL)
  }

  expression_data <- load_file(expression_path, "normalized expression")
  sample_metadata <- load_file(metadata_path, "sample metadata")

  if (is.null(expression_data)) {
    stop("ERROR: Normalized expression file is required for LE heatmaps", call. = FALSE)
  }
  if (is.null(sample_metadata)) {
    warning("Sample metadata not provided; LE heatmaps will use default grouping")
  }

  # Ensure gene names are present as a column named 'gene'.
  # gsea_vis_expression_data() resolves a gene column by name — rownames alone are not found.
  gene_col_present <- any(c("gene", "Gene", "gene_name", "GeneName", "symbol", "Symbol") %in% names(expression_data))
  if (!gene_col_present) {
    rn <- rownames(expression_data)
    if (!is.null(rn) && length(rn) > 0 && nzchar(rn[[1]]) && !grepl("^[0-9]+$", rn[[1]])) {
      expression_data <- data.frame(gene = rn, expression_data, check.names = FALSE,
                                    stringsAsFactors = FALSE)
      cat("  Added 'gene' column from rownames for heatmap compatibility.\n")
    }
  }

  cat(sprintf("  Expression matrix: %d genes x %d samples\n",
              nrow(expression_data), ncol(expression_data) - 1L))
  if (!is.null(sample_metadata)) {
    cat(sprintf("  Sample metadata: %d samples x %d columns\n",
                nrow(sample_metadata), ncol(sample_metadata)))
  }

  list(
    final_expression  = expression_data,
    normalized_counts = expression_data,
    metadata          = sample_metadata
  )
}

#' Find the GSEA Filter flat table in the default data asset mount folder
find_gsea_filter_file <- function(mount_dir) {
  candidates <- c(
    file.path(mount_dir, "filtered_gsea_results.rds"),
    file.path(mount_dir, "filtered_gsea_results.csv")
  )
  for (p in candidates) {
    if (file.exists(p)) return(p)
  }
  NULL
}

#' Find the first RDS or CSV file in a data asset mount folder
find_first_file <- function(mount_dir) {
  if (!dir.exists(mount_dir)) return(NULL)
  candidates <- c(
    list.files(mount_dir, pattern = "\\.rds$", ignore.case = TRUE, full.names = TRUE),
    list.files(mount_dir, pattern = "\\.csv$", ignore.case = TRUE, full.names = TRUE)
  )
  if (length(candidates) > 0) candidates[[1]] else NULL
}

#' Load a DEG table and build a ranked_stats list for GSEA_Visualization_Local.
#'
#' Supports two layouts:
#'   WIDE  — one row per gene, per-contrast columns named <contrast>_tstat /
#'            <contrast>_logFC / <contrast>_FC (auto-detected by suffix).
#'   LONG  — one row per gene × contrast, with a 'contrast' column and a
#'            single ranking column (stat, log2FoldChange, logFC, …).
#'
#' Returns a named list: contrast -> sorted named numeric vector (gene -> stat).
load_deg_table_as_ranked_stats <- function(path) {
  if (is.null(path)) return(NULL)

  cat(sprintf("Loading DEG table from: %s\n", path))
  if (!file.exists(path)) {
    stop(sprintf("ERROR: DEG table not found at: %s", path), call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  deg <- tryCatch({
    if (ext == "rds") {
      readRDS(path)
    } else {
      utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }, error = function(e) {
    stop(sprintf("ERROR: Failed to load DEG table: %s", e$message), call. = FALSE)
  })

  if (!is.data.frame(deg)) {
    stop("ERROR: DEG table must be a data frame (CSV or RDS containing a data frame)", call. = FALSE)
  }

  # ── Detect gene column ──────────────────────────────────────────────────────
  gene_col <- intersect(
    c("Gene", "gene", "gene_name", "gene_id", "gene_symbol", "symbol",
      "SYMBOL", "GENE"),
    names(deg)
  )
  if (length(gene_col) == 0) {
    rn <- rownames(deg)
    if (length(rn) > 0 && nzchar(rn[[1]]) && !grepl("^[0-9]+$", rn[[1]])) {
      deg[["gene"]] <- rn
      gene_col <- "gene"
    } else {
      gene_col <- names(deg)[[1]]
    }
  } else {
    gene_col <- gene_col[[1]]
  }
  cat(sprintf("  Gene column     : %s\n", gene_col))

  # ── Detect layout ───────────────────────────────────────────────────────────
  # Wide: look for columns ending in _tstat, _logFC, _logFoldChange, _log2FC, _FC
  wide_suffixes <- c("_tstat", "_logFC", "_logFoldChange", "_log2FC", "_FC")
  wide_pattern  <- paste0("(", paste(gsub("_", "_", wide_suffixes), collapse="|"), ")$")
  wide_cols     <- grep(wide_pattern, names(deg), value = TRUE, ignore.case = TRUE)

  if (length(wide_cols) > 0) {
    # ── WIDE format ─────────────────────────────────────────────────────────
    # Prefer _tstat > _logFC / _log2FC / _logFoldChange > _FC
    pref_suffix <- c("_tstat", "_logFC", "_log2FC", "_logFoldChange", "_FC")
    chosen_suffix <- NULL
    chosen_cols   <- character(0)
    for (sfx in pref_suffix) {
      hits <- grep(paste0(sfx, "$"), names(deg), value = TRUE, ignore.case = TRUE)
      if (length(hits) > 0) {
        chosen_suffix <- sfx
        chosen_cols   <- hits
        break
      }
    }

    contrasts <- sub(paste0(chosen_suffix, "$"), "", chosen_cols, ignore.case = TRUE)
    cat(sprintf("  Layout          : wide (%d contrasts, ranking = *%s)\n",
                length(contrasts), chosen_suffix))
    cat(sprintf("  Contrasts       : %s\n", paste(contrasts, collapse = ", ")))

    genes <- as.character(deg[[gene_col]])
    ranked_stats <- setNames(
      lapply(seq_along(chosen_cols), function(i) {
        vals <- suppressWarnings(as.numeric(deg[[chosen_cols[i]]]))
        s <- setNames(vals, genes)
        s <- s[!is.na(s) & !duplicated(names(s)) & nzchar(names(s))]
        sort(s, decreasing = TRUE)
      }),
      contrasts
    )

  } else {
    # ── LONG format ─────────────────────────────────────────────────────────
    rank_col <- intersect(
      c("stat", "statistic", "log2FoldChange", "log2FC", "logFC", "LFC", "t", "score"),
      names(deg)
    )
    if (length(rank_col) == 0) {
      stop(
        "ERROR: DEG table has no recognised ranking column (stat, log2FoldChange, logFC) ",
        "and no wide-format contrast columns (<contrast>_tstat / _logFC / _FC).",
        call. = FALSE
      )
    }
    rank_col <- rank_col[[1]]

    contrast_col <- intersect(
      c("contrast", "condition", "comparison", "group", "Contrast"),
      names(deg)
    )
    cat(sprintf("  Layout          : long (ranking = %s)\n", rank_col))

    genes <- as.character(deg[[gene_col]])
    stats <- suppressWarnings(as.numeric(deg[[rank_col]]))
    contrasts <- if (length(contrast_col) > 0) {
      cat(sprintf("  Contrast column : %s\n", contrast_col[[1]]))
      as.character(deg[[contrast_col[[1]]]])
    } else {
      rep("default", nrow(deg))
    }

    valid <- !is.na(genes) & nzchar(genes) & !is.na(stats) & !is.na(contrasts)
    genes     <- genes[valid]
    stats     <- stats[valid]
    contrasts <- contrasts[valid]

    if (length(genes) == 0) {
      stop("ERROR: DEG table produced no valid gene-statistic pairs after filtering NA values.",
           call. = FALSE)
    }

    ranked_stats <- lapply(
      split(seq_along(genes), contrasts),
      function(idx) {
        s <- setNames(stats[idx], genes[idx])
        s <- s[!duplicated(names(s))]
        sort(s, decreasing = TRUE)
      }
    )
    cat(sprintf("  Contrasts       : %s\n", paste(names(ranked_stats), collapse = ", ")))
  }

  cat(sprintf("  Genes (1st)     : %d\n", length(ranked_stats[[1]])))
  ranked_stats
}

#' Treat empty or whitespace-only string as NULL (for type=file params)
non_empty <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) NULL else x
}

#' Main execution
main <- function() {
  args <- get_args()
  
  cat("=== GSEA Enrichment Plot Generator ===\n\n")
  
  # Coerce empty-string file params (sent by App Panel when no file uploaded) to NULL
  args$gsea_filter_results   <- non_empty(args$gsea_filter_results)
  args$deg_table             <- non_empty(args$deg_table)
  args$normalized_expression <- non_empty(args$normalized_expression)
  args$sample_metadata       <- non_empty(args$sample_metadata)
  # Coerce empty-string column-mapping params to sensible defaults
  args$heatmap_gene_names_column   <- non_empty(args$heatmap_gene_names_column)   %||% "gene"
  args$heatmap_sample_names_column <- non_empty(args$heatmap_sample_names_column) %||% "Sample"
  args$heatmap_group_column        <- non_empty(args$heatmap_group_column)        %||% "Group"
  args$add_max_deviation_line      <- non_empty(args$add_max_deviation_line)      %||% "coordinate"
  args$rank_area_color             <- non_empty(args$rank_area_color)             %||% "grey"

  # ── 1. Load GSEA Filter Results (flat filtered table) ──────────────────────
  gsea_filter_path <- args$gsea_filter_results %||%
    find_gsea_filter_file("/data/gsea_filter_results")

  if (is.null(gsea_filter_path) || !file.exists(gsea_filter_path)) {
    stop(
      sprintf(
        "ERROR: GSEA Filter results not found.\n  Searched: %s\nUpload a file or attach the GSEA Filter Results data asset.",
        gsea_filter_path %||% "/data/gsea_filter_results"
      ),
      call. = FALSE
    )
  }

  gsea_filter <- load_file(gsea_filter_path, "GSEA Filter results")
  if (!is.data.frame(gsea_filter)) {
    stop(
      "ERROR: GSEA Filter results must be a flat data frame (filtered_gsea_results.csv or .rds).",
      call. = FALSE
    )
  }
  cat(sprintf("  Rows: %d, Columns: %d\n", nrow(gsea_filter), ncol(gsea_filter)))

  # ── 2. Load DEG Table → ranked_stats ───────────────────────────────────────
  deg_table_path <- args$deg_table %||%
    find_first_file("/data/deg_table")

  if (is.null(deg_table_path)) {
    stop(
      "ERROR: A DEG table is required.\nUpload a DEG table file, or attach/replace the DEG Table data asset.",
      call. = FALSE
    )
  }
  ranked_stats <- load_deg_table_as_ranked_stats(deg_table_path)
  # Wrap in the structure expected by gsea_vis_ranked_stats()
  gsea_preranked <- list(ranked_stats = ranked_stats)

  # ── 3. Optional GSDB ───────────────────────────────────────────────────────
  gsdb_result <- NULL
  if (!is.null(args$gsdb_result) && file.exists(args$gsdb_result)) {
    gsdb_result <- load_rds(args$gsdb_result, "GSDB result")
  }

  # ── 4. Load expression data for LE heatmaps (optional) ────────────────────
  batch_result <- NULL
  norm_expr_path <- args$normalized_expression %||%
    find_first_file("/data/normalized_expression")
  meta_path <- args$sample_metadata %||%
    find_first_file("/data/sample_metadata")

  if (!is.null(norm_expr_path) || !is.null(meta_path)) {
    batch_result <- construct_batch_result(norm_expr_path, meta_path)
  }

  # ── 5. Resolve plots_to_include ────────────────────────────────────────────
  plots_to_include <- args$plots_to_include
  if (grepl("LE", plots_to_include, fixed = TRUE) && is.null(batch_result)) {
    warning("LE heatmaps requested but no expression data provided; switching to ES+RNK.")
    plots_to_include <- "ES+RNK"
  }

  # top_n_pathways: 0 = plot all; positive integer = top N per contrast × collection
  top_n_raw      <- args$top_n_pathways %||% 1L
  plot_all_paths <- is.na(top_n_raw) || top_n_raw <= 0L
  top_n_pathways <- if (plot_all_paths) 1L else as.integer(top_n_raw)

  # max_plots_in_pdf: 0 = no limit; positive integer = global PDF ceiling
  max_raw          <- args$max_plots_in_pdf %||% 0L
  max_plots_unlimited <- is.na(max_raw) || max_raw <= 0L
  max_plots_in_pdf <- if (max_plots_unlimited) 10000L else as.integer(max_raw)

  # Contrast filtering: none (all) / keep (include only) / remove (exclude)
  contrast_filter_mode <- non_empty(args$contrast_filter) %||% "none"
  contrast_names       <- parse_comma_list(args$contrasts)

  plot_contrasts <- switch(
    contrast_filter_mode,
    "keep"   = contrast_names,
    "remove" = {
      if (length(contrast_names) == 0) {
        character(0)  # nothing to remove → all contrasts
      } else {
        all_contrasts <- unique(as.character(gsea_filter[["contrast"]] %||% character(0)))
        setdiff(all_contrasts, contrast_names)
      }
    },
    character(0)  # "none" or unrecognised → all contrasts
  )

  preview_contrasts <- character(0)  # PNG previews disabled

  cat(sprintf("\nConfiguration:\n"))
  cat(sprintf("  Contrast filter: %s\n", contrast_filter_mode))
  if (contrast_filter_mode != "none" && length(contrast_names) > 0) {
    cat(sprintf("  Contrasts (%s): %s\n", contrast_filter_mode, paste(contrast_names, collapse = ", ")))
  }
  cat(sprintf("  Resolved contrasts: %s\n",
              if (length(plot_contrasts) > 0) paste(plot_contrasts, collapse = ", ") else "all"))
  cat(sprintf("  Top N pathways: %s\n", if (plot_all_paths) "all" else top_n_pathways))
  cat(sprintf("  Max plots in PDF: %s\n", if (max_plots_unlimited) "no limit" else max_plots_in_pdf))
  cat(sprintf("  Plots to include: %s\n", plots_to_include))
  cat(sprintf("  Output directory: %s\n\n", args$output_dir))
  
  # Estimate plot count and warn before starting (actual count may vary after filtering)
  n_rows_selected <- nrow(gsea_filter)
  if (!plot_all_paths) {
    n_contrasts  <- max(1L, length(unique(gsea_filter[["contrast"]])))
    n_collection <- max(1L, length(unique(gsea_filter[["collection"]])))
    n_rows_selected <- min(n_rows_selected, top_n_pathways * n_contrasts * n_collection)
  }
  if (n_rows_selected > 50) {
    cat(sprintf(
      "NOTE: Estimated %d plots to generate. Large plot sets may produce a big PDF and take\n",
      n_rows_selected
    ))
    cat("      several minutes. Set 'Max plots in PDF' to a number to cap the output.\n\n")
  }
  if (!max_plots_unlimited && n_rows_selected > max_plots_val) {
    cat(sprintf(
      "NOTE: Max plots in PDF is set to %d — output will be truncated from ~%d plots.\n\n",
      max_plots_val, n_rows_selected
    ))
  }

  # Run visualization
  cat("Generating plots...\n")
  result <- GSEA_Visualization_Local(
    gsea_filter_result    = gsea_filter,
    gsea_preranked_result = gsea_preranked,
    gsdb_result           = gsdb_result,
    batch_result          = batch_result,
    plot_contrasts        = plot_contrasts,
    plot_all_pathways     = plot_all_paths,
    top_n_pathways        = top_n_pathways,
    preview_contrast      = "",
    preview_contrasts     = preview_contrasts,
    max_plots_in_pdf      = max_plots_in_pdf,
    stop_if_too_many_plots = FALSE,        # never stop; always truncate to cap
    plots_to_include      = plots_to_include,
    running_score_line_color = switch(
      args$running_score_line_color %||% "red/blue by NES",
      "red/blue by ES"  = "ES sign",
      "red/blue by NES" = "ES sign",   # backward compat
      "green"           = "green",
      args$running_score_line_color   # passthrough for raw values
    ),
    add_max_deviation_line   = switch(
      args$add_max_deviation_line %||% "x-coordinate",
      "x-coordinate"  = "coordinate",
      "y-coordinate"  = "horizontal",
      "xy-coordinate" = "both",
      "both"          = "both",
      "none"          = FALSE,
      args$add_max_deviation_line
    ),
    rank_area_color          = args$rank_area_color %||% "grey",
    heatmap_gene_names_column   = args$heatmap_gene_names_column,
    heatmap_sample_names_column = args$heatmap_sample_names_column,
    heatmap_group_column        = args$heatmap_group_column,
    heatmap_transform           = args$heatmap_transform,
    max_le_genes_heatmap        = args$max_le_genes_heatmap,
    cluster_le_heatmap_rows     = args$cluster_le_heatmap_rows,
    cluster_le_heatmap_columns  = args$cluster_le_heatmap_columns,
    show_le_heatmap_gene_names  = args$show_le_heatmap_gene_names,
    show_le_heatmap_sample_names = args$show_le_heatmap_sample_names,
    show_le_heatmap_rank_labels  = args$show_le_heatmap_rank_labels,
    order_le_heatmap_rows_by_rank = args$order_le_heatmap_rows_by_rank,
    pdf_width  = args$pdf_width,
    pdf_height = args$pdf_height,
    output_dir = args$output_dir
  )
  
  # Report results
  cat(sprintf("\n%s\n\n", result$message))
  
  if (is.data.frame(result$manifest) && nrow(result$manifest) > 0) {
    cat("Generated files:\n")
    if (!is.null(result$files$pdf) && file.exists(result$files$pdf)) {
      cat(sprintf("  PDF: %s\n", result$files$pdf))
    }
    if (!is.null(result$files$preview) && file.exists(result$files$preview)) {
      cat(sprintf("  Preview: %s\n", result$files$preview))
    }
    if (!is.null(result$files$running_es) && file.exists(result$files$running_es)) {
      cat(sprintf("  Running ES: %s\n", result$files$running_es))
    }

    # Copy key outputs to stable filenames so App Panel tabs resolve correctly
    stable <- list(
      pdf        = file.path(args$output_dir, "GSEA-Vis-Enrichment-Plots.pdf"),
      running_es = file.path(args$output_dir, "GSEA-Vis-RunningES.csv")
    )
    if (!is.null(result$files$pdf) && file.exists(result$files$pdf)) {
      file.copy(result$files$pdf, stable$pdf, overwrite = TRUE)
      cat(sprintf("  Stable PDF: %s\n", stable$pdf))
    }
    if (!is.null(result$files$running_es) && file.exists(result$files$running_es)) {
      file.copy(result$files$running_es, stable$running_es, overwrite = TRUE)
      cat(sprintf("  Stable Running ES: %s\n", stable$running_es))
    }

    cat(sprintf("\nSummary:\n"))
    cat(sprintf("  Total plots: %d\n", nrow(result$manifest)))
    cat(sprintf("  Contrasts: %d\n", length(unique(result$manifest$contrast))))
    cat(sprintf("  Pathways: %d\n", length(unique(result$manifest$pathway))))
  }
  
  if (length(result$skipped) > 0) {
    cat(sprintf("\nSkipped %d plot(s):\n", length(result$skipped)))
    for (skip in head(result$skipped, 5)) {
      cat(sprintf("  - %s\n", skip))
    }
    if (length(result$skipped) > 5) {
      cat(sprintf("  ... and %d more (see skipped file)\n", length(result$skipped) - 5))
    }
  }
  
  cat("\n=== GSEA Enrichment Plot generation complete ===\n")
  invisible(result)
}

# Execute main if run as script
if (!interactive()) {
  result <- tryCatch(
    main(),
    error = function(e) {
      cat(sprintf("\nERROR: %s\n", e$message))
      quit(status = 1)
    }
  )
  quit(status = 0)
}
