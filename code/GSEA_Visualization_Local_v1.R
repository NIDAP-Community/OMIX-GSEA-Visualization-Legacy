if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
}

gsea_vis_bool <- function(value, default = FALSE) {
  if (is.null(value) || length(value) == 0) {
    return(isTRUE(default))
  }

  if (is.logical(value)) {
    return(isTRUE(value[[1]]))
  }

  tolower(as.character(value[[1]])) %in% c("true", "t", "yes", "y", "1")
}

gsea_vis_vector <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(character(0))
  }

  values <- as.character(value)
  if (length(values) == 1) {
    text <- trimws(values[[1]])
    text <- sub("^c\\s*\\(", "", text)
    text <- sub("\\)$", "", text)
    values <- unlist(strsplit(text, "[\r\n,]+"), use.names = FALSE)
  }

  values <- trimws(values)
  values <- gsub("^['\"]|['\"]$", "", values)
  values[!is.na(values) & nzchar(values)]
}

gsea_vis_safe_file_name <- function(value) {
  value <- tolower(as.character(value %||% "plot"))
  value <- gsub("[^a-z0-9]+", "_", value)
  value <- gsub("^_+|_+$", "", value)
  if (!nzchar(value)) {
    value <- "plot"
  }

  value
}

gsea_vis_gene_vector <- function(value) {
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])) {
    return(character(0))
  }

  values <- unlist(strsplit(as.character(value[[1]]), "[;,]", perl = TRUE), use.names = FALSE)
  values <- trimws(values)
  unique(values[!is.na(values) & nzchar(values)])
}

gsea_vis_gsea_table <- function(gsea_filter_result) {
  if (is.data.frame(gsea_filter_result)) {
    return(gsea_filter_result)
  }

  if (is.list(gsea_filter_result)) {
    for (field in c("filtered_gsea", "gsea", "results", "data")) {
      if (is.data.frame(gsea_filter_result[[field]])) {
        return(gsea_filter_result[[field]])
      }
    }
  }

  stop("GSEA-Vis expected GSEA-F result data.", call. = FALSE)
}

gsea_vis_ranked_stats <- function(gsea_preranked_result) {
  ranked_stats <- NULL
  if (is.list(gsea_preranked_result) && !is.data.frame(gsea_preranked_result)) {
    ranked_stats <- gsea_preranked_result$ranked_stats
    if (!is.list(ranked_stats) && is.data.frame(gsea_preranked_result$gsea)) {
      ranked_stats <- attr(gsea_preranked_result$gsea, "ranked_stats", exact = TRUE)
    }
    if (!is.list(ranked_stats) && is.data.frame(gsea_preranked_result$results)) {
      ranked_stats <- attr(gsea_preranked_result$results, "ranked_stats", exact = TRUE)
    }
  } else if (is.data.frame(gsea_preranked_result)) {
    ranked_stats <- attr(gsea_preranked_result, "ranked_stats", exact = TRUE)
  }

  if (!is.list(ranked_stats) || length(ranked_stats) == 0) {
    stop(
      "GSEA-Vis requires ranked statistics from GSEA-P. Rerun GSEA-P before running GSEA-Vis.",
      call. = FALSE
    )
  }

  ranked_stats
}

gsea_vis_gsdb_path <- function(gsdb_result) {
  path <- NULL
  if (is.list(gsdb_result) && !is.data.frame(gsdb_result)) {
    path <- gsdb_result$dataset_path %||% gsdb_result$path
  }
  if (is.null(path) && exists("gsdb_parquet_dir", mode = "function")) {
    path <- gsdb_parquet_dir()
  }

  path
}

gsea_vis_read_gsdb_memberships <- function(rows, gsdb_result) {
  empty <- setNames(vector("list", 0), character(0))
  if (!is.data.frame(rows) || nrow(rows) == 0) {
    return(empty)
  }
  path <- gsea_vis_gsdb_path(gsdb_result)
  if (is.null(path) || !dir.exists(path)) {
    return(empty)
  }
  if (!requireNamespace("arrow", quietly = TRUE) || !requireNamespace("dplyr", quietly = TRUE)) {
    return(empty)
  }
  if (!exists("gsdb_open_dataset", mode = "function")) {
    return(empty)
  }

  species <- unique(as.character(rows$pathways_database_species %||% character(0)))
  species <- species[!is.na(species) & nzchar(species)]
  collections <- unique(as.character(rows$collection))
  pathways <- unique(as.character(rows$pathway))
  dataset <- gsdb_open_dataset(path)
  filtered <- dataset |>
    dplyr::select(
      species,
      collection,
      gene_set_name,
      gene_symbol
    ) |>
    dplyr::filter(
      .data$collection %in% !!collections,
      .data$gene_set_name %in% !!pathways
    )
  if (length(species) > 0) {
    filtered <- filtered |>
      dplyr::filter(.data$species %in% !!species)
  }

  membership <- tryCatch(
    as.data.frame(dplyr::collect(filtered)),
    error = function(error) data.frame()
  )
  if (!is.data.frame(membership) || nrow(membership) == 0) {
    return(empty)
  }

  keys <- paste(membership$collection, membership$gene_set_name, sep = "\r")
  split(as.character(membership$gene_symbol), keys)
}

gsea_vis_row_genes <- function(row, stats, membership_map = list()) {
  candidate_fields <- c("inPathway_orthologs", "inPathway", "leadingEdge_orthologs", "leadingEdge")
  for (field in candidate_fields) {
    if (field %in% names(row)) {
      genes <- gsea_vis_gene_vector(row[[field]])
      genes <- intersect(genes, names(stats))
      if (length(genes) > 0) {
        return(genes)
      }
    }
  }

  key <- paste(as.character(row$collection), as.character(row$pathway), sep = "\r")
  genes <- unique(as.character(membership_map[[key]] %||% character(0)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  intersect(genes, names(stats))
}

gsea_vis_select_rows <- function(
  gsea_table,
  plot_contrasts = character(0),
  plot_all_pathways = FALSE,
  top_n_pathways = 1
) {
  plot_contrasts <- gsea_vis_vector(plot_contrasts)
  rows <- gsea_table

  if (length(plot_contrasts) > 0 && "contrast" %in% names(rows)) {
    rows <- rows[as.character(rows$contrast) %in% plot_contrasts, , drop = FALSE]
  }
  if (nrow(rows) == 0) {
    return(rows)
  }

  key_cols <- intersect(c("contrast", "collection", "pathway"), names(rows))
  rows <- rows[!duplicated(rows[key_cols]), , drop = FALSE]
  rows <- gsea_vis_rank_rows(rows)
  if (!isTRUE(plot_all_pathways)) {
    top_n_pathways <- max(1L, as.integer(top_n_pathways %||% 1L))
    by_cols <- intersect(c("contrast", "collection"), names(rows))
    if (length(by_cols) > 0) {
      groups <- split(
        seq_len(nrow(rows)),
        do.call(interaction, c(as.list(rows[, by_cols, drop = FALSE]), list(drop = TRUE)))
      )
      keep <- unlist(lapply(groups, function(index) {
        utils::head(index, top_n_pathways)
      }), use.names = FALSE)
      rows <- rows[sort(keep), , drop = FALSE]
    } else {
      rows <- utils::head(rows, top_n_pathways)
    }
  }

  rownames(rows) <- NULL
  rows
}

gsea_vis_preview_row_index <- function(rows, preview_contrast = "") {
  if (!is.data.frame(rows) || nrow(rows) == 0) {
    return(NA_integer_)
  }

  preview_contrast <- gsea_vis_vector(preview_contrast)
  preview_contrast <- if (length(preview_contrast) > 0) preview_contrast[[1]] else ""
  if (nzchar(preview_contrast) && "contrast" %in% names(rows)) {
    matched <- which(as.character(rows$contrast) == preview_contrast)
    if (length(matched) > 0) {
      return(matched[[1]])
    }
  }

  1L
}

gsea_vis_preview_row_indices <- function(rows, preview_contrasts = character(0), preview_contrast = "") {
  if (!is.data.frame(rows) || nrow(rows) == 0 || !"contrast" %in% names(rows)) {
    return(integer(0))
  }

  preview_contrasts <- gsea_vis_vector(preview_contrasts)
  if (length(preview_contrasts) == 0) {
    preview_contrasts <- gsea_vis_vector(preview_contrast)
  }
  if (length(preview_contrasts) == 0) {
    preview_contrasts <- unique(as.character(rows$contrast))
  }

  preview_contrasts <- unique(preview_contrasts)
  preview_contrasts <- preview_contrasts[preview_contrasts %in% as.character(rows$contrast)]
  indices <- vapply(preview_contrasts, function(contrast) {
    matched <- which(as.character(rows$contrast) == contrast)
    if (length(matched) == 0) NA_integer_ else matched[[1]]
  }, integer(1))

  indices[!is.na(indices)]
}

gsea_vis_rank_rows <- function(rows) {
  if (!is.data.frame(rows) || nrow(rows) == 0) {
    return(rows)
  }

  rank_col <- intersect(c("pval", "padj"), names(rows))
  rank_values <- if (length(rank_col) > 0) {
    suppressWarnings(as.numeric(rows[[rank_col[[1]]]]))
  } else {
    rep(Inf, nrow(rows))
  }
  rank_values[is.na(rank_values)] <- Inf
  contrast_values <- if ("contrast" %in% names(rows)) as.character(rows$contrast) else rep("", nrow(rows))
  collection_values <- if ("collection" %in% names(rows)) as.character(rows$collection) else rep("", nrow(rows))
  pathway_values <- if ("pathway" %in% names(rows)) as.character(rows$pathway) else rep("", nrow(rows))

  rows[order(contrast_values, rank_values, collection_values, pathway_values), , drop = FALSE]
}

gsea_vis_plot_title <- function(row) {
  paste(as.character(row$pathway), as.character(row$contrast), sep = " | ")
}

gsea_vis_plot_subtitle <- function(row) {
  nes <- suppressWarnings(as.numeric(row$NES))
  parts <- c(
    as.character(row$collection),
    if (!is.na(nes)) paste0("NES=", formatC(nes, format = "f", digits = 2)),
    paste0("p=", signif(suppressWarnings(as.numeric(row$pval)), 3)),
    paste0("adj p=", signif(suppressWarnings(as.numeric(row$padj)), 3))
  )
  paste(parts[!is.na(parts) & nzchar(parts)], collapse = "; ")
}

gsea_vis_format_es <- function(row, running_data) {
  es <- suppressWarnings(as.numeric(row$ES))
  if (is.na(es)) {
    max_row <- running_data[running_data$max_deviation, , drop = FALSE]
    es <- if (nrow(max_row) > 0) max_row$running_es[[1]] else NA_real_
  }
  if (is.na(es)) {
    return("GSEA Enrichment Score")
  }

  paste0("GSEA Enrichment Score = ", formatC(es, format = "f", digits = 2))
}

gsea_vis_set_line_color <- function(plot, line_color) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(plot)
  }

  for (index in seq_along(plot$layers)) {
    if (inherits(plot$layers[[index]]$geom, "GeomLine")) {
      plot$layers[[index]]$aes_params$colour <- line_color
      plot$layers[[index]]$aes_params$color <- line_color
    }
  }

  plot
}

gsea_vis_plot_color <- function(row, running_score_line_color = "ES sign") {
  running_score_line_color <- as.character(running_score_line_color %||% "ES sign")
  if (identical(running_score_line_color, "green")) {
    return("green2")
  }

  nes <- suppressWarnings(as.numeric(row$NES))
  margins <- gsea_vis_color_margins()
  if (!is.na(nes) && nes < 0) margins["dn"] else margins["up"]
}

gsea_vis_color_margins <- function(
  dn_hew = 260,
  up_hew = 0,
  md_hew = 94,
  dn_c = 90,
  up_c = 180,
  md_c = 0,
  dn_l = 30,
  up_l = 30,
  md_l = 97,
  a = 1
) {
  if (requireNamespace("colorspace", quietly = TRUE)) {
    dn <- colorspace::sequential_hcl(
      1,
      h = dn_hew,
      c. = c(dn_c),
      l = c(dn_l),
      fixup = TRUE,
      alpha = a
    )
    up <- colorspace::sequential_hcl(
      1,
      h = up_hew,
      c. = c(up_c),
      l = c(up_l),
      fixup = TRUE,
      alpha = a
    )
    md <- colorspace::sequential_hcl(
      1,
      h = md_hew,
      c. = c(md_c),
      l = c(md_l),
      fixup = TRUE,
      alpha = a
    )
    return(c(up = up, dn = dn, md = md))
  }

  c(up = "#b33430", dn = "#1f5aa6", md = "#f7f7f7")
}

gsea_vis_group_colors <- function(group_values) {
  group_levels <- unique(as.character(group_values))
  group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]
  n_groups <- length(group_levels)
  if (n_groups == 0) {
    return(c())
  }

  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    colors <- if (n_groups <= 8) {
      RColorBrewer::brewer.pal(8, "Accent")[seq_len(n_groups)]
    } else if (n_groups <= 12) {
      RColorBrewer::brewer.pal(12, "Set3")[seq_len(n_groups)]
    } else {
      grDevices::rainbow(n_groups)
    }
  } else {
    colors <- grDevices::rainbow(n_groups)
  }
  names(colors) <- group_levels
  colors
}

gsea_vis_running_es_data <- function(stats, genes, gsea_param = 1) {
  stats <- stats[!is.na(stats)]
  genes <- intersect(unique(as.character(genes)), names(stats))
  n_genes <- length(stats)
  hits <- names(stats) %in% genes
  n_hits <- sum(hits)
  n_misses <- n_genes - n_hits

  if (n_genes == 0 || n_hits == 0 || n_misses == 0) {
    stop(
      "GSEA-Vis could not calculate a running enrichment score for this pathway.",
      call. = FALSE
    )
  }

  weights <- abs(stats)^as.numeric(gsea_param %||% 1)
  hit_denom <- sum(weights[hits])
  increments <- numeric(n_genes)
  increments[hits] <- weights[hits] / hit_denom
  increments[!hits] <- -1 / n_misses
  running_es <- cumsum(increments)
  max_index <- which.max(abs(running_es))

  data.frame(
    index = seq_along(stats),
    gene = names(stats),
    score = as.numeric(stats),
    in_pathway = hits,
    running_es = running_es,
    max_deviation = seq_along(stats) == max_index,
    stringsAsFactors = FALSE
  )
}

gsea_vis_row_leading_edge_genes <- function(row, stats, running_data) {
  # Always derive leading edge from the current running ES data so that the
  # LE boundary aligns with max_deviation (the peak of the ES curve).
  # The stored leadingEdge column in the GSEA-F results was computed with the
  # original GSEA-P ranking which may differ from the DEG table ranking used
  # here, causing the LE highlight to diverge from the x-coordinate guide line.
  max_index <- which(running_data$max_deviation)[[1]]
  nes <- suppressWarnings(as.numeric(row$NES))
  if (!is.na(nes) && nes < 0) {
    return(running_data$gene[running_data$in_pathway & running_data$index >= max_index])
  }
  running_data$gene[running_data$in_pathway & running_data$index <= max_index]
}

gsea_vis_needs_le_heatmap <- function(plots_to_include) {
  grepl("LE", as.character(plots_to_include %||% ""), fixed = TRUE)
}

gsea_vis_expression_table <- function(batch_result) {
  if (is.null(batch_result)) {
    return(NULL)
  }
  if (is.data.frame(batch_result)) {
    return(batch_result)
  }
  if (!is.list(batch_result)) {
    return(NULL)
  }

  for (field in c("final_expression", "batch_corrected", "normalized_counts", "counts", "data")) {
    if (is.data.frame(batch_result[[field]])) {
      return(batch_result[[field]])
    }
  }

  NULL
}

gsea_vis_metadata_table <- function(batch_result) {
  if (!is.list(batch_result) || is.data.frame(batch_result)) {
    return(NULL)
  }

  for (field in c("metadata", "sample_metadata", "samples")) {
    if (is.data.frame(batch_result[[field]])) {
      return(batch_result[[field]])
    }
  }

  NULL
}

gsea_vis_resolve_column <- function(df, requested = "", preferred = character(0), type = "any") {
  if (!is.data.frame(df) || ncol(df) == 0) {
    return(NULL)
  }

  requested <- as.character(requested %||% "")
  requested <- requested[nzchar(requested)]
  if (length(requested) > 0 && requested[[1]] %in% names(df)) {
    return(requested[[1]])
  }

  preferred <- preferred[preferred %in% names(df)]
  if (length(preferred) > 0) {
    return(preferred[[1]])
  }

  if (identical(type, "character")) {
    choices <- names(df)[vapply(df, function(value) is.character(value) || is.factor(value), logical(1))]
    if (length(choices) > 0) {
      return(choices[[1]])
    }
  }

  if (identical(type, "numeric")) {
    choices <- names(df)[vapply(df, is.numeric, logical(1))]
    if (length(choices) > 0) {
      return(choices[[1]])
    }
  }

  names(df)[[1]]
}

gsea_vis_unique_gene_matrix <- function(mat, genes) {
  keep <- !is.na(genes) & nzchar(genes)
  mat <- mat[keep, , drop = FALSE]
  genes <- as.character(genes[keep])
  if (nrow(mat) == 0) {
    return(mat)
  }

  if (anyDuplicated(genes)) {
    split_rows <- split(seq_along(genes), genes)
    mat <- do.call(rbind, lapply(split_rows, function(indices) {
      colMeans(mat[indices, , drop = FALSE], na.rm = TRUE)
    }))
  } else {
    rownames(mat) <- genes
  }

  mat
}

gsea_vis_expression_data <- function(
  batch_result,
  gene_names_column = "gene",
  sample_names_column = "Sample",
  heatmap_group_column = "Group"
) {
  expression_df <- gsea_vis_expression_table(batch_result)
  if (!is.data.frame(expression_df) || nrow(expression_df) == 0) {
    return(NULL)
  }

  metadata_df <- gsea_vis_metadata_table(batch_result)
  gene_column <- gsea_vis_resolve_column(
    expression_df,
    requested = gene_names_column,
    preferred = c("gene", "Gene", "GeneName", "gene_name", "symbol", "Symbol"),
    type = "character"
  )
  if (is.null(gene_column)) {
    return(NULL)
  }

  sample_column <- gsea_vis_resolve_column(
    metadata_df,
    requested = sample_names_column,
    preferred = c("Sample", "sample", "sample_id", "SampleID"),
    type = "character"
  )
  group_column <- gsea_vis_resolve_column(
    metadata_df,
    requested = heatmap_group_column,
    preferred = c("Group", "group", "Condition", "condition", "Treatment", "treatment"),
    type = "character"
  )

  numeric_columns <- setdiff(names(expression_df)[vapply(expression_df, is.numeric, logical(1))], gene_column)
  if (length(numeric_columns) == 0) {
    return(NULL)
  }

  if (!is.null(metadata_df) && !is.null(sample_column)) {
    metadata_samples <- as.character(metadata_df[[sample_column]])
    numeric_columns <- intersect(metadata_samples, numeric_columns)
    metadata_df <- metadata_df[match(numeric_columns, metadata_samples), , drop = FALSE]
    rownames(metadata_df) <- numeric_columns
  } else {
    metadata_df <- NULL
  }

  if (length(numeric_columns) == 0) {
    return(NULL)
  }

  mat <- as.matrix(expression_df[, numeric_columns, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat <- gsea_vis_unique_gene_matrix(mat, expression_df[[gene_column]])

  list(
    matrix = mat,
    metadata = metadata_df,
    gene_column = gene_column,
    sample_column = sample_column,
    group_column = group_column
  )
}

gsea_vis_transform_heatmap_matrix <- function(mat, heatmap_transform = "z-score") {
  heatmap_transform <- as.character(heatmap_transform %||% "z-score")
  if (identical(heatmap_transform, "none")) {
    return(mat)
  }

  if (identical(heatmap_transform, "center by row median")) {
    medians <- apply(mat, 1, stats::median, na.rm = TRUE)
    return(mat - medians)
  }

  if (identical(heatmap_transform, "center by row mean")) {
    means <- rowMeans(mat, na.rm = TRUE)
    return(mat - means)
  }

  means <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  sds[is.na(sds) | sds == 0] <- 1
  sweep(sweep(mat, 1, means, "-"), 1, sds, "/")
}

gsea_vis_row_leading_edge_expression_genes <- function(row, expression_genes) {
  candidate_fields <- c("leadingEdge", "leadingEdge_orthologs")
  for (field in candidate_fields) {
    if (field %in% names(row)) {
      genes <- intersect(gsea_vis_gene_vector(row[[field]]), expression_genes)
      if (length(genes) > 0) {
        return(genes)
      }
    }
  }

  character(0)
}

gsea_vis_contrast_levels <- function(contrast, available_levels) {
  available_levels <- as.character(available_levels %||% character(0))
  available_levels <- available_levels[!is.na(available_levels) & nzchar(available_levels)]
  if (length(available_levels) == 0) {
    return(character(0))
  }

  contrast <- as.character(contrast %||% "")
  tokens <- unlist(strsplit(contrast, "[+\\-*/()]+", perl = TRUE), use.names = FALSE)
  tokens <- trimws(gsub("^['\"]|['\"]$", "", tokens))
  tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
  matched <- available_levels[tolower(available_levels) %in% tolower(tokens)]
  unique(matched)
}

gsea_vis_filter_expression_for_contrast <- function(expression_data, row) {
  if (
    is.null(expression_data) ||
      !is.matrix(expression_data$matrix) ||
      !is.data.frame(expression_data$metadata) ||
      is.null(expression_data$group_column) ||
      !expression_data$group_column %in% names(expression_data$metadata)
  ) {
    return(expression_data)
  }

  metadata <- expression_data$metadata
  group_values <- as.character(metadata[[expression_data$group_column]])
  names(group_values) <- rownames(metadata)
  group_values <- group_values[colnames(expression_data$matrix)]
  contrast_levels <- gsea_vis_contrast_levels(as.character(row$contrast), unique(group_values))
  if (length(contrast_levels) == 0) {
    return(expression_data)
  }

  keep_samples <- names(group_values)[group_values %in% contrast_levels]
  keep_samples <- intersect(colnames(expression_data$matrix), keep_samples)
  if (length(keep_samples) == 0) {
    return(expression_data)
  }

  expression_data$matrix <- expression_data$matrix[, keep_samples, drop = FALSE]
  expression_data$metadata <- metadata[keep_samples, , drop = FALSE]
  expression_data
}

gsea_vis_empty_plot <- function(message) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(message, call. = FALSE)
  }

  ggplot2::ggplot() +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = 0, label = message), size = 3.2, color = "#555555") +
    ggplot2::theme_void()
}

gsea_vis_le_rank_labels <- function(genes, running_data, row) {
  genes <- as.character(genes %||% character(0))
  if (
    length(genes) == 0 ||
      !is.data.frame(running_data) ||
      !"gene" %in% names(running_data) ||
      !"score" %in% names(running_data)
  ) {
    return(genes)
  }

  rank_genes <- genes
  native_leading_edge <- gsea_vis_gene_vector(row[["leadingEdge"]])
  ortholog_leading_edge <- gsea_vis_gene_vector(row[["leadingEdge_orthologs"]])
  if (length(native_leading_edge) == length(ortholog_leading_edge) && length(native_leading_edge) > 0) {
    native_to_ortholog <- ortholog_leading_edge
    names(native_to_ortholog) <- native_leading_edge
    ortholog_to_native <- native_leading_edge
    names(ortholog_to_native) <- ortholog_leading_edge

    native_matches <- native_to_ortholog[genes]
    ortholog_matches <- ortholog_to_native[genes]
    if (sum(!is.na(native_matches)) >= sum(!is.na(ortholog_matches))) {
      rank_genes[!is.na(native_matches)] <- native_matches[!is.na(native_matches)]
    } else {
      rank_genes[!is.na(ortholog_matches)] <- ortholog_matches[!is.na(ortholog_matches)]
    }
  }

  score_lookup <- running_data$score
  names(score_lookup) <- running_data$gene
  le_scores <- data.frame(
    display_gene = genes,
    rank_gene = rank_genes,
    score = unname(score_lookup[rank_genes]),
    stringsAsFactors = FALSE
  )
  le_scores <- le_scores[is.finite(le_scores$score), , drop = FALSE]
  le_scores <- le_scores[!duplicated(le_scores$display_gene), , drop = FALSE]
  if (nrow(le_scores) == 0) {
    return(genes)
  }

  nes <- suppressWarnings(as.numeric(row$NES))
  if (!is.na(nes) && nes < 0) {
    le_scores <- le_scores[order(le_scores$score, le_scores$display_gene, na.last = TRUE), , drop = FALSE]
  } else {
    le_scores <- le_scores[order(-le_scores$score, le_scores$display_gene, na.last = TRUE), , drop = FALSE]
  }

  labels <- paste0("LE#", seq_len(nrow(le_scores)), " ", le_scores$display_gene)
  names(labels) <- le_scores$display_gene
  matched <- match(genes, names(labels))
  ranked_labels <- genes
  ranked_labels[!is.na(matched)] <- labels[matched[!is.na(matched)]]
  ranked_labels
}

gsea_vis_le_heatmap_plot <- function(
  row,
  expression_data,
  running_data = NULL,
  heatmap_transform = "z-score",
  max_le_genes_heatmap = 50,
  cluster_le_heatmap_rows = TRUE,
  cluster_le_heatmap_columns = FALSE,
  show_le_heatmap_sample_names = FALSE,
  show_le_heatmap_gene_names = TRUE,
  show_le_heatmap_rank_labels = TRUE,
  order_le_heatmap_rows_by_rank = TRUE
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("GSEA-Vis leading-edge heatmaps require the ComplexHeatmap package.", call. = FALSE)
  }

  if (is.null(expression_data) || !is.matrix(expression_data$matrix)) {
    return(gsea_vis_empty_plot("No expression matrix is available for the leading-edge heatmap."))
  }
  expression_data <- gsea_vis_filter_expression_for_contrast(expression_data, row)

  genes <- gsea_vis_row_leading_edge_expression_genes(row, rownames(expression_data$matrix))
  max_le_genes_heatmap <- max(1L, as.integer(max_le_genes_heatmap %||% 50L))
  if (length(genes) == 0) {
    return(gsea_vis_empty_plot("No leading-edge genes matched the expression matrix."))
  }

  mat <- expression_data$matrix[genes, , drop = FALSE]
  mat <- gsea_vis_transform_heatmap_matrix(mat, heatmap_transform = heatmap_transform)
  row_display_labels <- rownames(mat)
  ranked_row_labels <- row_display_labels
  if (
    gsea_vis_bool(show_le_heatmap_rank_labels, TRUE) ||
      gsea_vis_bool(order_le_heatmap_rows_by_rank, TRUE)
  ) {
    ranked_row_labels <- gsea_vis_le_rank_labels(row_display_labels, running_data, row)
  }
  if (gsea_vis_bool(show_le_heatmap_rank_labels, TRUE)) {
    row_display_labels <- ranked_row_labels
  }

  if (gsea_vis_bool(order_le_heatmap_rows_by_rank, TRUE)) {
    row_ranks <- rep(NA_integer_, length(ranked_row_labels))
    ranked_rows <- grepl("^LE#[0-9]+\\b", ranked_row_labels)
    row_ranks[ranked_rows] <- suppressWarnings(as.integer(sub("^LE#([0-9]+)\\b.*$", "\\1", ranked_row_labels[ranked_rows])))
    if (any(!is.na(row_ranks))) {
      row_order <- order(is.na(row_ranks), row_ranks, rownames(mat), na.last = TRUE)
      mat <- mat[row_order, , drop = FALSE]
      row_display_labels <- row_display_labels[row_order]
      cluster_le_heatmap_rows <- FALSE
    }
  }
  if (nrow(mat) > max_le_genes_heatmap) {
    mat <- mat[seq_len(max_le_genes_heatmap), , drop = FALSE]
    row_display_labels <- row_display_labels[seq_len(max_le_genes_heatmap)]
  }

  group_values <- NULL
  top_annotation <- NULL
  metadata <- expression_data$metadata
  group_column <- expression_data$group_column
  if (is.data.frame(metadata) && !is.null(group_column) && group_column %in% names(metadata)) {
    group_values <- as.character(metadata[[group_column]])
    names(group_values) <- colnames(expression_data$matrix)
    group_values <- group_values[colnames(mat)]
    group_colors <- gsea_vis_group_colors(group_values)
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      Group = group_values,
      col = if (length(group_colors) > 0) list(Group = group_colors) else NULL,
      annotation_height = grid::unit(0.3, "cm"),
      annotation_name_side = "left",
      annotation_legend_param = list(
        title_gp = grid::gpar(fontsize = 6),
        grid_width = grid::unit(0.25, "cm"),
        labels_gp = grid::gpar(fontsize = 6)
      )
    )
  }

  if (!is.null(group_values) && length(group_values) == ncol(mat)) {
    sample_order <- order(group_values, colnames(mat), na.last = TRUE)
    mat <- mat[, sample_order, drop = FALSE]
    if (!is.null(top_annotation)) {
      group_values <- group_values[colnames(mat)]
      group_colors <- gsea_vis_group_colors(group_values)
      top_annotation <- ComplexHeatmap::HeatmapAnnotation(
        Group = group_values,
        col = if (length(group_colors) > 0) list(Group = group_colors) else NULL,
        annotation_height = grid::unit(0.3, "cm"),
        annotation_name_side = "left",
        annotation_legend_param = list(
          title_gp = grid::gpar(fontsize = 6),
          grid_width = grid::unit(0.25, "cm"),
          labels_gp = grid::gpar(fontsize = 6)
        )
      )
    }
  }

  heatmap_colors <- if (requireNamespace("circlize", quietly = TRUE)) {
    limit <- stats::quantile(abs(as.matrix(mat)), 0.95, na.rm = TRUE)
    if (!is.finite(limit) || limit == 0) {
      limit <- 1
    }
    margins <- gsea_vis_color_margins()
    circlize::colorRamp2(c(-limit, 0, limit), c(margins["dn"], margins["md"], margins["up"]), space = "LAB")
  } else {
    margins <- gsea_vis_color_margins()
    grDevices::colorRampPalette(c(margins["dn"], margins["md"], margins["up"]))(101)
  }
  heatmap_legend <- if (identical(as.character(heatmap_transform), "z-score")) "SD" else "log2"

  heatmap <- ComplexHeatmap::Heatmap(
    mat,
    name = heatmap_legend,
    col = heatmap_colors,
    top_annotation = top_annotation,
    cluster_rows = isTRUE(cluster_le_heatmap_rows),
    cluster_columns = isTRUE(cluster_le_heatmap_columns),
    show_column_names = gsea_vis_bool(show_le_heatmap_sample_names, FALSE),
    show_row_names = gsea_vis_bool(show_le_heatmap_gene_names, TRUE),
    row_labels = row_display_labels,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 7),
    heatmap_legend_param = list(
      color_bar = "continuous",
      title_gp = grid::gpar(fontsize = 6),
      labels_gp = grid::gpar(fontsize = 5)
    )
  )

  grob <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      heatmap,
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    )
  )
  patchwork::wrap_elements(full = grob)
}

gsea_vis_es_plot <- function(
  row,
  running_data,
  line_color,
  add_max_deviation_line = "coordinate",
  show_x_label = FALSE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("GSEA-Vis requires the ggplot2 package.", call. = FALSE)
  }

  add_max_deviation_line <- as.character(add_max_deviation_line %||% "coordinate")
  # Title and guide line annotations are always black so they remain neutral
  # regardless of curve color (green or red/blue by NES sign).
  annotation_color <- "black"
  max_row <- running_data[running_data$max_deviation, , drop = FALSE]
  es_range <- range(running_data$running_es, na.rm = TRUE)
  es_span <- diff(es_range)
  if (!is.finite(es_span) || es_span == 0) {
    es_span <- 1
  }
  score_clip <- stats::quantile(abs(running_data$score), 0.95, na.rm = TRUE)
  if (!is.finite(score_clip) || score_clip == 0) {
    score_clip <- max(abs(running_data$score), na.rm = TRUE)
  }
  if (!is.finite(score_clip) || score_clip == 0) {
    score_clip <- 1
  }
  score_limited <- pmax(pmin(running_data$score, score_clip), -score_clip)
  all_ranks <- running_data
  all_ranks$score_limited <- score_limited
  margins <- gsea_vis_color_margins()
  strip_y <- es_range[[1]] - (0.08 * es_span)
  plot <- ggplot2::ggplot(running_data, ggplot2::aes(x = index, y = running_es)) +
    ggplot2::geom_tile(
      data = all_ranks,
      ggplot2::aes(x = index, y = strip_y, fill = score_limited),
      inherit.aes = FALSE,
      width = 1,
      height = 0.06 * es_span,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      high = margins["up"],
      low = margins["dn"],
      mid = margins["md"],
      midpoint = 0,
      limits = c(-1, 1) * score_clip
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "#555555") +
    ggplot2::geom_line(linewidth = 0.75, color = line_color) +
    ggplot2::geom_rug(
      data = running_data[running_data$in_pathway, , drop = FALSE],
      ggplot2::aes(x = index),
      sides = "b",
      inherit.aes = FALSE,
      length = grid::unit(0.03, "npc"),
      linewidth = 0.15,
      color = "#555555"
    ) +
    ggplot2::labs(
      title = gsea_vis_format_es(row, running_data),
      x = if (isTRUE(show_x_label)) "Gene Rank" else NULL,
      y = "Running Score"
    ) +
    ggplot2::annotate(
      geom = "text",
      x = 1L,
      y = strip_y,
      label = "+",
      hjust = 1.5,
      vjust = 0.5,
      size = 2.5,
      fontface = "bold",
      color = margins["up"]
    ) +
    ggplot2::annotate(
      geom = "text",
      x = max(running_data$index),
      y = strip_y,
      label = "-",
      hjust = -0.5,
      vjust = 0.5,
      size = 2.5,
      fontface = "bold",
      color = margins["dn"]
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(color = annotation_color, hjust = 0.5, size = 11),
      plot.margin = grid::unit(c(l = 0, r = 0.03, t = 0, b = 0), "npc")
    )

  if (nrow(max_row) > 0 && identical(add_max_deviation_line, "coordinate")) {
    plot <- plot + ggplot2::geom_vline(
      data = max_row,
      ggplot2::aes(xintercept = index),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.3,
      color = annotation_color
    )
  } else if (nrow(max_row) > 0 && identical(add_max_deviation_line, "horizontal")) {
    plot <- plot + ggplot2::geom_hline(
      yintercept = max_row$running_es[[1]],
      linetype = "dashed",
      linewidth = 0.3,
      color = annotation_color
    )
  } else if (nrow(max_row) > 0 && identical(add_max_deviation_line, "both")) {
    plot <- plot +
      ggplot2::geom_vline(
        data = max_row,
        ggplot2::aes(xintercept = index),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 0.3,
        color = annotation_color
      ) +
      ggplot2::geom_hline(
        yintercept = max_row$running_es[[1]],
        linetype = "dashed",
        linewidth = 0.3,
        color = annotation_color
      )
  }

  y_limits <- c(strip_y - 0.04 * es_span, es_range[[2]] + 0.04 * es_span)
  plot <- plot + ggplot2::coord_cartesian(ylim = y_limits, clip = "off")

  if (nrow(running_data) >= 1000) {
    plot <- plot + ggplot2::scale_x_continuous(
      position = "bottom",
      limits = c(0, max(running_data$index) + 1),
      labels = function(value) paste0(value / 1000, "K")
    )
  } else {
    plot <- plot + ggplot2::scale_x_continuous(
      position = "bottom",
      limits = c(0, max(running_data$index) + 1)
    )
  }

  plot
}

gsea_vis_rank_plot <- function(
  row,
  running_data,
  leading_edge_genes,
  line_color,
  rank_area_color = "grey",
  display_leading_edge_genes = FALSE,
  number_of_leading_edge_genes_to_display = 10,
  font_size_of_leading_edge_genes = 3
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("GSEA-Vis requires the ggplot2 package.", call. = FALSE)
  }

  leading_edge_genes <- intersect(leading_edge_genes, running_data$gene)
  ranks <- running_data
  nes <- suppressWarnings(as.numeric(row$NES))
  if (!is.na(nes) && nes < 0) {
    ranks$leading_edge <- ranks$gene %in% leading_edge_genes & ranks$score < 0
  } else {
    ranks$leading_edge <- ranks$gene %in% leading_edge_genes & ranks$score > 0
  }
  ranks$LE <- ifelse(ranks$leading_edge, "LE", "Outside")
  ranks$LE <- factor(ranks$LE, levels = c("Outside", "LE"))
  ranks$display_score <- ranks$score
  score_range <- range(ranks$score, na.rm = TRUE)
  score_span <- diff(score_range)
  if (!is.finite(score_span) || score_span == 0) {
    score_span <- 1
  }
  color_margin <- line_color
  # Text annotations always use black (neutral) so the LE count label
  # is readable regardless of curve color mode.
  annotation_color <- "black"
  le_count <- sum(ranks$leading_edge, na.rm = TRUE)
  leading_edge_indices <- which(ranks$leading_edge)
  le_xmin <- if (length(leading_edge_indices) > 0) min(leading_edge_indices) else NA_real_
  le_xmax <- if (length(leading_edge_indices) > 0) max(leading_edge_indices) else NA_real_
  le_ticks <- ranks[ranks$leading_edge & is.finite(ranks$score), c("index", "score"), drop = FALSE]
  le_text_x <- nrow(ranks) / 2
  le_text_hjust <- 0.5
  le_text_y <- score_range[[1]] + 0.12 * score_span
  le_label <- paste(le_count, "genes in the Leading Edge")
  if (!is.na(nes) && nes < 0) {
    le_text_y <- score_range[[2]] - 0.12 * score_span
  }
  area_layers <- {
    # "red_blue": per-gene gradient fill (red→white→blue by score value),
    #   identical palette to the ES strip rug — fully independent of line_color.
    # "grey"    : smooth grey area (default).
    if (identical(as.character(rank_area_color %||% "grey"), "red_blue")) {
      margins_rnk  <- gsea_vis_color_margins()
      clip_rnk     <- stats::quantile(abs(ranks$score), 0.95, na.rm = TRUE)
      if (!is.finite(clip_rnk) || clip_rnk == 0) clip_rnk <- max(abs(ranks$score), na.rm = TRUE)
      if (!is.finite(clip_rnk) || clip_rnk == 0) clip_rnk <- 1
      ranks$fill_score <- pmax(pmin(ranks$score, clip_rnk), -clip_rnk)
      list(
        ggplot2::geom_col(
          ggplot2::aes(y = display_score, fill = fill_score),
          width = if (nrow(ranks) >= 1000) 10 else 1,
          color = NA, show.legend = FALSE, na.rm = TRUE
        ),
        ggplot2::scale_fill_gradient2(
          high     = unname(margins_rnk["up"]),
          low      = unname(margins_rnk["dn"]),
          mid      = unname(margins_rnk["md"]),
          midpoint = 0,
          limits   = c(-1, 1) * clip_rnk
        )
      )
    } else {
      list(
        ggplot2::geom_area(fill = "#AAAAAA", color = NA, alpha = 0.85, na.rm = TRUE)
      )
    }
  }
  plot <- ggplot2::ggplot(ranks, ggplot2::aes(x = index, y = display_score)) +
    area_layers +
    {
      if (is.finite(le_xmin) && is.finite(le_xmax)) {
        ggplot2::annotate(
          geom = "rect",
          xmin = le_xmin,
          xmax = le_xmax,
          ymin = -Inf,
          ymax = Inf,
          fill = color_margin,
          alpha = 0.08
        )
      }
    } +
    ggplot2::scale_y_continuous(limits = score_range) +
    ggplot2::labs(x = "Gene Rank", y = "Gene Score") +
    ggplot2::annotate(
      geom = "text",
      x = le_text_x,
      y = le_text_y,
      label = le_label,
      color = annotation_color,
      hjust = le_text_hjust,
      size = 2.9
    )

  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 9),
      plot.margin = grid::unit(c(l = 0, r = 0.03, t = 0, b = 0), "npc")
    )

  if (nrow(ranks) >= 1000) {
    plot <- plot + ggplot2::scale_x_continuous(
      position = "bottom",
      limits = c(0, nrow(ranks)),
      labels = function(value) paste0(value / 1000, "K")
    )
  } else {
    plot <- plot + ggplot2::scale_x_continuous(position = "bottom", limits = c(0, nrow(ranks)))
  }

  plot
}

gsea_vis_pathway_plot <- function(
  row,
  stats,
  genes,
  expression_data = NULL,
  plots_to_include = "ES+RNK",
  running_score_line_color = "ES sign",
  add_max_deviation_line = "coordinate",
  rank_area_color = "grey",
  display_leading_edge_genes = FALSE,
  number_of_leading_edge_genes_to_display = 10,
  font_size_of_leading_edge_genes = 3,
  heatmap_transform = "z-score",
  max_le_genes_heatmap = 50,
  cluster_le_heatmap_rows = TRUE,
  cluster_le_heatmap_columns = FALSE,
  show_le_heatmap_sample_names = FALSE,
  show_le_heatmap_gene_names = TRUE,
  show_le_heatmap_rank_labels = TRUE,
  order_le_heatmap_rows_by_rank = TRUE
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("GSEA-Vis requires the patchwork package for ES+RNK plots.", call. = FALSE)
  }

  plots_to_include <- as.character(plots_to_include %||% "ES+RNK")
  plots_to_include <- if (plots_to_include %in% c("ES", "ES+RNK", "ES+RNK+LE", "LE")) {
    plots_to_include
  } else {
    "ES+RNK+LE"
  }
  line_color <- gsea_vis_plot_color(row, running_score_line_color)
  running_data <- gsea_vis_running_es_data(stats, genes)
  leading_edge_genes <- gsea_vis_row_leading_edge_genes(row, stats, running_data)
  include_es <- grepl("ES", plots_to_include, fixed = TRUE)
  include_rnk <- grepl("RNK", plots_to_include, fixed = TRUE)
  include_le <- grepl("LE", plots_to_include, fixed = TRUE)

  plot_parts <- list()
  if (isTRUE(include_es)) {
    plot_parts[[length(plot_parts) + 1]] <- gsea_vis_es_plot(
      row,
      running_data,
      line_color = line_color,
      add_max_deviation_line = add_max_deviation_line,
      show_x_label = !include_rnk   # label appears only when ES is the bottom panel
    )
  }

  if (isTRUE(include_rnk)) {
    rank_plot <- gsea_vis_rank_plot(
      row,
      running_data,
      leading_edge_genes = leading_edge_genes,
      line_color = line_color,
      rank_area_color = rank_area_color,
      display_leading_edge_genes = display_leading_edge_genes,
      number_of_leading_edge_genes_to_display = number_of_leading_edge_genes_to_display,
      font_size_of_leading_edge_genes = font_size_of_leading_edge_genes
    )
    plot_parts[[length(plot_parts) + 1]] <- rank_plot
  }

  if (isTRUE(include_le)) {
    le_plot <- gsea_vis_le_heatmap_plot(
      row,
      expression_data = expression_data,
      running_data = running_data,
      heatmap_transform = heatmap_transform,
      max_le_genes_heatmap = max_le_genes_heatmap,
      cluster_le_heatmap_rows = cluster_le_heatmap_rows,
      cluster_le_heatmap_columns = cluster_le_heatmap_columns,
      show_le_heatmap_sample_names = show_le_heatmap_sample_names,
      show_le_heatmap_gene_names = show_le_heatmap_gene_names,
      show_le_heatmap_rank_labels = show_le_heatmap_rank_labels,
      order_le_heatmap_rows_by_rank = order_le_heatmap_rows_by_rank
    )
  }

  if (isTRUE(include_le) && length(plot_parts) > 0) {
    left_plot <- do.call(patchwork::wrap_plots, c(plot_parts, list(ncol = 1)))
    plot <- (left_plot | le_plot) + patchwork::plot_layout(widths = c(1.7, 0.9))
  } else if (isTRUE(include_le)) {
    plot <- le_plot
  } else {
    plot <- do.call(patchwork::wrap_plots, c(plot_parts, list(ncol = 1)))
  }

  plot + patchwork::plot_annotation(
    title = gsea_vis_plot_title(row),
    subtitle = gsea_vis_plot_subtitle(row)
  )
}

gsea_vis_running_es_output <- function(row, running_data) {
  data.frame(
    contrast = as.character(row$contrast),
    collection = as.character(row$collection),
    pathway = as.character(row$pathway),
    geneRank = running_data$index,
    runningES = running_data$running_es,
    gene = running_data$gene,
    inPathway = running_data$in_pathway,
    geneScore = running_data$score,
    stringsAsFactors = FALSE
  )
}

gsea_vis_manifest_row <- function(row, genes, plot_id, pdf_page, preview = FALSE, genes_in_heatmap = 0L) {
  data.frame(
    plot_id = plot_id,
    pdf_page = pdf_page,
    preview = isTRUE(preview),
    contrast = as.character(row$contrast),
    collection = as.character(row$collection),
    pathway = as.character(row$pathway),
    NES = suppressWarnings(as.numeric(row$NES)),
    ES = suppressWarnings(as.numeric(row$ES)),
    pval = suppressWarnings(as.numeric(row$pval)),
    padj = suppressWarnings(as.numeric(row$padj)),
    genes_in_plot = length(genes),
    genes_in_heatmap = as.integer(genes_in_heatmap),
    stringsAsFactors = FALSE
  )
}

GSEA_Visualization_Local <- function(
  gsea_filter_result,
  gsea_preranked_result,
  gsdb_result = NULL,
  batch_result = NULL,
  plot_contrasts = character(0),
  plot_all_pathways = FALSE,
  top_n_pathways = 1,
  preview_contrast = "",
  preview_contrasts = character(0),
  max_plots_in_pdf = 50,
  stop_if_too_many_plots = TRUE,
  plots_to_include = "ES+RNK+LE",
  running_score_line_color = "ES sign",
  add_max_deviation_line = "coordinate",
  rank_plot_mode = NULL,
  rank_area_color = "grey",
  display_leading_edge_genes = FALSE,
  number_of_leading_edge_genes_to_display = 10,
  font_size_of_leading_edge_genes = 3,
  heatmap_gene_names_column = "gene",
  heatmap_sample_names_column = "Sample",
  heatmap_group_column = "Group",
  heatmap_transform = "z-score",
  max_le_genes_heatmap = 50,
  cluster_le_heatmap_rows = TRUE,
  cluster_le_heatmap_columns = FALSE,
  show_le_heatmap_sample_names = FALSE,
  show_le_heatmap_gene_names = TRUE,
  show_le_heatmap_rank_labels = TRUE,
  order_le_heatmap_rows_by_rank = TRUE,
  pdf_width = 8.5,
  pdf_height = 6.5,
  image_width = 8,
  image_height = 5.5,
  image_dpi = 160,
  output_dir = file.path(getwd(), "_workflow_runtime", "outputs", "gsea_vis")
) {
  gsea_table <- gsea_vis_gsea_table(gsea_filter_result)
  if (!is.data.frame(gsea_table) || nrow(gsea_table) == 0) {
    return(list(
      manifest = data.frame(),
      plots = list(),
      message = "GSEA-Vis did not generate plots because GSEA-F returned no rows."
    ))
  }

  ranked_stats <- gsea_vis_ranked_stats(gsea_preranked_result)
  selected <- gsea_vis_select_rows(
    gsea_table,
    plot_contrasts = plot_contrasts,
    plot_all_pathways = isTRUE(plot_all_pathways),
    top_n_pathways = top_n_pathways
  )
  if (nrow(selected) == 0) {
    return(list(
      manifest = data.frame(),
      plots = list(),
      message = "GSEA-Vis did not generate plots because no rows matched the selected contrasts."
    ))
  }

  max_plots_in_pdf <- max(1L, as.integer(max_plots_in_pdf %||% 50L))
  if (nrow(selected) > max_plots_in_pdf && isTRUE(stop_if_too_many_plots)) {
    stop(
      "GSEA-Vis selected ",
      nrow(selected),
      " plots, which is above Max plots in PDF (",
      max_plots_in_pdf,
      "). Narrow the GSEA-F filter or increase the cap.",
      call. = FALSE
    )
  }
  if (nrow(selected) > max_plots_in_pdf) {
    selected <- utils::head(selected, max_plots_in_pdf)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(output_dir, pattern = "^GSEA-Vis-", full.names = TRUE), force = TRUE)

  output_stamp <- gsub("[^0-9]", "", format(Sys.time(), "%Y%m%d%H%M%OS3"))
  membership_map <- gsea_vis_read_gsdb_memberships(selected, gsdb_result)
  expression_data <- NULL
  if (gsea_vis_needs_le_heatmap(plots_to_include)) {
    expression_data <- gsea_vis_expression_data(
      batch_result,
      gene_names_column = heatmap_gene_names_column,
      sample_names_column = heatmap_sample_names_column,
      heatmap_group_column = heatmap_group_column
    )
    if (is.null(expression_data)) {
      stop(
        "GSEA-Vis needs Batch Correction final expression and sample metadata to draw leading-edge heatmaps.",
        call. = FALSE
      )
    }
  }
  plot_rows <- list()
  manifest_rows <- list()
  running_es_rows <- list()
  skipped <- list()
  pdf_file <- file.path(output_dir, paste0("GSEA-Vis-Enrichment-Plots-", output_stamp, ".pdf"))
  preview_file <- file.path(output_dir, paste0("GSEA-Vis-Preview-", output_stamp, ".png"))
  preview_files <- list()

  grDevices::pdf(pdf_file, width = as.numeric(pdf_width), height = as.numeric(pdf_height))
  on.exit(grDevices::dev.off(), add = TRUE)

  preview_row_indices <- gsea_vis_preview_row_indices(
    selected,
    preview_contrasts = preview_contrasts,
    preview_contrast = preview_contrast
  )
  preview_done <- logical(0)

  for (row_index in seq_len(nrow(selected))) {
    row <- selected[row_index, , drop = FALSE]
    contrast <- as.character(row$contrast)
    stats <- ranked_stats[[contrast]]
    if (!is.numeric(stats) || length(stats) == 0) {
      skipped[[length(skipped) + 1]] <- paste0(contrast, " / ", row$pathway, ": missing ranked stats")
      next
    }

    genes <- gsea_vis_row_genes(row, stats, membership_map = membership_map)
    if (length(genes) == 0) {
      skipped[[length(skipped) + 1]] <- paste0(contrast, " / ", row$pathway, ": no pathway genes matched ranked stats")
      next
    }
    running_data <- tryCatch(
      gsea_vis_running_es_data(stats, genes),
      error = function(error) {
        skipped[[length(skipped) + 1]] <<- paste0(contrast, " / ", row$pathway, ": ", conditionMessage(error))
        NULL
      }
    )
    if (is.null(running_data)) {
      next
    }

    plot_id <- paste(
      gsea_vis_safe_file_name(contrast),
      gsea_vis_safe_file_name(row$pathway),
      sep = "__"
    )
    plot <- gsea_vis_pathway_plot(
      row,
      stats = stats,
      genes = genes,
      expression_data = expression_data,
      plots_to_include = plots_to_include,
      running_score_line_color = running_score_line_color,
      add_max_deviation_line = add_max_deviation_line,
      rank_area_color = rank_area_color,
      display_leading_edge_genes = isTRUE(display_leading_edge_genes),
      number_of_leading_edge_genes_to_display = number_of_leading_edge_genes_to_display,
      font_size_of_leading_edge_genes = font_size_of_leading_edge_genes,
      heatmap_transform = heatmap_transform,
      max_le_genes_heatmap = max_le_genes_heatmap,
      cluster_le_heatmap_rows = isTRUE(cluster_le_heatmap_rows),
      cluster_le_heatmap_columns = isTRUE(cluster_le_heatmap_columns),
      show_le_heatmap_sample_names = gsea_vis_bool(show_le_heatmap_sample_names, FALSE),
      show_le_heatmap_gene_names = gsea_vis_bool(show_le_heatmap_gene_names, TRUE),
      show_le_heatmap_rank_labels = gsea_vis_bool(show_le_heatmap_rank_labels, TRUE),
      order_le_heatmap_rows_by_rank = gsea_vis_bool(order_le_heatmap_rows_by_rank, TRUE)
    )
    print(plot)

    preview_name <- names(preview_row_indices)[match(row_index, preview_row_indices)]
    should_preview <- length(preview_name) == 1 && !is.na(preview_name) && nzchar(preview_name)
    if (should_preview) {
      local_preview_file <- file.path(
        output_dir,
        paste0("GSEA-Vis-Preview-", gsea_vis_safe_file_name(preview_name), "-", output_stamp, ".png")
      )
      ggplot2::ggsave(
        local_preview_file,
        plot = plot,
        width = as.numeric(image_width),
        height = as.numeric(image_height),
        bg = "white"
      )
      preview_files[[preview_name]] <- local_preview_file
      if (length(preview_done) == 0) {
        file.copy(local_preview_file, preview_file, overwrite = TRUE)
      }
      preview_done <- c(preview_done, preview_name)
    }

    manifest_rows[[length(manifest_rows) + 1]] <- gsea_vis_manifest_row(
      row,
      genes = genes,
      plot_id = plot_id,
      pdf_page = length(manifest_rows) + 1L,
      preview = isTRUE(should_preview),
      genes_in_heatmap = if (!is.null(expression_data)) {
        length(gsea_vis_row_leading_edge_expression_genes(row, rownames(expression_data$matrix)))
      } else {
        0L
      }
    )
    running_es_rows[[length(running_es_rows) + 1]] <- gsea_vis_running_es_output(row, running_data)
    plot_rows[[length(plot_rows) + 1]] <- row
  }

  grDevices::dev.off()
  on.exit(NULL, add = FALSE)

  if (length(preview_done) == 0 && length(plot_rows) > 0) {
    first_row <- plot_rows[[1]]
    first_stats <- ranked_stats[[as.character(first_row$contrast)]]
    first_genes <- gsea_vis_row_genes(first_row, first_stats, membership_map = membership_map)
    first_plot <- gsea_vis_pathway_plot(
      first_row,
      stats = first_stats,
      genes = first_genes,
      expression_data = expression_data,
      plots_to_include = plots_to_include,
      running_score_line_color = running_score_line_color,
      add_max_deviation_line = add_max_deviation_line,
      rank_area_color = rank_area_color,
      display_leading_edge_genes = isTRUE(display_leading_edge_genes),
      number_of_leading_edge_genes_to_display = number_of_leading_edge_genes_to_display,
      font_size_of_leading_edge_genes = font_size_of_leading_edge_genes,
      heatmap_transform = heatmap_transform,
      max_le_genes_heatmap = max_le_genes_heatmap,
      cluster_le_heatmap_rows = isTRUE(cluster_le_heatmap_rows),
      cluster_le_heatmap_columns = isTRUE(cluster_le_heatmap_columns),
      show_le_heatmap_sample_names = gsea_vis_bool(show_le_heatmap_sample_names, FALSE),
      show_le_heatmap_gene_names = gsea_vis_bool(show_le_heatmap_gene_names, TRUE),
      show_le_heatmap_rank_labels = gsea_vis_bool(show_le_heatmap_rank_labels, TRUE),
      order_le_heatmap_rows_by_rank = gsea_vis_bool(order_le_heatmap_rows_by_rank, TRUE)
    )
    ggplot2::ggsave(
      preview_file,
      plot = first_plot,
      width = as.numeric(image_width),
      height = as.numeric(image_height),
      bg = "white"
    )
    preview_files[[as.character(first_row$contrast)]] <- preview_file
  }

  manifest <- if (length(manifest_rows) > 0) {
    do.call(rbind, manifest_rows)
  } else {
    data.frame(
      plot_id = character(0),
      pdf_page = integer(0),
      preview = logical(0),
      contrast = character(0),
      collection = character(0),
      pathway = character(0),
      NES = numeric(0),
      ES = numeric(0),
      pval = numeric(0),
      padj = numeric(0),
      genes_in_plot = integer(0),
      genes_in_heatmap = integer(0),
      stringsAsFactors = FALSE
    )
  }

  if (nrow(manifest) == 0 && file.exists(pdf_file)) {
    unlink(pdf_file, force = TRUE)
  }

  manifest_file <- file.path(output_dir, paste0("GSEA-Vis-Manifest-", output_stamp, ".csv"))
  utils::write.csv(manifest, manifest_file, row.names = FALSE)
  running_es <- if (length(running_es_rows) > 0) {
    do.call(rbind, running_es_rows)
  } else {
    data.frame(
      contrast = character(0),
      collection = character(0),
      pathway = character(0),
      geneRank = integer(0),
      runningES = numeric(0),
      gene = character(0),
      inPathway = logical(0),
      geneScore = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  running_es_file <- file.path(output_dir, paste0("GSEA-Vis-RunningES-", output_stamp, ".csv"))
  utils::write.csv(running_es, running_es_file, row.names = FALSE)
  skipped_file <- NULL
  if (length(skipped) > 0) {
    skipped_file <- file.path(output_dir, paste0("GSEA-Vis-Skipped-", output_stamp, ".csv"))
    utils::write.csv(
      data.frame(reason = unlist(skipped, use.names = FALSE), stringsAsFactors = FALSE),
      skipped_file,
      row.names = FALSE
    )
  }

  list(
    manifest = manifest,
    running_es = running_es,
    plots = list(
      preview = if (file.exists(preview_file)) normalizePath(preview_file, winslash = "/", mustWork = FALSE) else NULL,
      previews = lapply(preview_files, function(path) normalizePath(path, winslash = "/", mustWork = FALSE)),
      pdf = if (file.exists(pdf_file)) normalizePath(pdf_file, winslash = "/", mustWork = FALSE) else NULL
    ),
    files = list(
      preview = if (file.exists(preview_file)) normalizePath(preview_file, winslash = "/", mustWork = FALSE) else NULL,
      previews = lapply(preview_files, function(path) normalizePath(path, winslash = "/", mustWork = FALSE)),
      pdf = if (file.exists(pdf_file)) normalizePath(pdf_file, winslash = "/", mustWork = FALSE) else NULL,
      manifest = normalizePath(manifest_file, winslash = "/", mustWork = FALSE),
      running_es = normalizePath(running_es_file, winslash = "/", mustWork = FALSE),
      skipped = if (!is.null(skipped_file)) normalizePath(skipped_file, winslash = "/", mustWork = FALSE) else NULL
    ),
    skipped = unlist(skipped, use.names = FALSE),
    message = paste0(
      "GSEA-Vis complete: generated ",
      nrow(manifest),
      " ",
      as.character(plots_to_include %||% "ES+RNK"),
      " plot(s)",
      if (length(skipped) > 0) paste0("; skipped ", length(skipped), ".") else "."
    )
  )
}
