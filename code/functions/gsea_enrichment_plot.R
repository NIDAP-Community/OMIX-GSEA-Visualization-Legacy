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

gsea_vis_numeric_scalar <- function(value, default = NA_real_) {
  values <- suppressWarnings(as.numeric(value))
  if (length(values) == 0) {
    return(as.numeric(default))
  }
  values[[1]]
}

gsea_vis_choice <- function(value, choices, default, parameter) {
  value <- as.character(value %||% default)
  value <- if (length(value) > 0) trimws(value[[1]]) else default
  if (!nzchar(value)) {
    value <- default
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
  if (is.null(value) || length(value) == 0) {
    return(character(0))
  }

  if (is.list(value) && length(value) == 1L) {
    value <- value[[1]]
  }
  values <- as.character(unlist(value, use.names = FALSE))
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(character(0))
  }

  values <- unlist(
    lapply(values, function(item) strsplit(item, "[;,]", perl = TRUE)[[1]]),
    use.names = FALSE
  )
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

gsea_vis_membership_candidate <- function(row, field, reference_genes) {
  genes <- if (field %in% names(row)) {
    gsea_vis_gene_vector(row[[field]])
  } else {
    character(0)
  }
  reference_genes <- as.character(reference_genes %||% character(0))
  mapped <- gsea_vis_map_genes_to_reference(
    genes,
    row = row,
    reference_genes = reference_genes
  )
  direct_index <- gsea_vis_match_gene_indices(genes, reference_genes)
  matched <- !is.na(mapped) & nzchar(mapped)

  list(
    field = field,
    genes = genes,
    mapped = mapped,
    matched_genes = unique(unname(mapped[matched])),
    missing_genes = unique(genes[!matched]),
    total = length(genes),
    matched_count = sum(matched),
    direct_count = sum(!is.na(direct_index)),
    overlap = if (length(genes) > 0) sum(matched) / length(genes) else NA_real_
  )
}

gsea_vis_choose_membership_candidate <- function(candidates, preferred_fields) {
  candidates <- candidates[vapply(candidates, function(candidate) candidate$total > 0L, logical(1))]
  if (length(candidates) == 0) {
    return(NULL)
  }

  fields <- vapply(candidates, `[[`, character(1), "field")
  matched <- vapply(candidates, `[[`, numeric(1), "matched_count")
  direct <- vapply(candidates, `[[`, numeric(1), "direct_count")
  overlap <- vapply(candidates, function(candidate) {
    if (is.finite(candidate$overlap)) candidate$overlap else -Inf
  }, numeric(1))
  priority <- match(fields, preferred_fields)
  priority[is.na(priority)] <- length(preferred_fields) + 1L
  candidates[[order(-matched, -direct, -overlap, priority)[[1]]]]
}

gsea_vis_authoritative_membership <- function(row, stats) {
  reference_genes <- names(stats)
  reference_genes <- as.character(reference_genes %||% character(0))

  pathway_fields <- c("inPathway_orthologs", "inPathway")
  pathway_candidates <- lapply(
    pathway_fields,
    function(field) gsea_vis_membership_candidate(row, field, reference_genes)
  )
  pathway <- gsea_vis_choose_membership_candidate(pathway_candidates, pathway_fields)

  result <- list(
    valid = FALSE,
    consistency_level = "failed_missing_pathway_membership",
    reason = "GSEA Filter results do not contain a non-empty inPathway membership list.",
    warning_message = NULL,
    pathway_field = NA_character_,
    leading_edge_field = NA_character_,
    pathway_genes_authoritative = character(0),
    pathway_genes_ranked = character(0),
    leading_edge_genes_authoritative = character(0),
    leading_edge_genes_ranked = character(0),
    pathway_total = 0L,
    pathway_matched = 0L,
    pathway_overlap = NA_real_,
    leading_edge_total = 0L,
    leading_edge_matched = 0L,
    missing_pathway_genes = character(0),
    missing_leading_edge_genes = character(0),
    leading_edge_not_in_pathway = character(0)
  )

  if (is.null(pathway)) {
    return(result)
  }

  result$pathway_field <- pathway$field
  result$pathway_genes_authoritative <- pathway$genes
  result$pathway_genes_ranked <- pathway$matched_genes
  result$pathway_total <- pathway$total
  result$pathway_matched <- pathway$matched_count
  result$pathway_overlap <- pathway$overlap
  result$missing_pathway_genes <- pathway$missing_genes

  preferred_leading_edge_fields <- if (identical(pathway$field, "inPathway_orthologs")) {
    c("leadingEdge_orthologs", "leadingEdge")
  } else {
    c("leadingEdge", "leadingEdge_orthologs")
  }
  leading_edge_candidates <- lapply(
    preferred_leading_edge_fields,
    function(field) gsea_vis_membership_candidate(row, field, reference_genes)
  )
  leading_edge <- gsea_vis_choose_membership_candidate(
    leading_edge_candidates,
    preferred_leading_edge_fields
  )

  if (is.null(leading_edge)) {
    result$consistency_level <- "failed_missing_leading_edge"
    result$reason <- "GSEA Filter results do not contain a non-empty leadingEdge list."
    return(result)
  }

  result$leading_edge_field <- leading_edge$field
  result$leading_edge_genes_authoritative <- leading_edge$genes
  result$leading_edge_genes_ranked <- leading_edge$matched_genes
  result$leading_edge_total <- leading_edge$total
  result$leading_edge_matched <- leading_edge$matched_count
  result$missing_leading_edge_genes <- leading_edge$missing_genes
  result$leading_edge_not_in_pathway <- setdiff(
    leading_edge$matched_genes,
    pathway$matched_genes
  )

  if (pathway$matched_count == 0L) {
    result$consistency_level <- "failed_no_pathway_gene_match"
    result$reason <- "No authoritative inPathway genes matched the ranked DEG input."
    return(result)
  }
  if (leading_edge$matched_count < leading_edge$total) {
    result$consistency_level <- "failed_incomplete_leading_edge"
    result$reason <- sprintf(
      "%d of %d authoritative leading-edge genes were absent from the ranked DEG input.",
      leading_edge$total - leading_edge$matched_count,
      leading_edge$total
    )
    return(result)
  }
  if (length(result$leading_edge_not_in_pathway) > 0) {
    result$consistency_level <- "failed_leading_edge_not_in_pathway"
    result$reason <- paste0(
      "Authoritative leading-edge genes were not present in the matched inPathway membership: ",
      paste(result$leading_edge_not_in_pathway, collapse = "; "),
      "."
    )
    return(result)
  }

  result$valid <- TRUE
  if (pathway$matched_count == pathway$total) {
    result$consistency_level <- "all_genes_strict"
    result$reason <- "Complete pathway and leading-edge match."
  } else {
    result$consistency_level <- "leading_edge_strict"
    result$reason <- sprintf(
      "%d of %d pathway genes matched; all %d leading-edge genes matched.",
      pathway$matched_count,
      pathway$total,
      leading_edge$total
    )
    result$warning_message <- result$reason
  }
  result
}

gsea_vis_row_genes <- function(row, stats) {
  membership <- gsea_vis_authoritative_membership(row, stats)
  membership$pathway_genes_ranked
}


gsea_vis_select_rows <- function(
  gsea_table,
  plot_contrasts = character(0),
  plot_all_pathways = FALSE,
  top_n_pathways = 1,
  top_n_by_sign = FALSE
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

    if (isTRUE(top_n_by_sign) && "ES" %in% names(rows)) {
      # Split each contrast × collection group by ES sign, then take top N from each sign.
      # Work on a plain data.frame copy to avoid mutating grouped tibbles.
      rows_df   <- as.data.frame(rows, stringsAsFactors = FALSE)
      es_vals   <- suppressWarnings(as.numeric(rows_df[["ES"]]))
      es_sign   <- ifelse(is.na(es_vals) | es_vals >= 0, "up", "down")
      group_factor <- if (length(by_cols) > 0) {
        interaction(
          c(as.list(rows_df[, by_cols, drop = FALSE]), list(es_sign)),
          drop = TRUE
        )
      } else {
        factor(es_sign)
      }
      groups <- split(seq_len(nrow(rows_df)), group_factor)
      keep   <- unlist(lapply(groups, function(idx) utils::head(idx, top_n_pathways)),
                       use.names = FALSE)
      rows <- rows[sort(keep), , drop = FALSE]
    } else if (length(by_cols) > 0) {
      groups <- split(
        seq_len(nrow(rows)),
        do.call(interaction, c(as.list(rows[, by_cols, drop = FALSE]), list(drop = TRUE)))
      )
      keep <- unlist(lapply(groups, function(idx) utils::head(idx, top_n_pathways)),
                     use.names = FALSE)
      rows <- rows[sort(keep), , drop = FALSE]
    } else {
      rows <- utils::head(rows, top_n_pathways)
    }
  }

  rownames(rows) <- NULL
  rows
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
  nes <- gsea_vis_numeric_scalar(row[["NES"]])
  parts <- c(
    as.character(row$collection),
    if (!is.na(nes)) paste0("NES=", formatC(nes, format = "f", digits = 2)),
    paste0("p=", signif(suppressWarnings(as.numeric(row$pval)), 3)),
    paste0("adj p=", signif(suppressWarnings(as.numeric(row$padj)), 3))
  )
  paste(parts[!is.na(parts) & nzchar(parts)], collapse = "; ")
}

gsea_vis_format_es <- function(row, running_data) {
  es <- gsea_vis_numeric_scalar(row[["ES"]])
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

  nes <- gsea_vis_numeric_scalar(row[["NES"]])
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

gsea_vis_reconstructed_leading_edge_genes <- function(running_data) {
  if (!is.data.frame(running_data) || nrow(running_data) == 0) {
    return(character(0))
  }
  peak_rows <- which(running_data$max_deviation)
  if (length(peak_rows) == 0) {
    return(character(0))
  }
  peak_index <- peak_rows[[1]]
  peak_es <- suppressWarnings(as.numeric(running_data$running_es[[peak_index]]))
  if (is.finite(peak_es) && peak_es < 0) {
    return(running_data$gene[running_data$in_pathway & running_data$index >= peak_index])
  }
  running_data$gene[running_data$in_pathway & running_data$index <= peak_index]
}

# Return the display-region boundary implied by the reconstructed running-ES
# curve. This is intentionally separate from authoritative leading-edge gene
# membership, which remains sourced from the upstream GSEA result.
gsea_vis_displayed_leading_edge_bounds <- function(running_data) {
  if (!is.data.frame(running_data) || nrow(running_data) == 0L) {
    return(c(xmin = NA_real_, xmax = NA_real_))
  }

  peak_rows <- which(running_data$max_deviation)
  if (length(peak_rows) == 0L) {
    return(c(xmin = NA_real_, xmax = NA_real_))
  }

  peak_index <- peak_rows[[1L]]
  peak_es <- suppressWarnings(as.numeric(running_data$running_es[[peak_index]]))
  if (!is.finite(peak_es)) {
    return(c(xmin = NA_real_, xmax = NA_real_))
  }

  if (peak_es < 0) {
    return(c(xmin = peak_index, xmax = max(running_data$index)))
  }

  c(xmin = min(running_data$index), xmax = peak_index)
}

gsea_vis_row_leading_edge_genes <- function(row, stats, running_data = NULL) {
  membership <- gsea_vis_authoritative_membership(row, stats)
  if (!isTRUE(membership$valid)) {
    stop(membership$reason, call. = FALSE)
  }
  membership$leading_edge_genes_ranked
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

gsea_vis_match_gene_indices <- function(query, reference) {
  query <- as.character(query %||% character(0))
  reference <- as.character(reference %||% character(0))
  if (length(query) == 0 || length(reference) == 0) {
    return(rep(NA_integer_, length(query)))
  }

  matched <- match(query, reference)
  missing <- is.na(matched)
  if (any(missing)) {
    matched[missing] <- match(tolower(query[missing]), tolower(reference))
  }
  matched
}

gsea_vis_identifier_pairs <- function(row) {
  pair_fields <- list(
    c("inPathway", "inPathway_orthologs"),
    c("leadingEdge", "leadingEdge_orthologs")
  )
  pairs <- list()

  for (fields in pair_fields) {
    left <- gsea_vis_gene_vector(row[[fields[[1]]]])
    right <- gsea_vis_gene_vector(row[[fields[[2]]]])
    if (length(left) > 0 && length(left) == length(right)) {
      pairs[[length(pairs) + 1L]] <- data.frame(
        left = left,
        right = right,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(pairs) == 0) {
    return(data.frame(left = character(0), right = character(0), stringsAsFactors = FALSE))
  }
  unique(do.call(rbind, pairs))
}

gsea_vis_map_genes_to_reference <- function(genes, row, reference_genes) {
  genes <- as.character(genes %||% character(0))
  reference_genes <- as.character(reference_genes %||% character(0))
  mapped <- rep(NA_character_, length(genes))
  names(mapped) <- genes
  if (length(genes) == 0 || length(reference_genes) == 0) {
    return(mapped)
  }

  direct_index <- gsea_vis_match_gene_indices(genes, reference_genes)
  direct <- !is.na(direct_index)
  mapped[direct] <- reference_genes[direct_index[direct]]

  pairs <- gsea_vis_identifier_pairs(row)
  if (nrow(pairs) == 0 || !any(is.na(mapped))) {
    return(mapped)
  }

  left_keys <- tolower(pairs$left)
  right_keys <- tolower(pairs$right)
  for (gene_index in which(is.na(mapped))) {
    key <- tolower(genes[[gene_index]])
    candidates <- unique(c(
      pairs$right[left_keys == key],
      pairs$left[right_keys == key]
    ))
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
    if (length(candidates) == 0) {
      next
    }
    candidate_index <- gsea_vis_match_gene_indices(candidates, reference_genes)
    candidate_index <- candidate_index[!is.na(candidate_index)]
    if (length(candidate_index) > 0) {
      mapped[[gene_index]] <- reference_genes[[candidate_index[[1]]]]
    }
  }

  mapped
}

gsea_vis_row_leading_edge_expression_genes <- function(
  row,
  leading_edge_genes,
  expression_genes
) {
  mapped <- gsea_vis_map_genes_to_reference(
    leading_edge_genes,
    row = row,
    reference_genes = expression_genes
  )
  unique(unname(mapped[!is.na(mapped) & nzchar(mapped)]))
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
  contrast_key <- as.character(row$contrast %||% "")
  contrast_key <- if (length(contrast_key) > 0) contrast_key[[1]] else ""

  if (
    is.list(expression_data) &&
      identical(expression_data$.gsea_vis_contrast, contrast_key)
  ) {
    return(expression_data)
  }

  if (
    is.null(expression_data) ||
      !is.matrix(expression_data$matrix) ||
      !is.data.frame(expression_data$metadata) ||
      is.null(expression_data$group_column) ||
      !expression_data$group_column %in% names(expression_data$metadata)
  ) {
    if (is.list(expression_data)) {
      expression_data$.gsea_vis_contrast <- contrast_key
    }
    return(expression_data)
  }

  metadata <- expression_data$metadata
  group_values <- as.character(metadata[[expression_data$group_column]])
  names(group_values) <- rownames(metadata)
  group_values <- group_values[colnames(expression_data$matrix)]
  contrast_levels <- gsea_vis_contrast_levels(contrast_key, unique(group_values))
  if (length(contrast_levels) > 0) {
    keep_samples <- names(group_values)[group_values %in% contrast_levels]
    keep_samples <- intersect(colnames(expression_data$matrix), keep_samples)
    if (length(keep_samples) > 0) {
      expression_data$matrix <- expression_data$matrix[, keep_samples, drop = FALSE]
      expression_data$metadata <- metadata[keep_samples, , drop = FALSE]
    }
  }

  expression_data$.gsea_vis_contrast <- contrast_key
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
  if (length(genes) == 0) {
    return(setNames(character(0), character(0)))
  }

  # Always return a named label vector keyed by the heatmap matrix row names.
  # ComplexHeatmap recommends keyed labels so correspondence is retained when
  # rows are reordered or clustered.
  labels <- setNames(genes, genes)
  if (
    !is.data.frame(running_data) ||
      !"gene" %in% names(running_data) ||
      !"score" %in% names(running_data)
  ) {
    labels[] <- paste0(seq_along(genes), ". ", genes)
    return(labels)
  }

  running_genes <- as.character(running_data$gene)
  running_scores <- suppressWarnings(as.numeric(running_data$score))
  mapped_running_genes <- gsea_vis_map_genes_to_reference(
    genes,
    row = row,
    reference_genes = running_genes
  )
  running_index <- gsea_vis_match_gene_indices(
    unname(mapped_running_genes),
    running_genes
  )

  matched_scores <- rep(NA_real_, length(genes))
  valid_matches <- !is.na(running_index)
  matched_scores[valid_matches] <- running_scores[running_index[valid_matches]]

  # Use score-based LE order whenever identifiers can be resolved. Direction
  # follows the current running-ES peak; unresolved genes follow in matrix order.
  peak_rows <- which(running_data$max_deviation)
  peak_es <- if (length(peak_rows) > 0) {
    suppressWarnings(as.numeric(running_data$running_es[[peak_rows[[1]]]]))
  } else {
    NA_real_
  }
  score_order <- which(is.finite(matched_scores))
  if (length(score_order) > 0) {
    if (is.finite(peak_es) && peak_es < 0) {
      score_order <- score_order[order(matched_scores[score_order], genes[score_order])]
    } else {
      score_order <- score_order[order(-matched_scores[score_order], genes[score_order])]
    }
  }

  unresolved <- setdiff(seq_along(genes), score_order)
  if (length(unresolved) > 0) {
    unresolved <- unresolved[order(unresolved, genes[unresolved], na.last = TRUE)]
  }

  rank_order <- c(score_order, unresolved)
  rank_number <- integer(length(genes))
  rank_number[rank_order] <- seq_along(rank_order)
  labels[] <- paste0(rank_number, ". ", genes)
  labels
}

gsea_vis_le_heatmap_plot <- function(
  row,
  expression_data,
  running_data = NULL,
  leading_edge_genes = character(0),
  heatmap_transform = "z-score",
  max_le_genes_heatmap = 50,
  heatmap_gene_order = "rank",
  heatmap_sample_order = "group",
  heatmap_gene_clustering_distance = "euclidean",
  heatmap_gene_clustering_method = "complete",
  heatmap_sample_clustering_distance = "euclidean",
  heatmap_sample_clustering_method = "complete",
  show_le_heatmap_sample_names = FALSE,
  show_le_heatmap_gene_names = TRUE,
  show_le_heatmap_rank_labels = TRUE,
  pdf_height = 6.5,
  pdf_width  = 8.5
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("GSEA-Vis leading-edge heatmaps require the ComplexHeatmap package.", call. = FALSE)
  }

  gene_order_mode <- gsea_vis_choice(
    heatmap_gene_order,
    choices = c("rank", "cluster", "input"),
    default = "rank",
    parameter = "heatmap_gene_order"
  )
  sample_order_mode <- gsea_vis_choice(
    heatmap_sample_order,
    choices = c("group", "cluster", "input"),
    default = "group",
    parameter = "heatmap_sample_order"
  )
  distance_choices <- c(
    "euclidean", "maximum", "manhattan", "canberra", "binary",
    "minkowski", "pearson", "spearman", "kendall"
  )
  method_choices <- c(
    "ward.D", "ward.D2", "single", "complete", "average",
    "mcquitty", "median", "centroid"
  )
  gene_clustering_distance <- gsea_vis_choice(
    heatmap_gene_clustering_distance,
    choices = distance_choices,
    default = "euclidean",
    parameter = "heatmap_gene_clustering_distance"
  )
  gene_clustering_method <- gsea_vis_choice(
    heatmap_gene_clustering_method,
    choices = method_choices,
    default = "complete",
    parameter = "heatmap_gene_clustering_method"
  )
  sample_clustering_distance <- gsea_vis_choice(
    heatmap_sample_clustering_distance,
    choices = distance_choices,
    default = "euclidean",
    parameter = "heatmap_sample_clustering_distance"
  )
  sample_clustering_method <- gsea_vis_choice(
    heatmap_sample_clustering_method,
    choices = method_choices,
    default = "complete",
    parameter = "heatmap_sample_clustering_method"
  )

  if (is.null(expression_data) || !is.matrix(expression_data$matrix)) {
    return(gsea_vis_empty_plot("No expression matrix is available for the leading-edge heatmap."))
  }
  expression_data <- gsea_vis_filter_expression_for_contrast(expression_data, row)

  genes <- gsea_vis_row_leading_edge_expression_genes(
    row,
    leading_edge_genes = leading_edge_genes,
    expression_genes = rownames(expression_data$matrix)
  )
  max_le_genes_heatmap <- max(1L, as.integer(max_le_genes_heatmap %||% 50L))
  if (length(genes) == 0) {
    return(gsea_vis_empty_plot("No leading-edge genes matched the expression matrix."))
  }

  mat <- expression_data$matrix[genes, , drop = FALSE]
  mat <- gsea_vis_transform_heatmap_matrix(mat, heatmap_transform = heatmap_transform)

  # Leading-edge ranks are independent of display order. Keep them keyed by matrix row
  # name so the correct label follows each gene through ordering and clustering.
  ranked_label_map <- gsea_vis_le_rank_labels(rownames(mat), running_data, row)
  ranked_row_labels <- unname(ranked_label_map[rownames(mat)])
  row_ranks <- suppressWarnings(
    as.integer(sub("^([0-9]+)\\.[[:space:]].*$", "\\1", ranked_row_labels))
  )

  # Rank mode uses the strongest LE genes first. Cluster and input modes retain
  # the incoming leading-edge subset before applying their requested ordering.
  if (identical(gene_order_mode, "rank") && any(!is.na(row_ranks))) {
    row_order <- order(is.na(row_ranks), row_ranks, rownames(mat), na.last = TRUE)
    mat <- mat[row_order, , drop = FALSE]
  }
  if (nrow(mat) > max_le_genes_heatmap) {
    mat <- mat[seq_len(max_le_genes_heatmap), , drop = FALSE]
  }

  show_rank_labels <- gsea_vis_bool(show_le_heatmap_rank_labels, TRUE)
  show_gene_names <- gsea_vis_bool(show_le_heatmap_gene_names, TRUE)
  row_display_labels <- if (show_rank_labels) {
    ranked_label_map[rownames(mat)]
  } else {
    setNames(rownames(mat), rownames(mat))
  }
  # A requested numbered rank label is itself a row label, so it must not be hidden
  # by a conflicting gene-name visibility setting.
  show_row_names <- show_gene_names || show_rank_labels
  cluster_rows <- identical(gene_order_mode, "cluster") && nrow(mat) > 1L

  group_values <- NULL
  top_annotation <- NULL
  metadata <- expression_data$metadata
  group_column <- expression_data$group_column
  if (is.data.frame(metadata) && !is.null(group_column) && group_column %in% names(metadata)) {
    group_values <- as.character(metadata[[group_column]])
    names(group_values) <- rownames(metadata)
    group_values <- group_values[colnames(mat)]
  }

  # Sample ordering is mutually exclusive: metadata group order, hierarchical
  # clustering, or the existing expression-matrix order.
  if (
    identical(sample_order_mode, "group") &&
      !is.null(group_values) &&
      length(group_values) == ncol(mat)
  ) {
    sample_order <- order(group_values, colnames(mat), na.last = TRUE)
    mat <- mat[, sample_order, drop = FALSE]
    group_values <- group_values[colnames(mat)]
  }
  cluster_columns <- identical(sample_order_mode, "cluster") && ncol(mat) > 1L

  if (!is.null(group_values) && length(group_values) == ncol(mat)) {
    group_values <- group_values[colnames(mat)]
    group_colors <- gsea_vis_group_colors(group_values)
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      Group = group_values,
      col = if (length(group_colors) > 0) list(Group = group_colors) else NULL,
      annotation_height = grid::unit(0.3, "cm"),
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        title_gp = grid::gpar(fontsize = 6),
        grid_width = grid::unit(0.25, "cm"),
        labels_gp = grid::gpar(fontsize = 6)
      )
    )
  }

  # ── Bounded-square cell sizing + auto font ────────────────────────────────
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  avail_h_pt <- max(36, (as.numeric(pdf_height) - 0.8) * 72)
  avail_w_body_pt <- max(20, as.numeric(pdf_width) * 0.346 * 72 - 128)
  cell_w_driven  <- avail_w_body_pt / n_cols
  cell_h_natural <- avail_h_pt / n_rows
  cell_size_pt   <- max(4, min(cell_w_driven, cell_h_natural, 14))
  row_font_size <- max(4, min(9, cell_size_pt * 0.65))
  col_font_size <- max(4, min(7, cell_w_driven * 0.50))

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
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    clustering_distance_rows = gene_clustering_distance,
    clustering_method_rows = gene_clustering_method,
    clustering_distance_columns = sample_clustering_distance,
    clustering_method_columns = sample_clustering_method,
    show_column_names = gsea_vis_bool(show_le_heatmap_sample_names, FALSE),
    show_row_names = show_row_names,
    row_labels = row_display_labels[rownames(mat)],
    row_names_gp = grid::gpar(fontsize = row_font_size),
    column_names_gp = grid::gpar(fontsize = col_font_size),
    column_names_rot = 45,
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
  add_max_deviation_line = "both",
  show_es_rank_bar = FALSE,
  show_es_le_highlight = TRUE,
  leading_edge_genes = character(0),
  le_count = NULL,
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

  # Align the display highlight with the peak/trough of the reconstructed
  # running-ES curve. Authoritative leading-edge membership is unchanged and
  # is still used for counts, heatmaps, exports, and diagnostics.
  highlight_bounds <- gsea_vis_displayed_leading_edge_bounds(running_data)
  le_xmin_es <- unname(highlight_bounds[["xmin"]])
  le_xmax_es <- unname(highlight_bounds[["xmax"]])

  plot <- ggplot2::ggplot(running_data, ggplot2::aes(x = index, y = running_es)) +
    {
      if (isTRUE(show_es_rank_bar)) {
        ggplot2::geom_tile(
          data = all_ranks,
          ggplot2::aes(x = index, y = strip_y, fill = score_limited),
          inherit.aes = FALSE,
          width = 1,
          height = 0.06 * es_span,
          show.legend = FALSE
        )
      }
    } +
    {
      if (isTRUE(show_es_rank_bar)) {
        ggplot2::scale_fill_gradient2(
          high = margins["up"],
          low = margins["dn"],
          mid = margins["md"],
          midpoint = 0,
          limits = c(-1, 1) * score_clip
        )
      }
    } +
    {
      if (isTRUE(show_es_le_highlight) && is.finite(le_xmin_es) && is.finite(le_xmax_es)) {
        ggplot2::annotate(
          geom  = "rect",
          xmin  = le_xmin_es,
          xmax  = le_xmax_es,
          ymin  = -Inf,
          ymax  = Inf,
          fill  = line_color,
          alpha = 0.08
        )
      }
    } +
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
      title = {
        t <- gsea_vis_format_es(row, running_data)
        if (!is.null(le_count) && is.finite(le_count)) {
          paste0(t, " \u00b7 ", as.integer(le_count), " LE genes")
        } else {
          t
        }
      },
      x = if (isTRUE(show_x_label)) "Gene Rank" else NULL,
      y = "Running Score"
    ) +
    {
      if (isTRUE(show_es_rank_bar)) {
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
        )
      }
    } +
    {
      if (isTRUE(show_es_rank_bar)) {
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
        )
      }
    } +
    {
      # When the gene score bar is hidden, keep +/- signs anchored to the
      # bottom of the panel (next to the rug tick carpet), in neutral black.
      if (!isTRUE(show_es_rank_bar)) {
        list(
          ggplot2::annotate(
            geom = "text",
            x = 1L,
            y = -Inf,
            label = "+",
            hjust = 1.5,
            vjust = -0.5,
            size = 2.5,
            fontface = "bold",
            color = "black"
          ),
          ggplot2::annotate(
            geom = "text",
            x = max(running_data$index),
            y = -Inf,
            label = "-",
            hjust = -0.5,
            vjust = -0.5,
            size = 2.5,
            fontface = "bold",
            color = "black"
          )
        )
      }
    } +
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
  rank_area_color = "red",
  show_rnk_le_highlight = TRUE,
  show_rnk_peak_line = TRUE,
  display_leading_edge_genes = FALSE,
  number_of_leading_edge_genes_to_display = 10,
  font_size_of_leading_edge_genes = 3
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("GSEA-Vis requires the ggplot2 package.", call. = FALSE)
  }

  leading_edge_genes <- intersect(leading_edge_genes, running_data$gene)
  ranks <- running_data
  ranks$leading_edge <- ranks$gene %in% leading_edge_genes
  ranks$LE <- ifelse(ranks$leading_edge, "LE", "Outside")
  ranks$LE <- factor(ranks$LE, levels = c("Outside", "LE"))
  ranks$display_score <- ranks$score
  score_range <- range(ranks$score, na.rm = TRUE)
  score_span <- diff(score_range)
  if (!is.finite(score_span) || score_span == 0) {
    score_span <- 1
  }
  color_margin <- line_color
  # Use the same reconstructed-curve boundary as the ES panel so both shaded
  # regions align with the displayed peak/trough and dashed guide line.
  highlight_bounds <- gsea_vis_displayed_leading_edge_bounds(running_data)
  le_xmin <- unname(highlight_bounds[["xmin"]])
  le_xmax <- unname(highlight_bounds[["xmax"]])
  area_layers <- {
    # "red_blue": per-gene gradient fill (red→white→blue by score value),
    #   identical palette to the ES strip rug — fully independent of line_color.
    # "grey"    : smooth grey area (default).
    if (grepl("^red", as.character(rank_area_color %||% "grey"), ignore.case = TRUE)) {
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
      if (isTRUE(show_rnk_le_highlight) && is.finite(le_xmin) && is.finite(le_xmax)) {
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
    ggplot2::labs(x = "Gene Rank", y = "Gene Score")

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

  if (isTRUE(show_rnk_peak_line)) {
    max_row_rnk <- running_data[running_data$max_deviation, , drop = FALSE]
    if (nrow(max_row_rnk) > 0) {
      plot <- plot + ggplot2::geom_vline(
        data = max_row_rnk,
        ggplot2::aes(xintercept = index),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 0.3,
        color = "#555555"
      )
    }
  }

  plot
}

gsea_vis_pathway_plot <- function(
  row,
  stats,
  genes,
  running_data = NULL,
  leading_edge_genes = NULL,
  expression_data = NULL,
  plots_to_include = "ES+RNK",
  running_score_line_color = "ES sign",
  add_max_deviation_line = "both",
  show_es_rank_bar = FALSE,
  show_es_le_highlight = TRUE,
  rank_area_color = "red",
  show_rnk_le_highlight = TRUE,
  show_rnk_peak_line = TRUE,
  display_leading_edge_genes = FALSE,
  number_of_leading_edge_genes_to_display = 10,
  font_size_of_leading_edge_genes = 3,
  heatmap_transform = "z-score",
  max_le_genes_heatmap = 50,
  heatmap_gene_order = "rank",
  heatmap_sample_order = "group",
  heatmap_gene_clustering_distance = "euclidean",
  heatmap_gene_clustering_method = "complete",
  heatmap_sample_clustering_distance = "euclidean",
  heatmap_sample_clustering_method = "complete",
  show_le_heatmap_sample_names = FALSE,
  show_le_heatmap_gene_names = TRUE,
  show_le_heatmap_rank_labels = TRUE,
  pdf_height = 6.5,
  pdf_width  = 8.5
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("GSEA-Vis requires the patchwork package for ES+RNK plots.", call. = FALSE)
  }

  plots_to_include <- as.character(plots_to_include %||% "ES+RNK")
  plots_to_include <- if (plots_to_include %in% c("ES", "ES+RNK", "ES+LE", "ES+RNK+LE", "LE")) {
    plots_to_include
  } else {
    "ES+RNK+LE"
  }
  line_color <- gsea_vis_plot_color(row, running_score_line_color)
  if (is.null(running_data)) {
    running_data <- gsea_vis_running_es_data(stats, genes)
  }
  if (is.null(leading_edge_genes)) {
    leading_edge_genes <- gsea_vis_row_leading_edge_genes(row, stats, running_data)
  }
  # Use the same leading-edge vector for the title count and all outputs.
  le_in_data <- intersect(leading_edge_genes, running_data$gene)
  le_count_title <- length(le_in_data)
  include_es  <- grepl("ES",  plots_to_include, fixed = TRUE)
  include_rnk <- grepl("RNK", plots_to_include, fixed = TRUE)
  include_le <- grepl("LE", plots_to_include, fixed = TRUE)

  plot_parts <- list()
  if (isTRUE(include_es)) {
    plot_parts[[length(plot_parts) + 1]] <- gsea_vis_es_plot(
      row,
      running_data,
      line_color = line_color,
      add_max_deviation_line = add_max_deviation_line,
      show_es_rank_bar = show_es_rank_bar,
      show_es_le_highlight = show_es_le_highlight,
      leading_edge_genes = leading_edge_genes,
      le_count = le_count_title,
      show_x_label = !include_rnk
    )
  }

  if (isTRUE(include_rnk)) {
    rank_plot <- gsea_vis_rank_plot(
      row,
      running_data,
      leading_edge_genes = leading_edge_genes,
      line_color = line_color,
      rank_area_color = rank_area_color,
      show_rnk_le_highlight = show_rnk_le_highlight,
      show_rnk_peak_line = show_rnk_peak_line,
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
      leading_edge_genes = leading_edge_genes,
      heatmap_transform = heatmap_transform,
      max_le_genes_heatmap = max_le_genes_heatmap,
      heatmap_gene_order = heatmap_gene_order,
      heatmap_sample_order = heatmap_sample_order,
      heatmap_gene_clustering_distance = heatmap_gene_clustering_distance,
      heatmap_gene_clustering_method = heatmap_gene_clustering_method,
      heatmap_sample_clustering_distance = heatmap_sample_clustering_distance,
      heatmap_sample_clustering_method = heatmap_sample_clustering_method,
      show_le_heatmap_sample_names = show_le_heatmap_sample_names,
      show_le_heatmap_gene_names = show_le_heatmap_gene_names,
      show_le_heatmap_rank_labels = show_le_heatmap_rank_labels,
      pdf_height = pdf_height,
      pdf_width  = pdf_width
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

gsea_vis_validate_pathway_consistency <- function(row, stats) {
  if (!is.numeric(stats) || length(stats) == 0) {
    membership <- gsea_vis_authoritative_membership(
      row,
      stats = setNames(numeric(0), character(0))
    )
    membership$valid <- FALSE
    membership$consistency_level <- "failed_missing_ranked_stats"
    membership$reason <- "Ranked DEG statistics were not available for this contrast."
    return(list(
      valid = FALSE,
      membership = membership,
      running_data = NULL,
      reconstructed_peak_es = NA_real_,
      upstream_es = gsea_vis_numeric_scalar(row[["ES"]]),
      es_direction_agrees = NA,
      reconstructed_leading_edge_genes = character(0),
      leading_edge_boundary_agrees = NA,
      leading_edge_overlap = NA_real_,
      missing_from_reconstructed_leading_edge = character(0),
      extra_in_reconstructed_leading_edge = character(0),
      boundary_warning_message = NULL,
      reason = membership$reason
    ))
  }

  membership <- gsea_vis_authoritative_membership(row, stats)
  if (!isTRUE(membership$valid)) {
    return(list(
      valid = FALSE,
      membership = membership,
      running_data = NULL,
      reconstructed_peak_es = NA_real_,
      upstream_es = gsea_vis_numeric_scalar(row[["ES"]]),
      es_direction_agrees = NA,
      reconstructed_leading_edge_genes = character(0),
      leading_edge_boundary_agrees = NA,
      leading_edge_overlap = NA_real_,
      missing_from_reconstructed_leading_edge = character(0),
      extra_in_reconstructed_leading_edge = character(0),
      boundary_warning_message = NULL,
      reason = membership$reason
    ))
  }

  running_data <- tryCatch(
    gsea_vis_running_es_data(stats, membership$pathway_genes_ranked),
    error = function(error) error
  )
  if (inherits(running_data, "error")) {
    membership$valid <- FALSE
    membership$consistency_level <- "failed_running_es"
    membership$reason <- conditionMessage(running_data)
    return(list(
      valid = FALSE,
      membership = membership,
      running_data = NULL,
      reconstructed_peak_es = NA_real_,
      upstream_es = gsea_vis_numeric_scalar(row[["ES"]]),
      es_direction_agrees = NA,
      reconstructed_leading_edge_genes = character(0),
      leading_edge_boundary_agrees = NA,
      leading_edge_overlap = NA_real_,
      missing_from_reconstructed_leading_edge = character(0),
      extra_in_reconstructed_leading_edge = character(0),
      boundary_warning_message = NULL,
      reason = membership$reason
    ))
  }

  peak_rows <- which(running_data$max_deviation)
  reconstructed_peak_es <- if (length(peak_rows) > 0) {
    suppressWarnings(as.numeric(running_data$running_es[[peak_rows[[1]]]]))
  } else {
    NA_real_
  }
  upstream_es <- gsea_vis_numeric_scalar(row[["ES"]])
  if (!is.finite(upstream_es)) {
    upstream_es <- gsea_vis_numeric_scalar(row[["NES"]])
  }
  tolerance <- sqrt(.Machine$double.eps)
  direction_comparable <- is.finite(upstream_es) && is.finite(reconstructed_peak_es) &&
    abs(upstream_es) > tolerance && abs(reconstructed_peak_es) > tolerance
  es_direction_agrees <- if (direction_comparable) {
    sign(upstream_es) == sign(reconstructed_peak_es)
  } else {
    NA
  }

  if (isFALSE(es_direction_agrees)) {
    membership$valid <- FALSE
    membership$consistency_level <- "failed_es_direction_mismatch"
    membership$reason <- sprintf(
      "Reconstructed running ES direction (%s) disagreed with the upstream GSEA ES direction (%s).",
      if (reconstructed_peak_es < 0) "negative" else "positive",
      if (upstream_es < 0) "negative" else "positive"
    )
    return(list(
      valid = FALSE,
      membership = membership,
      running_data = running_data,
      reconstructed_peak_es = reconstructed_peak_es,
      upstream_es = upstream_es,
      es_direction_agrees = FALSE,
      reconstructed_leading_edge_genes = gsea_vis_reconstructed_leading_edge_genes(running_data),
      leading_edge_boundary_agrees = NA,
      leading_edge_overlap = NA_real_,
      missing_from_reconstructed_leading_edge = character(0),
      extra_in_reconstructed_leading_edge = character(0),
      boundary_warning_message = NULL,
      reason = membership$reason
    ))
  }

  reconstructed_leading_edge_genes <- gsea_vis_reconstructed_leading_edge_genes(running_data)
  authoritative_leading_edge_genes <- membership$leading_edge_genes_ranked
  missing_from_reconstructed <- setdiff(
    authoritative_leading_edge_genes,
    reconstructed_leading_edge_genes
  )
  extra_in_reconstructed <- setdiff(
    reconstructed_leading_edge_genes,
    authoritative_leading_edge_genes
  )
  leading_edge_boundary_agrees <- setequal(
    authoritative_leading_edge_genes,
    reconstructed_leading_edge_genes
  )
  leading_edge_overlap <- if (length(authoritative_leading_edge_genes) > 0) {
    sum(authoritative_leading_edge_genes %in% reconstructed_leading_edge_genes) /
      length(authoritative_leading_edge_genes)
  } else {
    NA_real_
  }
  # Retain boundary-comparison diagnostics in the consistency table, but do
  # not emit a runtime warning for this expected legacy-ranking mismatch.
  boundary_warning_message <- NULL

  list(
    valid = TRUE,
    membership = membership,
    running_data = running_data,
    reconstructed_peak_es = reconstructed_peak_es,
    upstream_es = upstream_es,
    es_direction_agrees = es_direction_agrees,
    reconstructed_leading_edge_genes = reconstructed_leading_edge_genes,
    leading_edge_boundary_agrees = leading_edge_boundary_agrees,
    leading_edge_overlap = leading_edge_overlap,
    missing_from_reconstructed_leading_edge = missing_from_reconstructed,
    extra_in_reconstructed_leading_edge = extra_in_reconstructed,
    boundary_warning_message = boundary_warning_message,
    reason = membership$reason
  )
}

gsea_vis_consistency_diagnostic <- function(
  row,
  validation,
  heatmap_leading_edge_total = NA_integer_,
  heatmap_leading_edge_matched = NA_integer_,
  missing_heatmap_leading_edge_genes = character(0)
) {
  membership <- validation$membership
  data.frame(
    contrast = as.character(row$contrast),
    collection = as.character(row$collection),
    pathway = as.character(row$pathway),
    consistency_level = as.character(membership$consistency_level),
    pathway_membership_field = as.character(membership$pathway_field),
    leading_edge_field = as.character(membership$leading_edge_field),
    pathway_genes_total = as.integer(membership$pathway_total),
    pathway_genes_matched = as.integer(membership$pathway_matched),
    pathway_overlap = as.numeric(membership$pathway_overlap),
    leading_edge_total = as.integer(membership$leading_edge_total),
    leading_edge_matched = as.integer(membership$leading_edge_matched),
    upstream_ES = as.numeric(validation$upstream_es),
    reconstructed_peak_ES = as.numeric(validation$reconstructed_peak_es),
    es_direction_agrees = as.logical(validation$es_direction_agrees),
    reconstructed_leading_edge_total = as.integer(length(validation$reconstructed_leading_edge_genes %||% character(0))),
    leading_edge_boundary_agrees = as.logical(validation$leading_edge_boundary_agrees),
    authoritative_leading_edge_overlap_with_reconstructed = as.numeric(validation$leading_edge_overlap),
    missing_from_reconstructed_leading_edge = paste(
      validation$missing_from_reconstructed_leading_edge %||% character(0),
      collapse = ";"
    ),
    extra_in_reconstructed_leading_edge = paste(
      validation$extra_in_reconstructed_leading_edge %||% character(0),
      collapse = ";"
    ),
    heatmap_leading_edge_total = as.integer(heatmap_leading_edge_total),
    heatmap_leading_edge_matched = as.integer(heatmap_leading_edge_matched),
    missing_pathway_genes = paste(membership$missing_pathway_genes, collapse = ";"),
    missing_leading_edge_genes = paste(membership$missing_leading_edge_genes, collapse = ";"),
    missing_heatmap_leading_edge_genes = paste(missing_heatmap_leading_edge_genes, collapse = ";"),
    reason = as.character(validation$reason),
    stringsAsFactors = FALSE
  )
}

gsea_vis_running_es_output <- function(row, running_data, leading_edge_genes) {
  # The running score is calculated across the full ranked gene list, but the
  # exported table contains only genes that belong to the current pathway.
  pathway_data <- running_data[
    !is.na(running_data$in_pathway) & running_data$in_pathway,
    ,
    drop = FALSE
  ]

  is_leading_edge <- pathway_data$gene %in% as.character(
    leading_edge_genes %||% character(0)
  )

  data.frame(
    contrast = rep(as.character(row$contrast), nrow(pathway_data)),
    collection = rep(as.character(row$collection), nrow(pathway_data)),
    pathway = rep(as.character(row$pathway), nrow(pathway_data)),
    geneRank = pathway_data$index,
    runningES = pathway_data$running_es,
    gene = pathway_data$gene,
    isLeadingEdge = is_leading_edge,
    geneScore = pathway_data$score,
    stringsAsFactors = FALSE
  )
}

gsea_vis_manifest_row <- function(row, genes, plot_id, pdf_page, genes_in_heatmap = 0L) {
  data.frame(
    plot_id = plot_id,
    pdf_page = pdf_page,
    contrast = as.character(row$contrast),
    collection = as.character(row$collection),
    pathway = as.character(row$pathway),
    NES = gsea_vis_numeric_scalar(row[["NES"]]),
    ES = gsea_vis_numeric_scalar(row[["ES"]]),
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
  batch_result = NULL,
  plot_contrasts = character(0),
  plot_all_pathways = FALSE,
  top_n_pathways = 1,
  top_n_by_sign = FALSE,
  max_plots_in_pdf = 50,
  stop_if_too_many_plots = TRUE,
  plots_to_include = "ES+RNK+LE",
  running_score_line_color = "ES sign",
  add_max_deviation_line = "both",
  rank_plot_mode = NULL,
  show_es_rank_bar = FALSE,
  show_es_le_highlight = TRUE,
  rank_area_color = "red",
  show_rnk_le_highlight = TRUE,
  show_rnk_peak_line = TRUE,
  display_leading_edge_genes = FALSE,
  number_of_leading_edge_genes_to_display = 10,
  font_size_of_leading_edge_genes = 3,
  heatmap_gene_names_column = "gene",
  heatmap_sample_names_column = "Sample",
  heatmap_group_column = "Group",
  heatmap_transform = "z-score",
  max_le_genes_heatmap = 50,
  heatmap_gene_order = "rank",
  heatmap_sample_order = "group",
  heatmap_gene_clustering_distance = "euclidean",
  heatmap_gene_clustering_method = "complete",
  heatmap_sample_clustering_distance = "euclidean",
  heatmap_sample_clustering_method = "complete",
  show_le_heatmap_sample_names = FALSE,
  show_le_heatmap_gene_names = TRUE,
  show_le_heatmap_rank_labels = TRUE,
  pdf_width = 8.5,
  pdf_height = 6.5,
  output_dir = file.path(getwd(), "_workflow_runtime", "outputs", "gsea_vis"),
  selected_rows = NULL
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
  if (is.null(selected_rows)) {
    selected <- gsea_vis_select_rows(
      gsea_table,
      plot_contrasts = plot_contrasts,
      plot_all_pathways = isTRUE(plot_all_pathways),
      top_n_pathways = top_n_pathways,
      top_n_by_sign = isTRUE(top_n_by_sign)
    )
  } else {
    if (!is.data.frame(selected_rows)) {
      stop("GSEA-Vis selected_rows must be a data frame when provided.", call. = FALSE)
    }
    selected <- selected_rows
  }
  if (nrow(selected) == 0) {
    return(list(
      manifest = data.frame(),
      plots = list(),
      message = "GSEA-Vis did not generate plots because no rows matched the selected contrasts."
    ))
  }

  max_plots_value <- suppressWarnings(as.integer(max_plots_in_pdf %||% 0L))
  has_plot_cap <- length(max_plots_value) > 0 && !is.na(max_plots_value[[1]]) && max_plots_value[[1]] > 0L
  if (isTRUE(has_plot_cap) && nrow(selected) > max_plots_value[[1]] && isTRUE(stop_if_too_many_plots)) {
    stop(
      "GSEA-Vis selected ",
      nrow(selected),
      " plots, which is above Max plots in PDF (",
      max_plots_value[[1]],
      "). Narrow the GSEA-F filter or increase the cap.",
      call. = FALSE
    )
  }
  if (isTRUE(has_plot_cap) && nrow(selected) > max_plots_value[[1]]) {
    selected <- utils::head(selected, max_plots_value[[1]])
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(output_dir, pattern = "^GSEA-Vis-", full.names = TRUE), force = TRUE)

  expression_data <- NULL
  expression_by_contrast <- list()
  if (gsea_vis_needs_le_heatmap(plots_to_include)) {
    expression_data <- gsea_vis_expression_data(
      batch_result,
      gene_names_column = heatmap_gene_names_column,
      sample_names_column = heatmap_sample_names_column,
      heatmap_group_column = heatmap_group_column
    )
    if (is.null(expression_data)) {
      stop(
        "GSEA-Vis needs normalized expression and sample metadata to draw leading-edge heatmaps.",
        call. = FALSE
      )
    }

    selected_contrasts <- unique(as.character(selected$contrast))
    selected_contrasts <- selected_contrasts[!is.na(selected_contrasts) & nzchar(selected_contrasts)]
    expression_by_contrast <- setNames(
      lapply(selected_contrasts, function(contrast) {
        representative_row <- selected[match(contrast, as.character(selected$contrast)), , drop = FALSE]
        gsea_vis_filter_expression_for_contrast(expression_data, representative_row)
      }),
      selected_contrasts
    )
  }

  consistency_file <- file.path(output_dir, "GSEA-Vis-Input-Consistency.csv")
  skipped_file <- file.path(output_dir, "GSEA-Vis-Skipped.csv")
  validation_entries <- vector("list", nrow(selected))
  diagnostic_rows <- vector("list", nrow(selected))

  for (row_index in seq_len(nrow(selected))) {
    row <- selected[row_index, , drop = FALSE]
    contrast <- as.character(row$contrast)
    contrast <- if (length(contrast) > 0) contrast[[1]] else ""
    validation <- gsea_vis_validate_pathway_consistency(
      row,
      ranked_stats[[contrast]]
    )
    validation_entries[[row_index]] <- list(
      row = row,
      contrast = contrast,
      validation = validation
    )
    diagnostic_rows[[row_index]] <- gsea_vis_consistency_diagnostic(row, validation)
  }

  consistency_table <- do.call(rbind, diagnostic_rows)
  failed_validation <- !vapply(
    validation_entries,
    function(entry) isTRUE(entry$validation$valid),
    logical(1)
  )
  failed_count <- sum(failed_validation)
  selected_count <- length(validation_entries)

  # A strict majority of failed pathways indicates a systemic input mismatch.
  # Write the diagnostics first so the cause remains auditable in /results.
  if (failed_count > selected_count / 2) {
    #utils::write.csv(consistency_table, consistency_file, row.names = FALSE)
    failed_table <- consistency_table[failed_validation, , drop = FALSE]
    #utils::write.csv(failed_table, skipped_file, row.names = FALSE)
    stop(
      sprintf(
        paste0(
          "Input consistency checks failed for %d of %d selected pathways. ",
          "This suggests that the DEG ranking and GSEA Filter results do not come from a compatible analysis. ",
          "See %s."
        ),
        failed_count,
        selected_count,
        normalizePath(consistency_file, winslash = "/", mustWork = FALSE)
      ),
      call. = FALSE
    )
  }

  for (entry in validation_entries) {
    validation <- entry$validation
    pathway_label <- paste0(entry$contrast, " / ", as.character(entry$row$pathway))
    if (!isTRUE(validation$valid)) {
      message("WARNING: ", pathway_label, " — skipped: ", validation$reason)
    } else {
      if (identical(validation$membership$consistency_level, "leading_edge_strict")) {
        message("WARNING: ", pathway_label, " — ", validation$membership$warning_message)
      }
    }
  }

  manifest_rows <- list()
  skipped <- character(0)
  pdf_file <- file.path(output_dir, "GSEA-Vis-Enrichment-Plots.pdf")
  running_es_file <- file.path(output_dir, "GSEA-Vis-RunningES.csv")

  pdf_device <- NULL
  running_es_connection <- NULL
  running_es_header_written <- FALSE
  grDevices::pdf(pdf_file, width = as.numeric(pdf_width), height = as.numeric(pdf_height))
  pdf_device <- grDevices::dev.cur()
  on.exit({
    if (!is.null(running_es_connection)) {
      try(close(running_es_connection), silent = TRUE)
    }
    if (!is.null(pdf_device) && pdf_device %in% grDevices::dev.list()) {
      try(grDevices::dev.off(which = pdf_device), silent = TRUE)
    }
  }, add = TRUE)
  running_es_connection <- file(running_es_file, open = "wt")

  for (row_index in seq_along(validation_entries)) {
    entry <- validation_entries[[row_index]]
    validation <- entry$validation
    row <- entry$row
    contrast <- entry$contrast
    if (!isTRUE(validation$valid)) {
      skipped <- c(
        skipped,
        paste0(contrast, " / ", as.character(row$pathway), ": ", validation$reason)
      )
      next
    }

    membership <- validation$membership
    stats <- ranked_stats[[contrast]]
    genes <- membership$pathway_genes_ranked
    leading_edge_genes <- membership$leading_edge_genes_ranked
    running_data <- validation$running_data
    row_expression_data <- if (!is.null(expression_data)) {
      expression_by_contrast[[contrast]] %||% expression_data
    } else {
      NULL
    }

    heatmap_leading_edge_genes <- character(0)
    missing_heatmap_leading_edge_genes <- character(0)
    if (!is.null(row_expression_data) && is.matrix(row_expression_data$matrix)) {
      heatmap_mapping <- gsea_vis_map_genes_to_reference(
        leading_edge_genes,
        row = row,
        reference_genes = rownames(row_expression_data$matrix)
      )
      heatmap_leading_edge_genes <- unique(unname(
        heatmap_mapping[!is.na(heatmap_mapping) & nzchar(heatmap_mapping)]
      ))
      missing_heatmap_leading_edge_genes <- unique(
        names(heatmap_mapping)[is.na(heatmap_mapping) | !nzchar(heatmap_mapping)]
      )
      if (length(missing_heatmap_leading_edge_genes) > 0) {
        message(
          "WARNING: ",
          contrast,
          " / ",
          as.character(row$pathway),
          " — ",
          length(heatmap_leading_edge_genes),
          " of ",
          length(leading_edge_genes),
          " authoritative leading-edge genes matched the expression matrix."
        )
      }
    }

    diagnostic_rows[[row_index]] <- gsea_vis_consistency_diagnostic(
      row,
      validation,
      heatmap_leading_edge_total = if (!is.null(row_expression_data)) length(leading_edge_genes) else NA_integer_,
      heatmap_leading_edge_matched = if (!is.null(row_expression_data)) length(heatmap_leading_edge_genes) else NA_integer_,
      missing_heatmap_leading_edge_genes = missing_heatmap_leading_edge_genes
    )

    plot_id <- paste(
      gsea_vis_safe_file_name(contrast),
      gsea_vis_safe_file_name(row$pathway),
      sep = "__"
    )
    plot <- gsea_vis_pathway_plot(
      row,
      stats = stats,
      genes = genes,
      running_data = running_data,
      leading_edge_genes = leading_edge_genes,
      expression_data = row_expression_data,
      plots_to_include = plots_to_include,
      running_score_line_color = running_score_line_color,
      add_max_deviation_line = add_max_deviation_line,
      show_es_rank_bar = isTRUE(show_es_rank_bar),
      show_es_le_highlight = isTRUE(show_es_le_highlight),
      rank_area_color = rank_area_color,
      show_rnk_le_highlight = isTRUE(show_rnk_le_highlight),
      show_rnk_peak_line = isTRUE(show_rnk_peak_line),
      display_leading_edge_genes = isTRUE(display_leading_edge_genes),
      number_of_leading_edge_genes_to_display = number_of_leading_edge_genes_to_display,
      font_size_of_leading_edge_genes = font_size_of_leading_edge_genes,
      heatmap_transform = heatmap_transform,
      max_le_genes_heatmap = max_le_genes_heatmap,
      heatmap_gene_order = heatmap_gene_order,
      heatmap_sample_order = heatmap_sample_order,
      heatmap_gene_clustering_distance = heatmap_gene_clustering_distance,
      heatmap_gene_clustering_method = heatmap_gene_clustering_method,
      heatmap_sample_clustering_distance = heatmap_sample_clustering_distance,
      heatmap_sample_clustering_method = heatmap_sample_clustering_method,
      show_le_heatmap_sample_names = gsea_vis_bool(show_le_heatmap_sample_names, FALSE),
      show_le_heatmap_gene_names = gsea_vis_bool(show_le_heatmap_gene_names, TRUE),
      show_le_heatmap_rank_labels = gsea_vis_bool(show_le_heatmap_rank_labels, TRUE),
      pdf_height = as.numeric(pdf_height),
      pdf_width = as.numeric(pdf_width)
    )
    print(plot)

    manifest_rows[[length(manifest_rows) + 1]] <- gsea_vis_manifest_row(
      row,
      genes = genes,
      plot_id = plot_id,
      pdf_page = length(manifest_rows) + 1L,
      genes_in_heatmap = length(heatmap_leading_edge_genes)
    )

    running_es_output <- gsea_vis_running_es_output(
      row,
      running_data,
      leading_edge_genes = leading_edge_genes
    )
    utils::write.table(
      running_es_output,
      file = running_es_connection,
      sep = ",",
      row.names = FALSE,
      col.names = !running_es_header_written,
      quote = TRUE,
      qmethod = "double"
    )
    running_es_header_written <- TRUE
  }

  grDevices::dev.off(which = pdf_device)
  pdf_device <- NULL
  close(running_es_connection)
  running_es_connection <- NULL

  consistency_table <- do.call(rbind, diagnostic_rows)
  #utils::write.csv(consistency_table, consistency_file, row.names = FALSE)

  failed_table <- consistency_table[grepl("^failed_", consistency_table$consistency_level), , drop = FALSE]
  skipped_file_out <- NULL
  if (nrow(failed_table) > 0) {
    #utils::write.csv(failed_table, skipped_file, row.names = FALSE)
    skipped_file_out <- skipped_file
  }

  manifest <- if (length(manifest_rows) > 0) {
    do.call(rbind, manifest_rows)
  } else {
    data.frame(
      plot_id = character(0),
      pdf_page = integer(0),
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

  if (nrow(manifest) == 0) {
    if (file.exists(pdf_file)) unlink(pdf_file, force = TRUE)
    if (file.exists(running_es_file)) unlink(running_es_file, force = TRUE)
    stop(
      paste0(
        "GSEA-Vis generated no plots. ",
        if (length(skipped) > 0) {
          paste0("All selected pathways were skipped. First reason: ", skipped[[1]])
        } else {
          "No selected pathway produced plottable data."
        }
      ),
      call. = FALSE
    )
  }

  consistency_summary <- list(
    all_genes_strict = sum(consistency_table$consistency_level == "all_genes_strict"),
    leading_edge_strict = sum(consistency_table$consistency_level == "leading_edge_strict"),
    boundary_mismatch = sum(consistency_table$leading_edge_boundary_agrees %in% FALSE, na.rm = TRUE),
    failed = sum(grepl("^failed_", consistency_table$consistency_level))
  )

  list(
    manifest = manifest,
    running_es = NULL,
    plots = list(
      pdf = if (file.exists(pdf_file)) normalizePath(pdf_file, winslash = "/", mustWork = FALSE) else NULL
    ),
    files = list(
      pdf = if (file.exists(pdf_file)) normalizePath(pdf_file, winslash = "/", mustWork = FALSE) else NULL,
      running_es = if (file.exists(running_es_file)) normalizePath(running_es_file, winslash = "/", mustWork = FALSE) else NULL,
      consistency = normalizePath(consistency_file, winslash = "/", mustWork = FALSE),
      skipped = if (!is.null(skipped_file_out)) normalizePath(skipped_file_out, winslash = "/", mustWork = FALSE) else NULL
    ),
    consistency = consistency_summary,
    skipped = skipped,
    message = paste0(
      "GSEA-Vis complete: generated ",
      nrow(manifest),
      " ",
      as.character(plots_to_include %||% "ES+RNK"),
      " plot(s)",
      if (consistency_summary$leading_edge_strict > 0) {
        paste0("; ", consistency_summary$leading_edge_strict, " partial pathway match(es) with complete leading edge")
      } else {
        ""
      },
      if (length(skipped) > 0) paste0("; skipped ", length(skipped)) else "",
      "."
    )
  )
}

