#' Plot a tSNE or umap colored by feature.
#'
#' @param data input data
#' @param x x variable
#' @param y y variable
#' @param feature feature to color by
#' @param legend_name legend name to display, defaults to no name
#' @param c_cols character vector of colors to built color gradient
#' for continuous values. defaults to [`clustifyr::pretty_palette`]
#' @param d_cols character vector of colors for discrete values.
#' defaults to RColorBrewer paired palette
#' @param pt_size point size
#' @param scale_limits defaults to min = 0, max = max(data$x),
#' otherwise a two-element numeric vector indicating min and max to plot
#' @param do_label whether to label each cluster at median center
#' @param do_legend whether to draw legend
#' @return ggplot object, cells projected by dr, colored by feature
#' @export
plot_tsne <- function(data, x = "UMAP_1", y = "UMAP_2",
                      feature,
                      legend_name = "",
                      c_cols = pretty_palette2,
                      d_cols = NULL,
                      pt_size = 0.25,
                      scale_limits = NULL,
                      do_label = FALSE,
                      do_legend = TRUE) {
  if ((length(unique(data[[feature]])) > 12) & identical(c_cols, pretty_palette2) & (typeof(data[[feature]]) %in% c(
    "character",
    "logical",
    "factor"
  ))) {
    c_cols <- pretty_palette_ramp_d(length(unique(data[[feature]])))
    d_cols <- pretty_palette_ramp_d(length(unique(data[[feature]])))
  }
  # sort data to avoid plotting null values over colors
  data <- dplyr::arrange(data, !!dplyr::sym(feature))

  # Add backticks to allow special characters in column names
  x <- paste0("`", x, "`")
  y <- paste0("`", y, "`")

  p <- ggplot2::ggplot(data, ggplot2::aes_string(x, y)) +
    geom_point(ggplot2::aes_string(color = paste0("`", feature, "`")),
      size = pt_size
    )

  if (typeof(data[[feature]]) %in% c(
    "character",
    "logical"
  ) | is.factor(data[[feature]])) {
    if (!is.null(d_cols)) {
      # use colors provided
      p <- p + scale_color_manual(
        values = d_cols,
        name = legend_name
      )
    } else {
      p <- p + scale_color_brewer(
        palette = "Paired",
        name = legend_name
      )
    }
  } else {
    # continuous values
    if (is.null(scale_limits)) {
      scale_limits <- c(
        ifelse(min(data[[feature]]) < 0,
          min(data[[feature]]),
          0
        ),
        max(data[[feature]])
      )
    }
    p <- p + scale_color_gradientn(
      colors = c_cols,
      name = legend_name,
      limits = scale_limits
    )
  }

  if (do_label) {
    centers <- dplyr::group_by_(data, .dots = feature)
    centers <- dplyr::summarize(centers, t1 = mean(!!dplyr::sym(x)), t2 = mean(!!dplyr::sym(y)))
    p <- p +
      geom_point(data = centers, mapping = aes(x = !!dplyr::sym("t1"), y = !!dplyr::sym("t2")), size = 0, alpha = 0) +
      geom_text(data = centers, mapping = aes(x = !!dplyr::sym("t1"), y = !!dplyr::sym("t2"), label = centers[[feature]]))
  }

  p <- p + cowplot::theme_cowplot()

  if (!do_legend) {
    p <- p + theme(legend.position = "none")
  }

  p
}

#' Color palette for plotting continous variables
#' @return vector of colors
#' @export
pretty_palette <- rev(scales::brewer_pal(palette = "RdGy")(6))

#' Color palette for plotting continous variables, starting at gray
#' @return vector of colors
#' @export
pretty_palette2 <- scales::brewer_pal(palette = "Reds")(9)

#' black and white palette for plotting continous variables
#' @return vector of colors
#' @export
not_pretty_palette <- scales::brewer_pal(palette = "Greys")(9)

#' Expanded color palette ramp for plotting discrete variables
#' @param n number of colors to use
#' @return color ramp
#' @export
pretty_palette_ramp_d <- grDevices::colorRampPalette(scales::brewer_pal(palette = "Paired")(12))

#' Plot similarity measures on a tSNE or umap
#'
#' @param cor_matrix input similarity matrix
#' @param metadata input metadata with per cell tsne or umap coordinates and cluster ids
#' @param data_to_plot colname of data to plot, defaults to all
#' @param cluster_col colname of clustering data in metadata, defaults to rownames of the
#' metadata if not supplied.
#' @param x metadata column name with 1st axis dimension.
#' defaults to "UMAP_1".
#' @param y metadata column name with 2nd axis dimension.
#' defaults to "UMAP_2".
#' @param scale_legends if TRUE scale all legends to maximum values in entire
#' correlation matrix. if FALSE scale legends to maximum for each plot. A
#' two-element numeric vector can also be passed to supply custom values i.e. c(0, 1)
#' @param ... passed to plot_tsne
#' @return list of ggplot objects, cells projected by dr, colored by cor values
#' @examples
#' res <- clustify(
#'   input = pbmc_matrix_small,
#'   metadata = pbmc_meta,
#'   ref_mat = pbmc_bulk_matrix,
#'   query_genes = pbmc_vargenes,
#'   cluster_col = "classified"
#' )
#'
#' plts <- plot_cor(
#'   cor_matrix = res,
#'   metadata = pbmc_meta,
#'   data_to_plot = colnames(res)[1:2],
#'   cluster_col = "classified",
#'   x = "UMAP_1",
#'   y = "UMAP_2"
#' )
#'
#' plts
#' @export
plot_cor <- function(cor_matrix,
                     metadata,
                     data_to_plot = colnames(cor_matrix),
                     cluster_col = NULL,
                     x = "UMAP_1",
                     y = "UMAP_2",
                     scale_legends = FALSE,
                     ...) {
  if (!any(data_to_plot %in% colnames(cor_matrix))) {
    stop("cluster ids not shared between metadata and correlation matrix", call. = FALSE)
  }

  if (is.null(cluster_col)) {
    cluster_col <- "rownames"
    metadata <- tibble::rownames_to_column(metadata, cluster_col)
  }

  cor_df <- as.data.frame(cor_matrix)
  cor_df <- tibble::rownames_to_column(cor_df, cluster_col)
  cor_df_long <- tidyr::gather(
    cor_df,
    "ref_cluster",
    "expr",
    -dplyr::matches(cluster_col)
  )

  # If cluster_col is factor, convert to character
  if (is.factor(metadata[, cluster_col])) {
    metadata[, cluster_col] <- as.character(metadata[, cluster_col])
  }

  # checks matrix rownames, 2 branches for cluster number (avg) or cell bar code (each cell)
  if (cor_df[[cluster_col]][1] %in% metadata[[cluster_col]]) {
    plt_data <- dplyr::left_join(cor_df_long,
      metadata,
      by = cluster_col
    )
  } else {
    plt_data <- dplyr::left_join(cor_df_long,
      metadata,
      by = structure(names = cluster_col, "rn")
    )
  }

  # determine scaling method, either same for all plots, or per plot (default)
  if (typeof(scale_legends) == "logical" && scale_legends) {
    scale_limits <- c(
      ifelse(min(plt_data$expr) < 0,
        min(plt_data$expr),
        0
      ),
      max(max(plt_data$expr))
    )
  } else if (typeof(scale_legends) == "logical" && !scale_legends) {
    scale_limits <- NULL
  } else {
    scale_limits <- scale_legends
  }

  plts <- vector("list", length(data_to_plot))

  for (i in seq_along(data_to_plot)) {
    tmp_data <- dplyr::filter(
      plt_data,
      !!dplyr::sym("ref_cluster") == data_to_plot[i]
    )
    plts[[i]] <- plot_tsne(tmp_data,
      x = x,
      y = y,
      feature = "expr",
      legend_name = data_to_plot[i],
      scale_limits = scale_limits,
      ...
    )
  }
  plts
}

#' Plot gene expression on to tSNE or umap
#'
#' @param expr_mat input single cell matrix
#' @param metadata data.frame with tSNE or umap coordinates
#' @param genes gene(s) to color tSNE or umap
#' @param cell_col column name in metadata containing cell ids, defaults
#' to rownames if not supplied
#' @param ... additional arguments passed to `[clustifyr::plot_tsne()]`
#' @return list of ggplot object, cells projected by dr, colored by gene expression
#' @examples
#' genes <- c(
#'   "RP11-314N13.3",
#'   "ARF4"
#' )
#'
#' plts <- plot_gene(
#'   expr_mat = pbmc_matrix_small,
#'   metadata = tibble::rownames_to_column(pbmc_meta, "rn"),
#'   genes = genes,
#'   cell_col = "rn"
#' )
#'
#' plts
#' @export
plot_gene <- function(expr_mat,
                      metadata,
                      genes,
                      cell_col = NULL,
                      ...) {
  genes_to_plot <- genes[genes %in% rownames(expr_mat)]
  genes_missing <- setdiff(genes, genes_to_plot)

  if (length(genes_missing) != 0) {
    warning(paste0(
      "the following genes were not present in the input matrix ",
      paste(genes_missing, collapse = ",")
    ))
  }

  if (length(genes_to_plot) == 0) {
    stop("no genes present to plot", call. = FALSE)
  }
  expr_dat <- t(as.matrix(expr_mat[genes_to_plot, , drop = FALSE]))
  expr_dat <- tibble::rownames_to_column(as.data.frame(expr_dat), "cell")

  if (is.null(cell_col)) {
    mdata <- tibble::rownames_to_column(metadata, "cell")
    cell_col <- "cell"
  } else {
    mdata <- metadata
  }

  if (!cell_col %in% colnames(mdata)) {
    stop("please supply a cell_col that is present in metadata", call. = FALSE)
  }

  plt_dat <- dplyr::left_join(expr_dat, mdata,
    by = c("cell" = cell_col)
  )

  lapply(
    genes_to_plot,
    function(gene) {
      plot_tsne(plt_dat,
        feature = gene,
        legend_name = gene,
        ...
      )
    }
  )
}

#' Plot called clusters on a tSNE or umap, for each reference cluster given
#'
#' @param cor_matrix input similarity matrix
#' @param metadata input metadata with tsne or umap coordinates and cluster ids
#' @param data_to_plot colname of data to plot, defaults to all
#' @param ... passed to plot_tsne
#' @return list of ggplot object, cells projected by dr, colored by cell type classification
#' @examples
#' res <- clustify(
#'   input = pbmc_matrix_small,
#'   metadata = pbmc_meta,
#'   ref_mat = pbmc_bulk_matrix,
#'   query_genes = pbmc_vargenes,
#'   cluster_col = "classified"
#' )
#'
#' plts <- plot_call(
#'   cor_matrix = res,
#'   metadata = pbmc_meta,
#'   data_to_plot = colnames(res)[1:2],
#'   cluster_col = "classified"
#' )
#'
#' plts
#' @export
plot_call <- function(cor_matrix,
                      metadata,
                      data_to_plot = colnames(cor_matrix),
                      ...) {
  df_temp <- as.data.frame(t(apply(cor_matrix, 1, function(x) x - max(x))))
  df_temp[df_temp == 0] <- "1"
  df_temp[df_temp != "1"] <- "0"
  plot_cor(
    df_temp,
    metadata,
    data_to_plot,
    ...
  )
}

#' Plot best calls for each cluster on a tSNE or umap
#'
#' @param cor_matrix input similarity matrix
#' @param metadata input metadata with tsne or umap coordinates and cluster ids
#' @param cluster_col metadata column, can be cluster or cellid
#' @param collapse_to_cluster if a column name is provided, takes the most frequent call of entire cluster to color in plot
#' @param threshold minimum correlation coefficent cutoff for calling clusters
#' @param x x variable
#' @param y y variable
#' @param plot_r whether to include second plot of cor eff for best call
#' @param per_cell whether the cor_matrix was generate per cell or per cluster
#' @param ... passed to plot_tsne
#' @return ggplot object, cells projected by dr, colored by cell type classification
#' @examples
#' res <- clustify(
#'   input = pbmc_matrix_small,
#'   metadata = pbmc_meta,
#'   ref_mat = pbmc_bulk_matrix,
#'   query_genes = pbmc_vargenes,
#'   cluster_col = "classified"
#' )
#'
#' plts <- plot_best_call(
#'   cor_matrix = res,
#'   metadata = pbmc_meta,
#'   cluster_col = "classified"
#' )
#'
#' plts
#' @export
plot_best_call <- function(cor_matrix,
                           metadata,
                           cluster_col = "cluster",
                           collapse_to_cluster = FALSE,
                           threshold = 0,
                           x = "UMAP_1", y = "UMAP_2",
                           plot_r = FALSE,
                           per_cell = FALSE,
                           ...) {
  col_meta <- colnames(metadata)
  if ("type" %in% col_meta | "type2" %in% col_meta) {
    warning('metadata column name clash of "type"/"type2"')
    return()
  }
  df_temp <- cor_to_call(
    cor_matrix,
    metadata = metadata,
    cluster_col = cluster_col,
    threshold = threshold
  )

  df_temp_full <<- call_to_metadata(
    df_temp,
    metadata = metadata,
    cluster_col = cluster_col,
    per_cell = per_cell
  )

  if (collapse_to_cluster != FALSE) {
    df_temp_full <- collapse_to_cluster(df_temp_full,
      metadata,
      collapse_to_cluster,
      threshold = threshold
    )
  }

  g <- plot_tsne(df_temp_full,
    feature = "type",
    x = x, y = y,
    ...
  )

  if (plot_r) {
    l <- list()
    l[[1]] <- g
    l[[2]] <- plot_tsne(df_temp_full,
      feature = "r",
      x = x, y = y,
      ...
    )
  } else {
    l <- g
  }

  l
}

#' Plot variable median per cluster from reference metadata vs new assigned metadata, for visually evaluating classification
#'
#' @param metadata clustify-ed metadata
#' @param cluster_col metadata column for original clustering
#' @param cluster_col_called metadata column for clustifyr assignments
#' @param plot_col metadata column to plot
#' @param metadata_ref reference single cell RNA seq metadata
#' @param cluster_col_ref metadata column for cluster for reference
#' @param plot_col_ref metadata column to plot for reference
#' @return ggplot scatterplot
#' @examples
#' plt <- plot_cols(
#'   metadata = pbmc_meta,
#'   cluster_col = "seurat_clusters",
#'   cluster_col_called = "classified",
#'   plot_col = "UMAP_1",
#'   metadata_ref = pbmc_meta,
#'   cluster_col_ref = "classified",
#'   plot_col_ref = "UMAP_1"
#' )
#'
#' plt
#' @export
plot_cols <- function(metadata,
                      cluster_col,
                      cluster_col_called,
                      plot_col,
                      metadata_ref,
                      cluster_col_ref,
                      plot_col_ref) {
  temp1 <- dplyr::group_by_at(metadata, vars(cluster_col, cluster_col_called))
  temp1 <- dplyr::summarise(temp1, med = median(!!dplyr::sym(plot_col), na.rm = TRUE))
  colnames(temp1) <- c("original_cluster", "type", paste(plot_col, "query", sep = "_"))

  temp2 <- dplyr::group_by_at(metadata_ref, cluster_col_ref)
  temp2 <- dplyr::summarise(temp2, med = median(!!dplyr::sym(plot_col_ref), na.rm = TRUE))
  colnames(temp2) <- c("type", paste(plot_col, "ref", sep = "_"))

  temp <- dplyr::left_join(temp1,
    temp2,
    by = "type"
  )
  temp[is.na(temp)] <- 0
  temp[["full"]] <- stringr::str_c(temp[["original_cluster"]], temp[["type"]], sep = "->")
  xmax <- max(temp[, 4])
  ymax <- max(temp[, 3])

  ggplot2::ggplot(
    temp,
    ggplot2::aes_string(
      x = colnames(temp)[4],
      y = colnames(temp)[3],
      label = "full"
    )
  ) +
    geom_point(alpha = 0.23) +
    geom_label(
      alpha = 0.23,
      aes(color = !!dplyr::sym("type")),
      vjust = "inward",
      hjust = "inward"
    ) +
    cowplot::theme_cowplot() +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, xmax * 1.1)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, ymax * 1.1)
    )
}

#' Plot similarity measures on heatmap
#'
#' @param cor_matrix input similarity matrix
#' @param metadata input metadata with per cell tsne or umap cooordinates and cluster ids
#' @param cluster_col colname of clustering data in metadata, defaults to rownames of the
#' metadata if not supplied.
#' @param col color ramp to use
#' @param legend_title legend title to pass to Heatmap
#' @param ... passed to Heatmap
#' @return complexheatmap object
#' @examples
#' res <- clustify(
#'   input = pbmc_matrix_small,
#'   metadata = pbmc_meta,
#'   ref_mat = pbmc_bulk_matrix,
#'   query_genes = pbmc_vargenes,
#'   cluster_col = "classified",
#'   per_cell = FALSE
#' )
#'
#' plt <- plot_cor_heatmap(res)
#'
#' plt
#' @export
plot_cor_heatmap <- function(cor_matrix,
                             metadata = NULL,
                             cluster_col = NULL,
                             col = not_pretty_palette,
                             legend_title = NULL,
                             ...) {
  ComplexHeatmap::Heatmap(cor_matrix,
    col = col,
    heatmap_legend_param = list(title = legend_title),
    ...
  )
}
