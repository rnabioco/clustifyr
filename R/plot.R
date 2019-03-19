#' Plot a tSNE colored by feature.
#'
#' @param data input data
#' @param x x variable
#' @param y y variable
#' @param feature feature to color by
#' @param legend_name legend name to display, defaults to no name
#' @param c_cols character vector of colors to built color gradient
#' for continuous values. defaults to [`clustifyR::pretty_palette`]
#' @param d_cols character vector of colors for discrete values.
#' defaults to RColorBrewer paired palette
#' @param pt_size point size
#' @param scale_limits defaults to min = 0, max = max(data$x),
#' otherwise a two-element numeric vector indicating min and max to plot
#' @param do.label whether to label each cluster at median center
#' @param do.legend whether to draw legend
#' @export
plot_tsne <- function(data, x = "tSNE_1", y = "tSNE_2",
                      feature,
                      legend_name = "",
                      c_cols = pretty_palette,
                      d_cols = NULL,
                      pt_size = 0.25,
                      scale_limits = NULL,
                      do.label = FALSE,
                      do.legend = TRUE) {
  if ((length(unique(data[[feature]])) > 12) & identical(c_cols, pretty_palette)) {
    print("using ramp colors instead")
    c_cols <- pretty_palette_ramp_d(length(unique(data[[feature]])))
    d_cols <- pretty_palette_ramp_d(length(unique(data[[feature]])))
  }
  # sort data to avoid plotting null values over colors
  data <- dplyr::arrange(data, !!sym(feature))

  p <- ggplot(data, aes_string(x, y)) +
    geom_point(aes_string(color = paste0("`", feature,"`")), # backticks protect special character gene names
      size = pt_size
    )

  if (typeof(data[[feature]]) %in% c(
    "character",
    "logical"
  ) | is.factor(data[[feature]])) {
    if(!is.null(d_cols)) {
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
    if (is.null(scale_limits)){
      scale_limits <- c(ifelse(min(data[[feature]]) < 0,
                             min(data[[feature]]),
                             0),
                      max(data[[feature]]))
    }
    p <- p + scale_color_gradientn(
      colors = c_cols,
      name = legend_name,
      limits = scale_limits
    )
  }

  if (do.label) {
    centers <- dplyr::group_by_(data, .dots = feature)
    centers <- dplyr::summarize(centers, t1 = mean(!!sym(x)), t2 = mean(!!sym(y)))
    p <- p +
      geom_point(data = centers, mapping = aes(x = t1, y = t2), size = 0, alpha = 0) +
      geom_text(data = centers, mapping = aes(x = t1, y = t2, label = centers[[feature]]))
  }

  p <- p + cowplot::theme_cowplot()

  if (do.legend == FALSE) {
    p <- p + theme(legend.position="none")
  }

  p
}

#' Color palette for plotting continous variables
#' @export
pretty_palette <- rev(RColorBrewer::brewer.pal(11, "RdGy")[c(1:5, 7)])

#' Expanded color palette ramp for plotting discrete variables
#' @export
pretty_palette_ramp_d <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))

#' Plot similarity measures on a tSNE
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with per cell tsne coordinates and cluster ids
#' @param ref_data_to_plot colname of data to plot, defaults to all
#' @param cluster_col colname of clustering data in metadata, defaults to rownames of the
#' metadata if not supplied.
#' @param x metadata column name with 1st axis dimension.
#' defaults to "tSNE_1".
#' @param y metadata column name with 2nd axis dimension.
#' defaults to "tSNE_2".
#' @param scale_legends if TRUE scale all legends to maximum values in entire
#' correlation matrix. if FALSE scale legends to maximum for each plot. A
#' two-element numeric vector can also be passed to supply custom values i.e. c(0, 1)
#' @param ... passed to plot_tsne
#'
#' @export
plot_cor <- function(correlation_matrix,
                     metadata,
                     ref_data_to_plot = colnames(correlation_matrix),
                     cluster_col = NULL,
                     x = "tSNE_1",
                     y = "tSNE_2",
                     scale_legends = FALSE,
                     ...) {
  if (!any(ref_data_to_plot %in% colnames(correlation_matrix))) {
    stop("cluster ids not shared between metadata and correlation matrix")
  }

  if(is.null(cluster_col)){
    cluster_col <- "rownames"
    metadata <- tibble::rownames_to_column(metadata, cluster_col)
  }

  cor_df <- as.data.frame(correlation_matrix)
  cor_df <- tibble::rownames_to_column(cor_df, cluster_col)
  cor_df_long <- tidyr::gather(
    cor_df,
    ref_cluster,
    expr,
    -dplyr::matches(cluster_col)
  )

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
  if (typeof(scale_legends) == "logical" && scale_legends){
    scale_limits <- c(ifelse(min(plt_data$expr) < 0,
                               min(plt_data$expr),
                               0),
                        max(max(plt_data$expr)))
  } else if (typeof(scale_legends) == "logical" && !scale_legends){
    scale_limits = NULL
  } else {
    scale_limits = scale_legends
  }

  plts <- vector("list", length(ref_data_to_plot))
  for(i in seq_along(ref_data_to_plot)){
      tmp_data <- dplyr::filter(
        plt_data,
        ref_cluster == ref_data_to_plot[i]
      )
      plts[[i]] <- plot_tsne(tmp_data,
                             x = x,
                             y = y,
                             feature = "expr",
                             legend_name = ref_data_to_plot[i],
                             scale_limits = scale_limits,
                             ...)
  }
  plts
}

#' Plot gene expression on to tSNE
#'
#'
#' @param expr_mat input single cell matrix
#' @param metadata data.frame with tSNE coordinates
#' @param genes gene(s) to color tSNE
#' @param cell_col column name in metadata containing cell ids, defaults
#' to rownames if not supplied
#' @param ... additional arguments passed to `[clustifyR::plot_tsne()]`
#' @export
plot_gene <- function(expr_mat,
                      metadata,
                      genes,
                      cell_col = NULL,
                      ...) {

  genes_to_plot <- genes[genes %in% rownames(expr_mat)]
  genes_missing <- setdiff(genes, genes_to_plot)

  if (length(genes_missing) != 0) {
    warning(paste0("the following genes were not present in the input matrix ",
            paste(genes_missing, collapse = ",")))
  }

  if (length(genes_to_plot) == 0) {
    stop("no genes present to plot")
  }
  expr_dat <- t(as.matrix(expr_mat[genes_to_plot, , drop = FALSE]))
  expr_dat <- tibble::rownames_to_column(as.data.frame(expr_dat), "cell")

  if(is.null(cell_col)){
    mdata <- tibble::rownames_to_column(metadata, "cell")
    cell_col <- "cell"
  } else {
    mdata <- metadata
  }

  if(!cell_col %in% colnames(mdata)) {
    stop("please supply a cell_col that is present in metadata")
  }

  plt_dat <- dplyr::left_join(expr_dat, mdata,
                       by = c("cell" = cell_col))

  lapply(genes_to_plot,
         function(gene){
           plot_tsne(plt_dat,
                     feature = gene,
                     legend_name = gene,
                     ...)
         })
}

#' Plot called clusters on a tSNE, for each reference cluster given
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne coordinates and cluster ids
#' @param ref_data_to_plot colname of data to plot, defaults to all
#' @param ... passed to plot_tsne
#'
#' @export
plot_call <- function(correlation_matrix,
                      metadata,
                      ref_data_to_plot = colnames(correlation_matrix),
                      ...) {
  df_temp <- as.data.frame(t(apply(correlation_matrix, 1, function(x) x - max(x))))
  df_temp[df_temp == 0] <- "1"
  df_temp[df_temp != "1"] <- "0"
  plot_cor(df_temp,
           metadata,
           ref_data_to_plot,
           ...)
}

#' Plot best calls for each cluster on a tSNE
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne coordinates and cluster ids
#' @param col metadata column, can be cluster or cellid
#' @param collapse_to_cluster if a column name is provided, takes the most frequent call of entire cluster to color in plot
#' @param threshold minimum correlation coefficent cutoff for calling clusters
#' @param x x variable
#' @param y y variable
#' @param plot_r whether to include second plot of cor eff for best call
#' @param ... passed to plot_tsne
#'
#' @export
plot_best_call <- function(correlation_matrix,
                           metadata,
                           col = "cluster",
                           collapse_to_cluster = FALSE,
                           threshold = 0,
                           x = "tSNE_1", y = "tSNE_2",
                           plot_r = FALSE,
                           ...) {
  col_meta <- colnames(metadata)
  if("type" %in% col_meta | "type2" %in% col_meta){
    warning('metadata column name clash of "type"/"type2"')
    return()
  }
  df_temp <- tibble::as_tibble(correlation_matrix, rownames = col)
  df_temp <- tidyr::gather(df_temp, key = type, value = r, -!!col)
  df_temp[["type"]][df_temp$r < threshold] <- paste0("r<", threshold,", unassigned")
  df_temp <- dplyr::top_n(dplyr::group_by_at(df_temp, 1), 1, r)
  if (nrow(df_temp) != nrow(correlation_matrix)) {
    clash <- dplyr::group_by_at(df_temp, 1)
    clash <- dplyr::summarize(clash, n = n())
    clash <- dplyr::filter(clash, n > 1)
    clash <- dplyr::pull(clash, 1)
    df_temp[lapply(df_temp[,1], FUN = function(x) x %in% clash)[[1]],2] <- paste0(df_temp[["type"]][lapply(df_temp[,1], FUN = function(x) x %in% clash)[[1]]], "-CLASH!")
    df_temp <- dplyr::distinct(df_temp, exclude = "type", .keep_all =T)
  }
  df_temp_full <- dplyr::left_join(metadata, df_temp, by = col)

  if(collapse_to_cluster != FALSE){
    df_temp_full2 <- dplyr::mutate(df_temp_full, type2 = metadata[[collapse_to_cluster]])
    df_temp_full2 <- dplyr::group_by(df_temp_full2, type, type2)
    df_temp_full2 <- dplyr::summarize(df_temp_full2, sum = sum(r), n = n())
    df_temp_full2 <- dplyr::group_by(df_temp_full2, type2)
    df_temp_full2 <- dplyr::arrange(df_temp_full2, desc(n), desc(sum))
    df_temp_full2 <- dplyr::filter(df_temp_full2, type != paste0("r<", threshold,", unassigned"))
    df_temp_full2 <- dplyr::slice(df_temp_full2, 1)
    df_temp_full2 <- dplyr::right_join(df_temp_full2, select(df_temp_full, -type), by = stats::setNames(collapse_to_cluster, "type2"))
    df_temp_full <- dplyr::mutate(df_temp_full2, type = tidyr::replace_na(type, paste0("r<", threshold,", unassigned")))
  }

  g <- plot_tsne(df_temp_full,
            feature = "type",
            x = x, y = y,
            ...)

  if (plot_r == T) {
    l <- list()
    l[[1]] <- g
    l[[2]] <- plot_tsne(df_temp_full,
              feature = "r",
              x = x, y = y,
              ...)
  } else {
    l <- g
  }

  l
}
