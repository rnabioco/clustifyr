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
#' @export
plot_tsne <- function(data, x = "tSNE_1", y = "tSNE_2",
                      feature,
                      legend_name = "",
                      c_cols = pretty_palette,
                      d_cols = NULL,
                      pt_size = 0.25,
                      scale_limits = NULL) {

  # sort data to avoid plotting null values over colors
  data <- arrange(data, !!sym(feature))

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

  p + cowplot::theme_cowplot()
}

#' Color palette for plotting continous variables
#' @export
pretty_palette <- rev(RColorBrewer::brewer.pal(11, "RdGy")[c(1:5, 7)])

#' Plot similarity measures on a tSNE
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with per cell tsne coordinates and cluster ids
#' @param bulk_data_to_plot colname of data to plot, defaults to all
#' @param cluster_col colname of clustering data in metadata, defaults to rownames of the
#' metadata if not supplied.
#' @param dim1_col metadata column name with 1st axis dimension.
#' defaults to "tSNE_1".
#' @param dim2_col metadata column name with 2nd axis dimension.
#' defaults to "tSNE_2".
#' @param scale_legends if TRUE scale all legends to maximum values in entire
#' correlation matrix. if FALSE scale legends to maximum for each plot. A
#' two-element numeric vector can also be passed to supply custom values i.e. c(0, 1)
#' @param ... passed to plot_tsne
#'
#' @export
plot_cor <- function(correlation_matrix,
                     metadata,
                     bulk_data_to_plot = colnames(correlation_matrix),
                     cluster_col = NULL,
                     dim1_col = "tSNE_1",
                     dim2_col = "tSNE_2",
                     scale_legends = FALSE,
                     ...) {
  if (!any(bulk_data_to_plot %in% colnames(correlation_matrix))) {
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
    bulk_cluster,
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

  plts <- vector("list", length(bulk_data_to_plot))
  for(i in seq_along(bulk_data_to_plot)){
      tmp_data <- dplyr::filter(
        plt_data,
        bulk_cluster == bulk_data_to_plot[i]
      )
      plts[[i]] <- plot_tsne(tmp_data,
                             x = dim1_col,
                             y = dim2_col,
                             feature = "expr",
                             legend_name = bulk_data_to_plot[i],
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
  genes_missing <- setdiff(genes_to_plot, genes)

  if (length(genes_missing) != 0) {
    warning(paste0("the following genes were not present in the input matrix ",
            paste(genes_missing, collapse = ",")))
  }

  if (length(genes_to_plot) == 0) {
    stop("no genes present to plot")
  }
  expr_dat <- t(as.matrix(expr_mat[genes_to_plot, ]))
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

  plt_dat <- left_join(expr_dat, metadata,
                       by = c("cell" = cell_col))

  lapply(genes,
         function(gene){
           plot_tsne(plt_dat,
                     feature = gene,
                     legend_name = gene,
                     ...)
         })
}

#' Plot called clusters on a tSNE
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne coordinates and cluster ids
#' @param bulk_data_to_plot colname of data to plot, defaults to all
#' @param ... passed to plot_tsne
#'
#' @export
plot_call <- function(correlation_matrix,
                      metadata,
                      bulk_data_to_plot = colnames(correlation_matrix),
                      ...) {
  df_temp <- as.data.frame(t(apply(correlation_matrix, 1, function(x) x - max(x))))
  df_temp[df_temp == 0] <- "1"
  df_temp[df_temp != "1"] <- "0"
  plot_cor(df_temp,
           metadata,
           bulk_data_to_plot,
           ...)
}

#' Plot called clusters on a tSNE
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne coordinates and cluster ids
#' @param feature feature name, defaults to "type"
#' @param ... passed to plot_tsne
#'
#' @export
plot_best_call <- function(correlation_matrix,
                           metadata,
                           feature = "type",
                           ...) {
  df_temp <- tibble::as_tibble(correlation_matrix, rownames = "cluster")
  df_temp <- tidyr::gather(df_temp, key = !!feature, value = r, -cluster)
  df_temp <- dplyr::top_n(dplyr::group_by(df_temp, cluster), 1, r)
  df_temp_full <- left_join(metadata, df_temp, by = "cluster")

  plot_tsne(df_temp_full,
            feature = feature,
            ...)
}

plot_best_call <- function(correlation_matrix,
                           metadata,
                           feature = "type") {
  df_temp <- tibble::as_tibble(correlation_matrix, rownames = "cluster")
  df_temp <- tidyr::gather(df_temp, key = !!feature, value = r, -cluster)
  df_temp <- dplyr::top_n(dplyr::group_by(df_temp, cluster), 1, r)
  df_temp_full <- left_join(metadata, df_temp, by = "cluster")

  plot_tsne(df_temp_full, feature = feature)
}
