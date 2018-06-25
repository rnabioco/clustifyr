#' Plot a tSNE colored by feature.
#'
#' @param data input data
#' @param x x variable
#' @param y y variable
#' @param feature feature to color by
#' @param legend_name legend name to display, defaults to no name
#' @param cols character vector of colors to built color gradient
#' for continuous values. defaults to [`clustifyR::pretty_palette`]
#' @param pt_size point size
#' @export
plot_tsne <- function(data, x = "tSNE_1", y = "tSNE_2",
                      feature,
                      legend_name = "",
                      cols = pretty_palette,
                      pt_size = 0.25) {

  p <- ggplot(data, aes_string(x, y)) +
    geom_point(aes_string(color = feature),
               size = pt_size)

  if(typeof(data[[feature]]) %in% c("character",
                                    "logical") | is.factor(data[[feature]])){
    p <- p + scale_color_brewer(palette = "Paired",
                                name = legend_name)

  } else {
    p <- p + scale_color_gradientn(colors = cols,
                                   name = legend_name)
  }

  p + cowplot::theme_cowplot()
}

#' Color palette for plotting continous variables
#' @export
pretty_palette <- rev(RColorBrewer::brewer.pal(11, "RdGy")[c(1:5, 7)])

#' Plot similarity measures on a tSNE
#'
#' @param correlation_matrix input similarity matrix
#' @param metadata input metadata with tsne coordinates and cluster ids
#' @param bulk_data_to_plot colname of data to plot, defaults to all
#' @param cluster_col colname of cluster IDs, defaults to cluster
#' @param ... passed to plot_tsne
#'
#' @export
plot_cor <- function(correlation_matrix,
                     metadata,
                     bulk_data_to_plot = colnames(correlation_matrix),
                     cluster_col = "cluster") {

  if (!any(bulk_data_to_plot %in% colnames(correlation_matrix))){
    stop("cluster ids not shared between metadata and correlation matrix")
  }

  cor_df <- as.data.frame(correlation_matrix)
  cor_df <- tibble::rownames_to_column(cor_df, cluster_col)
  cor_df_long <- tidyr::gather(cor_df, bulk_cluster,
                               expr, -dplyr::matches(cluster_col))

  # checks matrix rownames, 2 branches for cluster number (avg) or cell bar code (each cell)
  if(cor_df[[cluster_col]][1] %in% metadata[[cluster_col]]){
    plt_data <- dplyr::left_join(cor_df_long,
                                 metadata,
                                 by = cluster_col)
  } else {
    plt_data <- dplyr::left_join(cor_df_long,
                                 metadata,
                                 by = structure(names = cluster_col, "rn"))
  }

  lapply(bulk_data_to_plot,
         function(x){
           tmp_data <- dplyr::filter(plt_data,
                                     bulk_cluster == x)
           plot_tsne(tmp_data, feature = "expr", legend_name = x)
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
                      bulk_data_to_plot = colnames(correlation_matrix)) {
  df_temp <- as.data.frame(t(apply(correlation_matrix, 1, function(x) x - max(x))))
  df_temp[df_temp==0]="1"
  df_temp[df_temp!="1"]="0"
  plot_cor(df_temp, metadata, bulk_data_to_plot)
}
