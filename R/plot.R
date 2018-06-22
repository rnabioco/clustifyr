#' Plot a t-SNE colored by feature.
#'
#' @param data input data
#' @param x x variable
#' @param y y variable
#' @param feature feature to color by
#'
#' @export
plot_tsne <- function(data, x = "tSNE_1", y = "tSNE_2",
                      feature, legend_name = "",
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

#' @noRd
pretty_palette <- rev(RColorBrewer::brewer.pal(11, "RdGy")[c(1:5, 7)])

#' Plot similarity measures on a tSNE
#'
#' @param correlation_matrix input similarity matrix
#' @param meta_data input metadata with tsne coordinates and cluster ids
#' @param bulk_data_to_plot colname of data to plot, defaults to all
#' @param ... passed to plot_tsne
#'
#' @export
plot_cor <- function(correlation_matrix,
                     meta_data,
                     bulk_data_to_plot = colnames(correlation_matrix)) {

  if (!any(bulk_data_to_plot %in% colnames(correlation_matrix))){
    stop("cluster ids not shared between meta_data and correlation matrix")
  }

  cor_df <- as.data.frame(correlation_matrix)
  cor_df <- tibble::rownames_to_column(cor_df, "cluster")
  cor_df_long <- tidyr::gather(cor_df, bulk_cluster, expr, -cluster)

  # checks matrix rownames, 2 branches for cluster number (avg) or cell bar code (each cell)
  if(cor_df$cluster[1] %in% meta_data$cluster){
    plt_data <- dplyr::left_join(cor_df_long,
                                 meta_data,
                                 by = "cluster")
  } else {
    plt_data <- dplyr::left_join(cor_df_long,
                                 meta_data,
                                 by = c("cluster" = "rn"))
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
#' @param meta_data input metadata with tsne coordinates and cluster ids
#' @param bulk_data_to_plot colname of data to plot, defaults to all
#' @param ... passed to plot_tsne
#'
#' @export
plot_call <- function(correlation_matrix,
                      meta_data,
                      bulk_data_to_plot = colnames(correlation_matrix)) {
  df_temp <- as.data.frame(t(apply(correlation_matrix, 1, function(x) x - max(x))))
  df_temp[df_temp==0]="1"
  df_temp[df_temp!="1"]="0"
  plot_cor(df_temp, meta_data, bulk_data_to_plot)
}
