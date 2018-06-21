#' Plot a t-SNE colored by feature.
#'
#' @param data input data
#' @param x x variable
#' @param y y variable
#' @param feature feature to color by
#'
#' @export
plot_tsne <- function(data, x = "tSNE_1", y = "tSNE_2", feature) {
  ggplot(data, aes(x = !!sym(x), y = !!sym(y))) +
    geom_point(aes(color = !!enquo(feature)))
}
