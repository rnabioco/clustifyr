#' Heatmap of average expression of most predictive genes across cell types
#'
#' @param important importance measurements of genes (predictors)
#' from the random forest model
#' @param avg_matrix average expression matrix
#' @param meta contains cluster info and classified cell types
#' @param MDA_thresh MeanDecreaseAccuracy threshold
#' @param MDG_thresh MeanDecreaseGini threshold
#'
#' @export
importance_heatmap <- function(importance, avg_matrix, meta, MDA_thresh, MDG_thresh) {

  # Index contains cluster numbers and classified cell types
  index <- unique(meta[(ncol(meta) - 1):ncol(meta)])

  # A list of genes that are most predictive of cell types
  gene_list <- rownames(important[(important$MeanDecreaseAccuracy >= MDA_thresh & important$MeanDecreaseGini >= MDG_thresh), ])

  # Reduce the average expression matrix
  avg_matrix <- avg_matrix[rownames(avg_matrix) %in% gene_list, ]

  # change column names
  for (i in 0:(nrow(index) - 1)) {
    names(avg_matrix)[names(avg_matrix) == i] <- index$classified[index$cluster == i]
  }

  # log transformation of the average expression values
  log_mat <- as.matrix(log(avg_matrix + 1))
  bk1 <- seq(min(log_mat), max(log_mat), length.out = 11)
  colorCount <- ncol(log_mat)
  getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  color <- getPalette(colorCount)

  ComplexHeatmap::Heatmap(
    log_mat,
    col = color,
    name = "heatmap",
    show_heatmap_legend = TRUE
  )
}
