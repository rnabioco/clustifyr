#' Heatmap of average expression of most predictive genes across cell types
#'
#' @param important importance measurements of genes (predictors) from the random forest model
#' @param avg_matrix average expression matrix
#' @param meta contains cluster info and classified cell types
#' @param MDA_thresh MeanDecreaseAccuracy threshold
#' @param MDG_thresh MeanDecreaseGini threshold
#'
#' @export

important_heatmap <- function(important, avg_matrix, meta, MDA_thresh, MDG_thresh){
  # Index contains cluster numbers and classified cell types
  index <-unique(meta[(ncol(meta)-1):ncol(meta)])
  # A list of genes that are most predictive of cell types
  gene_list <- rownames(important[(important$MeanDecreaseAccuracy >= MDA_thresh & important$MeanDecreaseGini >= MDG_thresh),])
  # Reduce the average expression matrix
  avg_matrix <- avg_matrix[rownames(avg_matrix) %in% gene_list, ]
  # change column names
  for(i in 0:(nrow(index)-1)) {
    names(avg_matrix)[names(avg_matrix) == i] <- index$classified[index$cluster == i]
  }
  # log transformation of the average expression values
  log_mat  <- as.matrix(log(avg_matrix+1))
  # rc <- rainbow(nrow(log_mat))
  # cc <- rainbow(ncol(log_mat))
  # colorside = gray(1:10/10)
  bk1 = seq(min(log_mat),max(log_mat),length.out=11)
  color <- colorpanel(ncol(log_mat),"azure","steelblue1","royalblue4") 
  # Plot the heatmap
  heatmap<-heatmap.2( log_mat, 
             key=TRUE, 
             key.title = NA,
             key.xlab = 'log(avg_expression)',
             trace="none",
             # ColSideColors=colorside,
             cexRow=1,
             cexCol=1,
             breaks=bk1,
             col=color,
             sepwidth=c(0.1,0.1),
             sepcolor="white",
             colsep=1:ncol(log_mat),
             rowsep=1:nrow(log_mat),
             margins=c(10,6))
}

