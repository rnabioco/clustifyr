
#' Refine clustering using reference data
#'
#'
#' @param expr_mat input single cell data after gene subsetting
#' @param bulk_mat input bulk data after gene subsetting
#' @param sc_cluster initial clustering assignment.
#' Note that the number of bulk cell types must be no more than
#' the number of clusters. for additional clusters,
#' there is an assignment of default_similiarity
#' @param lambda control parameters
#' @param epsilon control parameters
#' @param if_compute_sigma whether compute the standard deviation
#'  of each gene's expression within each cluster.
#' If not, the function assumes that all standard deviations equal to 1.
#' @param disagreement_freq maximal fraction of changes in cluster
#'  assignment after one iteration to be considered as converged.
#' @param num_iteration maximal number of iteration before exiting
#' the algorithm
#' @param default_similiarity default similarity score if there is
#'  no reference data to compare to
#' @param if_plot whether to plot the clustering result at each iteration.
#' @param tsne_coord t-SNE coordinate of every cells for plotting
#' @param cluster_name Name of each cluster for plotting
#' @param compute_method corr_coef
#' @param ... control parameters for calculating similarity score
#' @export
refine_clusters <- function(expr_mat, sc_cluster, bulk_mat, lambda=0.5, epsilon=1,
                            if_compute_sigma=TRUE, num_iteration=100, default_similiarity=-5,
                            disagreement_freq=0.001, if_plot=FALSE, tsne_coord=NULL,
                            cluster_names=NULL, compute_method, ...) {
  # step 1. check input conditions
  num_cells <- ncol(expr_mat); num_cluster <- max(sc_cluster); num_bulk <- ncol(bulk_mat);
  if (num_cluster < num_bulk) {
    error('refine_clusters: there are more bulk samples than clusters');
  }

  # step 2. iteratively re-assign the clusters
  curr_cluster <- sc_cluster; # a numeric value for each cell, corr. to the column of the bulk expr.
  if (if_plot)
    print(plot_tsne(data.frame(tSNE_1=tsne_coord[,'tSNE_1'], tSNE_2=tsne_coord[,'tSNE_2'], Classification=cluster_names[curr_cluster]), feature="Classification", legend_name='Initial'));
  curr_num_iteration <- 0;
  bulk_info <- compute_bulk_score(expr_mat, bulk_mat, num_cluster, default_similiarity, compute_method, ...); # the bulk info has nothing to do with current clusters
  while (1) {
    curr_num_iteration <- curr_num_iteration + 1;
    print(paste('Iteration: ', curr_num_iteration, sep=''));
    # compute centroid info
    centroid_info <- compute_centroid(expr_mat, curr_cluster);
    if (!if_compute_sigma) {
      centroid_info$sigma[,] <- 1;
    }

    # compute log probability (up to a constant) for each single cell and each cluster
    log_score <- lambda*epsilon*bulk_info; # similarity prob between each single cell and each bulk sample
    for (j in 1:num_bulk) { # for each cluster
      id <- which(centroid_info$sigma[, j]>0);
      log_score[,j] <- log_score[,j] + sapply(1:num_cells, function(i) -(1-lambda)/2*mean((expr_mat[id,i]-centroid_info$mu[id,j])^2/centroid_info$sigma[id,j]^2) - (1-lambda)*mean(log(centroid_info$sigma[id,j])))
    }

    # compute the new assignment
    new_cluster <- apply(log_score, 1, which.max);
    if (length(unique(new_cluster)) < num_cluster) {
      warning(paste('refine_clusters: iteration ', curr_num_iteration, ' some clusters are degenerated.', sep=''));
    }
    if (mean(curr_cluster!=new_cluster) <= disagreement_freq) {
      message(paste("refine_clusters: converged after ", curr_num_iteration, " iterations.", sep=''));
      return(new_cluster);
    }
    if (curr_num_iteration >= num_iteration) {
      warning("refine_clusters: reach maximal number of iterations. Clustering fails to converge.");
      return(new_cluster);
    }
    curr_cluster <- new_cluster;
    if (if_plot)
      print(plot_tsne(data.frame(tSNE_1=tsne_coord[,'tSNE_1'], tSNE_2=tsne_coord[,'tSNE_2'], Classification=cluster_names[new_cluster]), feature="Classification", legend_name=paste('Iteration ', curr_num_iteration, sep='')));
  }
}

#' @noRd
compute_centroid <- function(expr_mat, sc_cluster) {
  num_cluster <- max(sc_cluster); num_genes <- nrow(expr_mat);
  mu_info <- matrix(NA, nrow=num_genes, ncol=num_cluster);
  sigma_info <- matrix(NA, nrow=num_genes, ncol=num_cluster);
  for (i in 1:num_cluster) {
    curr_cluster_expr <- expr_mat[,sc_cluster==i];
    mu_info[,i] <- Matrix::rowMeans(curr_cluster_expr);
    sigma_info[,i] <- sqrt(Matrix::rowMeans(curr_cluster_expr^2) - mu_info[,i]^2);
  }
  return(list(mu=mu_info, sigma=sigma_info))
}

#' @noRd
compute_bulk_score <- function(expr_mat, bulk_mat, num_cluster, default_similiarity, compute_method, ...) {
  # for clusters corresponding to bulk data: compute similarity score
  num_cells <- ncol(expr_mat); num_bulk <- ncol(bulk_mat);
  bulk_score <- matrix(NA, nrow=num_cells, ncol=num_bulk);
  for (i in 1:num_cells) {
    curr_cell_expr <- expr_mat[,i];
    bulk_score[i,] <- sapply(1:num_bulk, function(x) compute_similarity(curr_cell_expr, bulk_mat[,x], compute_method, ...))
  }
  # for clusters without bulk data, use a default value
  if (num_cluster > num_bulk) {
    for (j in (num_bulk+1):num_cluster) {
      bulk_score <- cbind(bulk_score, rep(default_similiarity, num_cells));
    }
  }
  return(bulk_score);
}


