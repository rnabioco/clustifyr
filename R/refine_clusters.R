# math model
#   P(x belongs to cluster i)
# = P(x is similar to bulk i)^\lambda  # must similar to a pre-defined cell type
#   * P(x is close to centroid of cluster i)^(1-\lambda)  # must make sure points within a cluster are "closed" to each other

# Algorithm description
# 1. Initiate a cluster using k-means, with k = known cell-types = # bulk data
# 2. For each cluster i and each cell (expr vector x), compute P(x, i). Assign x to the cluster which maximizes P(x,i)
# 3. Create a new cluster and update the centroid of each cluster
# 4. Repeat Step 2-3 until convergence (!!need to prove)

# implementation issues
# 1. log-scale to avoid issue of machine precision
# 2. P(x is similar to bulk i): make use of the similarity score (\in [-1, 1]).
#    Use exponential scoring function P(similarity) = C exp(\epsilon * similarity_score), where C=\epsilon/(exp(\epsilon) - exp(-\epsilon)).
#TODO: Consider changing similarity score to [0, 1] to simplify the expression of C. Probably can ignore because we are comparing relative probabilities
# 3. P(x is closed to centroid of cluster i): make use of Euclidean distance between x and centroid \mu_i
#    Use Euclidean distance because of better change to have convergence (k-means using Euclidean distance is guaranteed to converge to the unique point)
#    Assume something like random walk (i.e. diffusion) -> Gaussian distribution
#    P = 1/sqrt(2*\pi*\sigma^2) exp(-|x-\mu_i|^2/2/\sigma^2)
#TODO: In principle it's a multivariate Gaussian distribution, so \sigma^2 should be replaced by the covariance matrix \Sigma
#    If we assume that all cells/centroids have the same set of genes (which should be!), and every gene has the same amount of error (after normalizing genes!), this expression should be ok.
#    Otherwise, consider model \Sigma explicitly. As a first-order approximation, can assume \Sigma is a diagonal matrix,
#    where each diagonal element represents the stochasticity of each gene (which can be estimate from looking at cluster-wise/data-wise variance)
#    In that case, P \propto exp(-0.5*\sum_j (x_j-\mu_{i,j})^2/\sigma_j^2)/\prod_j \sigma_j, where j iterates over every gene. Note that \sigma_j should at least constant per cluster, in which case \sigma_j -> \sigma_{i,j}

# Simplification using logarithm scale
# ln P(x similar to bulk i)^\lambda = \lambda \epsilon *similarity_score + \lambda ln C   (EQ 1)
# ln P(x closed to centroid i)^(1-\lambda) = -(1-\lambda) 0.5*\sum_j (x_j-\mu_{i,j})^2/\sigma_{i,j}^2 - (1-\lambda)\sum_j ln \sigma_{i,j} - 0.5*(1-\lambda)*n ln(2*\pi)  (EQ 2)
# overall probability
# ln P = (EQ 1) + (EQ 2) = \lambda \epsilon *similarity_score - (1-\lambda) 0.5*\sum_j (x_j-\mu_{i,j})^2/\sigma_{i,j}^2 - (1-\lambda)\sum_j ln \sigma_{i,j} + some constant
# similarity score can be computed by compute_similarity function
# \mu_{i,j} and \sigma_{i,j} can be computed per cluster

# Numerical consideration
# What happened if a cluster runs out of data points? Is it possible that everything will be converge to one cluster, possibly due to any normalization issue? Need a constant penalty factor if that is the case (i.e. P -> P/(P+some constant))

# Testing consideration
# Can we have authentic data? Running permutation definitely works.

# refine_clusters: refine clustering using data
# sc_expr, bulk_expr: data after gene subsetting
# sc_cluster: initial clustering assignment.
# the number of bulk cell types must be no more than the number of clusters. for additional clusters, there is an assignment of default_similiarity
# lambda, epsilon: control parameters
# num_iteration: maximal number of iteration before exiting the algorithm
# compute_method, ...: control parameters for calculating similarity score
refine_clusters <- function(sc_expr, sc_cluster, bulk_expr, lambda=0.5, epsilon=1, if_compute_sigma=TRUE, num_iteration=100, default_similiarity=-5, disagreement_freq=0.001, if_plot=FALSE, tsne_coord=NULL, cluster_names=NULL, compute_method, ...) {
  # step 1. check input conditions
  num_cells <- ncol(sc_expr); num_cluster <- max(sc_cluster); num_bulk <- ncol(bulk_expr);
  if (num_cluster < num_bulk) {
    error('refine_clusters: there are more bulk samples than clusters');
  }

  # step 2. iteratively re-assign the clusters
  curr_cluster <- sc_cluster; # a numeric value for each cell, corr. to the column of the bulk expr.
  print(plot_tsne(data.frame(tSNE_1=tsne_coord[,'tSNE_1'], tSNE_2=tsne_coord[,'tSNE_2'], Classification=cluster_names[curr_cluster]), feature="Classification", legend_name='Initial'));
  curr_num_iteration <- 0;
  while (1) {
    curr_num_iteration <- curr_num_iteration + 1;
    print(paste('Iteration: ', curr_num_iteration, sep=''));
    # compute centroid info
    centroid_info <- compute_centroid(sc_expr, curr_cluster);
    if (!if_compute_sigma) {
      centroid_info$sigma[,] <- 1;
    }
    bulk_info <- compute_bulk_score(sc_expr, bulk_expr, num_cluster, default_similiarity, compute_method, ...);

    # compute log probability (up to a constant) for each single cell and each cluster
    log_score <- lambda*epsilon*bulk_info; # similarity prob between each single cell and each bulk sample
    for (j in 1:num_bulk) { # for each cluster
      id <- which(centroid_info$sigma[, j]>0);
      log_score[,j] <- log_score[,j] + sapply(1:num_cells, function(i) -(1-lambda)/2*mean((sc_expr[id,i]-centroid_info$mu[id,j])^2/centroid_info$sigma[id,j]^2)/length(id) - (1-lambda)*mean(log(centroid_info$sigma[id,j]))/length(id))
    }

    # compute the new assignment
    new_cluster <- apply(log_score, 1, which.max);
    if (length(unique(new_cluster)) < num_cluster) {
      warning(paste('refine_clusters: iteration ', curr_num_iteration, ' some clusters are degenerated.', sep=''));
    }
    if (mean(curr_cluster!=new_cluster) <= disagreement_freq) {
      message(paste("refine_clusters: converged after ", curr_num_iteration, " runs.", sep=''));
      return(new_cluster);
    }
    if (curr_num_iteration >= num_iteration) {
      warning("refine_clusters: reach maximal number of iterations. Clustering fails to converge.");
      return(new_cluster);
    }
    curr_cluster <- new_cluster;
    print(plot_tsne(data.frame(tSNE_1=tsne_coord[,'tSNE_1'], tSNE_2=tsne_coord[,'tSNE_2'], Classification=cluster_names[new_cluster]), feature="Classification", legend_name=paste('Iteration ', curr_num_iteration, sep='')));
  }
}

compute_centroid <- function(sc_expr, sc_cluster) {
  num_cluster <- max(sc_cluster); num_genes <- nrow(sc_expr);
  mu_info <- matrix(NA, nrow=num_genes, ncol=num_cluster);
  sigma_info <- matrix(NA, nrow=num_genes, ncol=num_cluster);
  for (i in 1:num_cluster) {
    curr_cluster_expr <- sc_expr[,sc_cluster==i];
    mu_info[,i] <- Matrix::rowMeans(curr_cluster_expr);
    sigma_info[,i] <- sqrt(Matrix::rowMeans(curr_cluster_expr^2) - mu_info[,i]^2);
  }
  return(list(mu=mu_info, sigma=sigma_info))
}

compute_bulk_score <- function(sc_expr, bulk_expr, num_cluster, default_similiarity, compute_method, ...) {
  # for clusters corresponding to bulk data: compute similarity score
  num_cells <- ncol(sc_expr); num_bulk <- ncol(bulk_expr);
  bulk_score <- matrix(NA, nrow=num_cells, ncol=num_bulk);
  for (i in 1:num_cells) {
    curr_cell_expr <- sc_expr[,i];
    bulk_score[i,] <- sapply(1:num_bulk, function(x) compute_similarity(curr_cell_expr, bulk_expr[,x], compute_method, ...))
  }
  # for clusters without bulk data, use a default value
  if (num_cluster > num_bulk) {
    for (j in (num_bulk+1):num_cluster) {
      bulk_score <- cbind(bulk_score, rep(default_similiarity, num_cells));
    }
  }
  return(bulk_score);
}


