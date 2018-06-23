# evaluate whether refine clusters works
#data("pbmc4k_meta"); data("pbmc4k_matrix"); data("pbmc_bulk_matrix"); data("pbmc4k_markers_M3Drop")
#
## construct sc_cluster
#gene_constraints <- list(rownames(pbmc4k_matrix), rownames(pbmc_bulk_matrix), pbmc4k_markers_M3Drop[,'Gene']);
#sc_expr <- select_gene_subset(pbmc4k_matrix, gene_constraints);
#bulk_expr <- select_gene_subset(pbmc_bulk_matrix, gene_constraints);
#sc_tsne_coord <- pbmc4k_meta[,c('tSNE_1', 'tSNE_2')];
#sc_cluster <- rep(0, ncol(sc_expr));
#sc_cluster[which(pbmc4k_meta[,'classified']=="CD14+ Monocytes")] <- 1; # 1 monocytes
#sc_cluster[which(pbmc4k_meta[,'classified']=="FCGR3A+ Monocytes")] <- 1;
#sc_cluster[which(pbmc4k_meta[,'classified']=="B cells")] <- 2; # 2 B cells
#sc_cluster[which(pbmc4k_meta[,'classified']=="CD4 T cells, 1")] <- 3; # 3 T cells
#sc_cluster[which(pbmc4k_meta[,'classified']=="CD4 T cells, 2")] <- 3;
#sc_cluster[which(pbmc4k_meta[,'classified']=="CD8 T cells")] <- 3;
#sc_cluster[which(pbmc4k_meta[,'classified']=="Dendritic cells, 1")] <- 4; # 4 DC
#sc_cluster[which(pbmc4k_meta[,'classified']=="Dendritic Cells, 2?")] <- 4;
#sc_cluster[which(pbmc4k_meta[,'classified']=="NK cells")] <- 5; # 5 NK cells
#sc_cluster[which(pbmc4k_meta[,'classified']=="Megakaryocytes")] <- 6; # 6 megakaryocytes (not found)
#bulk_expr <- bulk_expr[, c(3, 1, 7, 2, 5)]
#colnames(bulk_expr) <- c("Monocytes", "B Cells", "T Cells", "Dendritic Cells", "NK Cells");
#cluster_names <- c("Monocytes", "B Cells", "T Cells", "Dendritic Cells", "NK Cells", "Megakaryocytes");
#
## re-assign the cluster information
#res <- refine_clusters(sc_expr, sc_cluster, bulk_expr, lambda=0, epsilon=10, if_compute_sigma=TRUE, num_iteration=5, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sc_cluster, bulk_expr, lambda=0.5, epsilon=10, if_compute_sigma=TRUE, num_iteration=5, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sc_cluster, bulk_expr, lambda=1, epsilon=10, if_compute_sigma=TRUE, num_iteration=5, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sample(sc_cluster, length(sc_cluster)), bulk_expr, lambda=0, epsilon=10, if_compute_sigma=TRUE, num_iteration=10, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sample(sc_cluster, length(sc_cluster)), bulk_expr, lambda=0.5, epsilon=10, if_compute_sigma=TRUE, num_iteration=10, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sample(sc_cluster, length(sc_cluster)), bulk_expr, lambda=1, epsilon=10, if_compute_sigma=TRUE, num_iteration=10, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#
#res <- refine_clusters(sc_expr, sc_cluster, bulk_expr, lambda=0, epsilon=10, if_compute_sigma=FALSE, num_iteration=5, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sc_cluster, bulk_expr, lambda=0.5, epsilon=10, if_compute_sigma=FALSE, num_iteration=5, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sc_cluster, bulk_expr, lambda=1, epsilon=10, if_compute_sigma=FALSE, num_iteration=5, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sample(sc_cluster, length(sc_cluster)), bulk_expr, lambda=0, epsilon=10, if_compute_sigma=FALSE, num_iteration=10, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sample(sc_cluster, length(sc_cluster)), bulk_expr, lambda=0.5, epsilon=10, if_compute_sigma=FALSE, num_iteration=10, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#res <- refine_clusters(sc_expr, sample(sc_cluster, length(sc_cluster)), bulk_expr, lambda=1, epsilon=10, if_compute_sigma=FALSE, num_iteration=10, default_similiarity=-1, disagreement_freq=0.0001, if_plot=TRUE, tsne_coord=sc_tsne_coord, cluster_names=cluster_names, compute_method=corr_coef);
#
