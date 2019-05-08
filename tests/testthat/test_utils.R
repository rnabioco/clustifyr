context("utils")

test_that("get_vargenes works for both matrix and dataframe form", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  var1 <- get_vargenes(pbmc4k_mm)
  var2 <- get_vargenes(pbmc4k_markers)

  expect_equal(var1[1], var2[1])
})

test_that("matrixize_markers with remove_rp option", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  pbmc4k_mm2 <- matrixize_markers(pbmc4k_markers,
    remove_rp = T
  )

  expect_true(nrow(pbmc4k_mm) != nrow(pbmc4k_mm2))
})

test_that("matrixize_markers to turn matrix into ranked list", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers, n = 50)
  pbmc4k_mm2 <- matrixize_markers(pbmc4k_mm, ranked = T, unique = T)

  expect_true(nrow(pbmc4k_mm) < nrow(pbmc4k_mm2))
})

test_that("matrixize_markers uses supplied labels", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers, n = 50, labels = pbmc4k_meta)
  pbmc4k_mm2 <- matrixize_markers(pbmc4k_mm, labels = unique(pbmc4k_meta$classified), ranked = T)

  expect_true(nrow(pbmc4k_mm) < nrow(pbmc4k_mm2))
})

test_that("average_clusters works as intended", {
  pbmc4k_avg2 <- average_clusters(pbmc4k_matrix,
    pbmc4k_meta,
    log_scale = F
  )
  expect_equal(nrow(pbmc4k_avg2), 2663)
})

test_that("average_clusters detects wrong cluster ident", {
  expect_error(pbmc4k_avg2 <- average_clusters(pbmc4k_matrix,
                                               matrix(5,5),
                                               log_scale = F
  ))
})

test_that("average_clusters able to coerce factors", {
  col <- factor(pbmc4k_meta$cluster)
  pbmc4k_avg2 <- average_clusters(pbmc4k_matrix,
    col,
    log_scale = F
  )
  expect_equal(nrow(pbmc4k_avg2), 2663)
})

test_that("average_clusters works with median option", {
  pbmc4k_avg2 <- average_clusters(pbmc4k_matrix,
    pbmc4k_meta,
    method = "median"
  )
  expect_equal(nrow(pbmc4k_avg2), 2663)
})

test_that("average_clusters works when one cluster contains only 1 cell", {
  pbmc4k_meta2 <- pbmc4k_meta
  pbmc4k_meta2[1, "cluster"] <- 15
  pbmc4k_avg2 <- average_clusters(
    pbmc4k_matrix,
    pbmc4k_meta2
  )
  expect_equal(ncol(pbmc4k_avg2), 10 + 1)
})

test_that("average_clusters works when low cell number clusters should be removed", {
  pbmc4k_meta2 <- pbmc4k_meta
  pbmc4k_meta2[1, "cluster"] <- 15
  pbmc4k_avg2 <- average_clusters(pbmc4k_matrix,
    pbmc4k_meta2,
    low_threshold = 2
  )
  expect_equal(ncol(pbmc4k_avg2), 10)
})

test_that("average_clusters works when cluster info contains NA", {
  pbmc4k_meta2 <- pbmc4k_meta
  pbmc4k_meta2[1, "cluster"] <- NA
  pbmc4k_avg2 <- average_clusters(pbmc4k_matrix,
    pbmc4k_meta2,
    low_threshold = 2
  )
  expect_equal(ncol(pbmc4k_avg2), 10)
})

test_that("average_clusters_filter works on strings", {
  avg1 <- average_clusters_filter(pbmc4k_matrix, pbmc4k_meta,
    filter_on = "cluster",
    filter_method = "==",
    filter_value = "1"
  )
  remove_background(pbmc4k_matrix, avg1, 1)
  expect_equal(class(avg1), "numeric")
})

test_that("cor_to_call threshold works as intended", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )
  call1 <- cor_to_call(res,
    metadata = pbmc4k_meta,
    col = "cluster",
    collapse_to_cluster = FALSE,
    threshold = 0.5
  )

  expect_true("r<0.5, unassigned" %in% call1$type)
})

test_that("cor_to_call threshold works as intended, on per cell and collapsing", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    per_cell = T
  )
  call1 <- cor_to_call(res,
    metadata = pbmc4k_meta,
    col = "rn",
    collapse_to_cluster = "cluster",
    threshold = 0.1
  )

  expect_true(!any(is.na(call1)))
})

test_that("assign_ident works with equal length vectors and just 1 ident", {
  m1 <- assign_ident(pbmc4k_meta,
    ident_col = "classified",
    clusters = c("1", "2"),
    idents = c("whatever1", "whatever2")
  )
  m2 <- assign_ident(pbmc4k_meta,
    ident_col = "classified",
    clusters = c("1", "2"),
    idents = "whatever1"
  )
  expect_true(nrow(m1) == nrow(m2))
})

test_that("cor_to_call_topn works as intended", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )
  call1 <- cor_to_call_topn(res,
    metadata = pbmc4k_meta,
    col = "cluster",
    collapse_to_cluster = FALSE,
    threshold = 0.5
  )

  expect_true(nrow(call1) == 2 * nrow(res))
})

test_that("cor_to_call_topn works as intended on collapse to cluster option", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    per_cell = T
  )
  call1 <- cor_to_call_topn(res,
    metadata = pbmc4k_meta,
    col = "rn",
    collapse_to_cluster = "cluster",
    threshold = 0
  )

  expect_true(nrow(call1) == 2 * nrow(res))
})

test_that("gene_pct and gene_pct_markerm work as intended", {
  res <- gene_pct(
    pbmc4k_matrix,
    cbmc_m$B,
    pbmc4k_meta$cluster
  )

  res2 <- gene_pct_markerm(pbmc4k_matrix,
    cbmc_m,
    pbmc4k_meta,
    cluster_col = "cluster"
  )
  expect_error(res2 <- gene_pct_markerm(pbmc4k_matrix,
                                        cbmc_m,
                                        matrix(5,5),
                                        cluster_col = "cluster"
  ))
  expect_true(nrow(res2) == 10)
})

test_that("gene_pct_markerm norm options work", {
  res <- gene_pct_markerm(pbmc4k_matrix,
                           cbmc_m,
                           pbmc4k_meta,
                           cluster_col = "cluster",
                           norm = NULL
  )
  res2 <- gene_pct_markerm(pbmc4k_matrix,
                           cbmc_m,
                           pbmc4k_meta,
                           cluster_col = "cluster",
                           norm = "divide"
  )
  res3 <- gene_pct_markerm(pbmc4k_matrix,
                          cbmc_m,
                          pbmc4k_meta,
                          cluster_col = "cluster",
                          norm = 0.3
  )

  expect_true(nrow(res2) == 10)
})

test_that("clustify_nudge works with options and seruat2", {
  res <- clustify_nudge(
    input = s_small,
    ref_mat = cbmc_ref,
    marker = cbmc_m,
    cluster_col = "res.1",
    threshold = 0.8,
    seurat_out = F
  )
  expect_true(nrow(res) == 4)
})

test_that("clustify_nudge works with seurat_out option", {
  res <- clustify_nudge(
    input = s_small,
    ref_mat = cbmc_ref,
    marker = cbmc_m,
    cluster_col = "res.1",
    threshold = 0.8,
    seurat_out = F,
    marker_inmatrix = F
  )
  expect_true(nrow(res) == 4)
})

test_that("clustify_nudge works with list of markers", {
  res <- clustify_nudge(
    input = pbmc4k_matrix,
    ref_mat = average_clusters(pbmc4k_matrix, pbmc4k_meta),
    metadata = pbmc4k_meta,
    marker = pbmc4k_markers,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    threshold = 0.8,
    call = F,
    marker_inmatrix = F
  )
  expect_true(nrow(res) == 10)
})

test_that("overcluster_test works with ngenes option", {
  g <- overcluster_test(pbmc4k_matrix,
    pbmc4k_meta,
    cbmc_ref,
    cluster_col = "cluster"
  )
  g2 <- overcluster_test(pbmc4k_matrix,
                         pbmc4k_meta,
                         cbmc_ref,
                         cluster_col = "cluster",
                         ngenes = 100
  )
  expect_true(ggplot2::is.ggplot(g))
})

test_that("overcluster_test works with defined other cluster column", {
  g <- overcluster_test(pbmc4k_matrix,
    pbmc4k_meta,
    cbmc_ref,
    cluster_col = "cluster",
    newclustering = "classified",
    do.label = F
  )
  expect_true(ggplot2::is.ggplot(g))
})

test_that("ref_feature_select chooses the correct number of features", {
  pbmc4k_avg <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  res <- ref_feature_select(pbmc4k_avg[1:100, ], 5)
  expect_true(length(res) == 5)
})

test_that("feature_select_PCA will log transform", {
  res <- feature_select_PCA(pbmc_bulk_matrix, log_scale = F)
  res2 <- feature_select_PCA(pbmc_bulk_matrix, log_scale = T)
  expect_true(length(res) > 0)
})

test_that("feature_select_PCA can handle precalculated PCA", {
  pcs <- prcomp(t(as.matrix(pbmc_bulk_matrix)))$rotation
  res <- feature_select_PCA(pbmc_bulk_matrix, log_scale = T)
  res2 <- feature_select_PCA(pcs = pcs, log_scale = T)
  expect_true(all.equal(rownames(res), rownames(res2)))
})

test_that("downsample_matrix sets seed correctly", {
  mat1 <- downsample_matrix(pbmc4k_matrix,
                            cluster_info = pbmc4k_meta$cluster,
                            n = 0.5,
                            keep_cluster_proportions = T,
                            set_seed = 41)
  mat2 <- downsample_matrix(pbmc4k_matrix,
                            cluster_info = pbmc4k_meta$cluster,
                            n = 0.5,
                            keep_cluster_proportions = T,
                            set_seed = 41)
  expect_true(all.equal(colnames(mat1), colnames(mat2)))
})

test_that("downsample_matrix can select same number of cells per cluster", {
  mat1 <- downsample_matrix(pbmc4k_matrix,
                            cluster_info = pbmc4k_meta$cluster,
                            n = 30,
                            keep_cluster_proportions = T,
                            set_seed = 41)
  mat2 <- downsample_matrix(pbmc4k_matrix,
                            cluster_info = pbmc4k_meta$cluster,
                            n = 30,
                            keep_cluster_proportions = F,
                            set_seed = 41)

  expect_true(all.equal(ncol(mat1), 30*length(unique(pbmc4k_meta$cluster))))
})

test_that("percent_clusters works with defaults", {
  res <- percent_clusters(pbmc4k_matrix,
                          pbmc4k_meta)
  expect_equal(nrow(res), nrow(pbmc4k_matrix))
})

test_that("get_best_str finds correct values", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    per_cell = F
  )
  a <- get_best_str("0", get_best_match_matrix(res), res)
  a2 <- get_best_str("0", get_best_match_matrix(res), res, carry_cor = F)

  expect_equal(stringr::str_sub(a, 1, 3), stringr::str_sub(a2, 1, 3))
})

test_that("use_seurat_comp gets correct averages", {
  avg <- use_seurat_comp(s_small,
                         cluster_col = "res.1",
                         var.genes_only = T)
  avg2 <- use_seurat_comp(s_small,
                         cluster_col = "res.1",
                         var.genes_only = "PCA")
  expect_true(ncol(avg) == 4)
})

test_that("use_seurat_comp gets other assay slots", {
  avg <- use_seurat_comp(s_small,
                         cluster_col = "res.1",
                         assay_name = "ADT",
                         var.genes_only = T)
  avg2 <- use_seurat_comp(s_small,
                         cluster_col = "res.1",
                         assay_name = c("ADT","ADT2"),
                         var.genes_only = T)
  expect_true(nrow(avg2) - nrow(avg) == 2)
})

test_that("use_seurat_comp gets correct averages with seurat3 object", {
  avg <- use_seurat_comp(s_small3,
                         cluster_col = "RNA_snn_res.1",
                         assay_name = c("ADT","ADT2"),
                         var.genes_only = T)
  avg <- use_seurat3_comp(s_small3,
                         cluster_col = "RNA_snn_res.1",
                         assay_name = c("ADT"),
                         var.genes_only = T)
  avg2 <- use_seurat_comp(s_small3,
                          cluster_col = "RNA_snn_res.1",
                          assay_name = c("ADT","ADT2"),
                          var.genes_only = "PCA")
  expect_true(nrow(avg2) - nrow(avg) == 2)
})

test_that("clustify_intra works on test data", {
  pbmc4k_meta2 <- pbmc4k_meta
  pbmc4k_meta2$sample <- c(rep("A",150), rep("B",150))
  pbmc4k_meta2$cluster <- c(pbmc4k_meta2$classified[1:150], pbmc4k_meta2$cluster[151:300])
  res <- clustify_intra(pbmc4k_matrix,
                        pbmc4k_meta2,
                        query_genes = pbmc4k_vargenes,
                        cluster_col = "cluster",
                        sample_col = "sample",
                        sample_id = "A")
  expect_true(ncol(res) == length(unique(pbmc4k_meta2$classified[1:150])))
})
