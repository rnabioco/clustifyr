context("clustify")

test_that("output is correctly formatted", {
  res <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    verbose = TRUE
  )
  n_clusters <- length(unique(pbmc_meta$classified))
  n_ref_samples <- ncol(pbmc_bulk_matrix)

  expect_equal(ncol(res), n_ref_samples)
  expect_equal(n_clusters, nrow(res))
})

test_that("clustify takes accidental dataframe as well", {
  res <- clustify(
    input = as.data.frame(as.matrix(pbmc_matrix_small)),
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    verbose = TRUE
  )
  n_clusters <- length(unique(pbmc_meta$classified))
  n_ref_samples <- ncol(pbmc_bulk_matrix)

  expect_equal(ncol(res), n_ref_samples)
  expect_equal(n_clusters, nrow(res))
})

test_that("clustify takes factor for metadata", {
  res <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta$classified,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    verbose = TRUE
  )
  n_clusters <- length(unique(pbmc_meta$classified))
  n_ref_samples <- ncol(pbmc_bulk_matrix)
  
  expect_equal(ncol(res), n_ref_samples)
  expect_equal(n_clusters, nrow(res))
})

test_that("run all correlation functions", {
  results <- lapply(
    clustifyr_methods,
    function(x) {
      clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = pbmc_bulk_matrix,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        compute_method = x
      )
    }
  )
  nrows <- lapply(results, nrow) %>% unlist()
  ncols <- lapply(results, ncol) %>% unlist()

  expect_equal(1, length(unique(nrows)))
  expect_equal(1, length(unique(ncols)))
})

test_that("test bad inputs", {
  expect_error(clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    compute_method = "foo"
  ))
})

test_that("test per cell", {
  res <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cell_col = "rn",
    per_cell = TRUE
  )

  expect_equal(nrow(res), ncol(pbmc_matrix_small))
})

test_that("test permutation", {
  res1 <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified"
  )

  res_full <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    n_perm = 2, return_full = TRUE
  )

  expect_equal(res1, res_full$score)
  expect_equal(length(res_full), 2)
  expect_true(all(res_full$p_val >= 0 | res_full$p_val <= 0))
})

test_that("seurat object clustifying", {
  res <- clustify(s_small,
    pbmc_bulk_matrix,
    cluster_col = "res.1",
    dr = "tsne"
  )
  res <- clustify(s_small,
    pbmc_bulk_matrix,
    cluster_col = "res.1",
    seurat_out = FALSE,
    per_cell = TRUE,
    dr = "tsne"
  )
  res <- clustify(s_small,
    pbmc_bulk_matrix,
    cluster_col = "res.1",
    seurat_out = FALSE,
    dr = "tsne"
  )
  g <- plot_best_call(res,
    seurat_meta(s_small,
      dr = "tsne"
    ),
    cluster_col = "res.1",
    plot_r = TRUE,
    x = "tSNE_1",
    y = "tSNE_2"
  )
  expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("clustify reinserts seurat metadata correctly", {
  res <- clustify(s_small,
    pbmc_bulk_matrix,
    cluster_col = "res.1",
    seurat_out = TRUE,
    per_cell = TRUE,
    dr = "tsne"
  )
  res2 <- clustify(s_small,
    pbmc_bulk_matrix,
    cluster_col = "res.1",
    seurat_out = TRUE,
    dr = "tsne"
  )
  expect_true(class(res) %in% c("matrix", "seurat"))
})

test_that("seurat3 object clustifying", {
  res <- clustify(s_small3,
    pbmc_bulk_matrix,
    cluster_col = "RNA_snn_res.1",
    dr = "tsne"
  )
  res <- clustify(s_small3,
    pbmc_bulk_matrix,
    cluster_col = "RNA_snn_res.1",
    seurat_out = FALSE,
    per_cell = TRUE,
    dr = "tsne"
  )
  res <- clustify(s_small3,
    pbmc_bulk_matrix,
    cluster_col = "RNA_snn_res.1",
    seurat_out = FALSE,
    dr = "tsne"
  )
  g <- plot_best_call(res,
    seurat_meta(s_small3,
      dr = "tsne"
    ),
    cluster_col = "RNA_snn_res.1",
    plot_r = TRUE
  )
  expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("clustify reinserts seurat3 metadata correctly", {
  res <- clustify(s_small3,
    pbmc_bulk_matrix,
    cluster_col = "RNA_snn_res.1",
    seurat_out = TRUE,
    per_cell = TRUE,
    dr = "tsne"
  )
  res2 <- clustify(s_small3,
    pbmc_bulk_matrix,
    cluster_col = "RNA_snn_res.1",
    seurat_out = TRUE,
    dr = "tsne"
  )
  expect_true(class(res) %in% c("matrix", "Seurat"))
})

test_that("get_similarity handles NA entries", {
  pbmc_meta2 <- pbmc_meta
  pbmc_meta2[1, "classified"] <- NA
  res <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta2,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified"
  )
  n_clusters <- length(unique(pbmc_meta$classified))
  n_ref_samples <- ncol(pbmc_bulk_matrix)

  expect_equal(ncol(res), n_ref_samples)
  expect_equal(n_clusters + 1, nrow(res))
})

test_that("get_similarity can exclude 0s as missing data", {
  res <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    per_cell = TRUE,
    rm0 = TRUE
  )

  expect_equal(ncol(pbmc_matrix_small), nrow(res))
})

test_that("permute_similarity runs per cell", {
  res <- permute_similarity(
    pbmc_matrix_small[c("PPBP", "LYZ", "S100A9"), c(7, 11)],
    cbmc_ref[c("PPBP", "LYZ", "S100A9"), 1:3],
    colnames(pbmc_matrix_small[c("PPBP", "LYZ", "S100A9"), c(7, 11)]),
    n_perm = 2,
    per_cell = TRUE,
    compute_method = "spearman"
  )
  expect_equal(length(res), 2)
})

test_that("error for unsupported method", {
  expect_error(res <- permute_similarity(
    pbmc_matrix_small[c("RBM28", "CCDC136", "TNPO3"), c(7, 11)],
    cbmc_ref[c("RBM28", "CCDC136", "TNPO3"), 1:3],
    pbmc_meta$rn[c(7, 11)],
    n_perm = 2,
    per_cell = TRUE,
    compute_method = "a"
  ))
})

test_that("cor throws readable error when mat has 0 rows", {
  expect_error(res <- clustify(
    input = pbmc_matrix_small[0,],
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    verbose = TRUE
  ))
})

test_that("cor throws readable error when mat has wrong number of cols", {
  expect_error(res <- clustify(
    input = pbmc_matrix_small[,1:2630],
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    verbose = TRUE
  ))
})

test_that("cor throws readable error when ref_mat has 0 rows", {
  expect_error(res <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix[0,],
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    verbose = TRUE
  ))
})

test_that("cor throws readable error when ref_mat has 0 cols", {
  expect_error(res <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix[,0],
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    verbose = TRUE
  ))
})