context("clustify")

test_that("output is correctly formatted", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    verbose = TRUE
  )
  n_clusters <- length(unique(pbmc4k_meta$cluster))
  n_ref_samples <- ncol(pbmc_bulk_matrix)

  expect_equal(ncol(res), n_ref_samples)
  expect_equal(n_clusters, nrow(res))
})

test_that("clustify takes accidental dataframe as well", {
  res <- clustify(
    input = as.data.frame(as.matrix(pbmc4k_matrix)),
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    verbose = TRUE
  )
  n_clusters <- length(unique(pbmc4k_meta$cluster))
  n_ref_samples <- ncol(pbmc_bulk_matrix)

  expect_equal(ncol(res), n_ref_samples)
  expect_equal(n_clusters, nrow(res))
})

test_that("run all correlation functions", {
  results <- lapply(
    clustifyr_methods,
    function(x) {
      clustify(
        input = pbmc4k_matrix,
        metadata = pbmc4k_meta,
        ref_mat = pbmc_bulk_matrix,
        query_genes = pbmc4k_vargenes,
        cluster_col = "cluster",
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
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    compute_method = "foo"
  ))
})

test_that("test per cell", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cell_col = "rn",
    per_cell = TRUE
  )

  expect_equal(nrow(res), ncol(pbmc4k_matrix))
})

test_that("test permutation", {
  res1 <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )

  res_full <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    num_perm = 2, return_full = TRUE
  )

  expect_equal(res1, res_full$score)
  expect_equal(length(res_full), 2)
  expect_true(all(res_full$p_val >= 0 | res_full$p_val <= 0))
})

test_that("seurat object clustifying", {
  res <- clustify(s_small,
    pbmc_bulk_matrix,
    cluster_col = "res.1"
  )
  res <- clustify(s_small,
    pbmc_bulk_matrix,
    cluster_col = "res.1",
    seurat_out = FALSE,
    per_cell = TRUE
  )
  res <- clustify(s_small,
    pbmc_bulk_matrix,
    cluster_col = "res.1",
    seurat_out = FALSE
  )
  g <- plot_best_call(res,
    use_seurat_meta(s_small),
    col = "res.1",
    plot_r = TRUE
  )
  expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("clustify reinserts seurat metadata correctly", {
  res <- clustify(s_small,
                  pbmc_bulk_matrix,
                  cluster_col = "res.1",
                  seurat_out = TRUE,
                  per_cell = TRUE
  )
  res2 <- clustify(s_small,
                  pbmc_bulk_matrix,
                  cluster_col = "res.1",
                  seurat_out = TRUE
  )
  expect_true(class(res) %in% c("matrix", "seurat"))
})

test_that("seurat3 object clustifying", {
  res <- clustify(s_small3,
                  pbmc_bulk_matrix,
                  cluster_col = "RNA_snn_res.1"
  )
  res <- clustify(s_small3,
                  pbmc_bulk_matrix,
                  cluster_col = "RNA_snn_res.1",
                  seurat_out = FALSE,
                  per_cell = TRUE
  )
  res <- clustify(s_small3,
                  pbmc_bulk_matrix,
                  cluster_col = "RNA_snn_res.1",
                  seurat_out = FALSE
  )
  g <- plot_best_call(res,
                      use_seurat_meta(s_small3),
                      col = "RNA_snn_res.1",
                      plot_r = TRUE
  )
  expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("clustify reinserts seurat3 metadata correctly", {
  res <- clustify(s_small3,
                  pbmc_bulk_matrix,
                  cluster_col = "RNA_snn_res.1",
                  seurat_out = TRUE,
                  per_cell = TRUE
  )
  res2 <- clustify(s_small3,
                   pbmc_bulk_matrix,
                   cluster_col = "RNA_snn_res.1",
                   seurat_out = TRUE
  )
  expect_true(class(res) %in% c("matrix", "Seurat"))
})

test_that("get_similarity handles NA entries", {
  pbmc4k_meta2 <- pbmc4k_meta
  pbmc4k_meta2[1, "cluster"] <- NA
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta2,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )
  n_clusters <- length(unique(pbmc4k_meta$cluster))
  n_ref_samples <- ncol(pbmc_bulk_matrix)

  expect_equal(ncol(res), n_ref_samples)
  expect_equal(n_clusters + 1, nrow(res))
})

test_that("get_similarity can exclude 0s as missing data", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    per_cell = TRUE,
    rm0 = TRUE
  )

  expect_equal(ncol(pbmc4k_matrix), nrow(res))
})

test_that("permute_similarity runs per cell", {
  res <- permute_similarity(
    pbmc4k_matrix[c("RBM28","CCDC136","TNPO3"),c(7,11)],
    cbmc_ref[c("RBM28","CCDC136","TNPO3"),1:3],
    colnames(pbmc4k_matrix[c("RBM28","CCDC136","TNPO3"),c(7,11)]),
    num_perm = 2,
    per_cell = TRUE,
    compute_method = "spearman"
  )
  expect_equal(length(res), 2)
})

test_that("error for unsupported method", {
  expect_error(res <- permute_similarity(
    pbmc4k_matrix[c("RBM28","CCDC136","TNPO3"),c(7,11)],
    cbmc_ref[c("RBM28","CCDC136","TNPO3"),1:3],
    pbmc4k_meta$rn[c(7,11)],
    num_perm = 2,
    per_cell = TRUE,
    compute_method = "a"
  ))
})
