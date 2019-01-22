context("clustify")

test_that("output is correctly formatted", {
  res <- clustify(
           input = pbmc4k_matrix,
           metadata = pbmc4k_meta,
           bulk_mat = pbmc_bulk_matrix,
           query_genes = pbmc4k_vargenes,
           cluster_col = "cluster"
  )
  n_clusters <- length(unique(pbmc4k_meta$cluster))
  n_bulk_samples <- ncol(pbmc_bulk_matrix)

  expect_equal(ncol(res), n_bulk_samples)
  expect_equal(n_clusters, nrow(res))
  })

test_that("run all correlation functions", {
  results <- lapply(clustifyr_methods,
         function(x) {
           clustify(
           expr_mat = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    compute_method = x
  )})
  nrows <- lapply(results, nrow) %>% unlist()
  ncols <- lapply(results, ncol) %>% unlist()

  expect_equal(1, length(unique(nrows)))
  expect_equal(1, length(unique(ncols)))
})

test_that("test bad inputs", {
  expect_error(clustify(
    expr_mat = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    compute_method = "foo"
  ))
})

test_that("test per cell", {
  res <- clustify(
    expr_mat = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cell_col = "rn",
    per_cell = T
  )

  expect_equal(nrow(res), ncol(pbmc4k_matrix))
})

test_that("test permutation", {

  res1 <- clustify(
    expr_mat = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )

  res_full <- clustify(
    expr_mat = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    num_perm = 2, return_full = T
  )

  expect_equal(res1, res_full$score)
  expect_equal(length(res_full), 2)
  expect_true(all(res_full$p_val >= 0 | res_full$p_val <= 0))
})
