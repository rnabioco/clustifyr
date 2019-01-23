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
           input = pbmc4k_matrix,
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
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    compute_method = "foo"
  ))
})

test_that("test per cell", {
  res <- clustify(
    input = pbmc4k_matrix,
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
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )

  res_full <- clustify(
    input = pbmc4k_matrix,
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

test_that("run all gene list functions", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  pbmc4k_avgb <- binarize_expr(pbmc4k_avg)
  gene_list_methods <- c("spearman", "hyper", "jaccard", "gsea")
  results <- lapply(gene_list_methods,
                    function(x) {
                      compare_lists(pbmc4k_avgb,
                                    pbmc4k_mm,
                                    metric = x)})

  expect_equal(4, length(results))
})

test_that("gene list function options", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  pbmc4k_avgb <- binarize_expr(pbmc4k_avg)
  expect_error(res <- compare_lists(pbmc4k_avgb,
                       pbmc4k_mm,
                       metric = "hyper",
                       output_high = F,
                       n = 5))
})
