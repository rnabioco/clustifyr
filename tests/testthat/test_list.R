context("compare_list")

test_that("warning if matrix is not binarized", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  pbmc4k_avg <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  pbmc4k_avgb <- binarize_expr(pbmc4k_avg)
  gene_list_methods <- c("hyper")
  results <- lapply(
    gene_list_methods,
    function(x) {
      compare_lists(pbmc4k_avg,
                    pbmc4k_mm,
                    metric = x
      )
    }
  )

  expect_true(results[[1]][1,1] > 0)
})


test_that("run all gene list functions", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  pbmc4k_avg <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  pbmc4k_avgb <- binarize_expr(pbmc4k_avg)
  gene_list_methods <- c("spearman", "hyper", "jaccard", "gsea")
  results <- lapply(
    gene_list_methods,
    function(x) {
      compare_lists(pbmc4k_avgb,
        pbmc4k_mm,
        metric = x
      )
    }
  )

  expect_equal(4, length(results))
})

test_that("gene list function options", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  pbmc4k_avg <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  pbmc4k_avgb <- binarize_expr(pbmc4k_avg)
  expect_error(suppressWarnings(res <- compare_lists(pbmc4k_avgb,
    pbmc4k_mm,
    metric = "hyper",
    output_high = F,
    n = 5
  )))
})

test_that("run all gene list functions in clustify_lists", {
  gene_list_methods <- c("spearman", "hyper", "jaccard", "gsea")
  results <- lapply(
    gene_list_methods,
    function(x) {
      clustify_lists(pbmc4k_matrix,
        per_cell = F,
        cluster_info = pbmc4k_meta,
        marker = pbmc4k_markers,
        marker_inmatrix = F,
        metric = x
      )
    }
  )

  expect_equal(4, length(results))
})

test_that("gsea outputs in cor matrix format", {
  res <- clustify_lists(pbmc4k_matrix,
                        per_cell = F,
                        cluster_info = pbmc4k_meta,
                        marker = pbmc4k_markers,
                        marker_inmatrix = F,
                        metric = "gsea")
  res2 <- cor_to_call(res)

  expect_equal(10, nrow(res2))
})

test_that("seurat object clustify_lists-ing", {
  res <- clustify_lists(s_small,
    per_cell = F,
    marker = pbmc4k_markers,
    marker_inmatrix = F,
    metric = "jaccard",
    cluster_col = "res.1",
    seurat_out = F
  )
  res <- clustify_lists(s_small,
    per_cell = F,
    marker = pbmc4k_markers,
    marker_inmatrix = F,
    metric = "jaccard",
    cluster_col = "res.1",
    seurat_out = F
  )
  g <- plot_best_call(res,
    use_seurat_meta(s_small),
    col = "res.1",
    plot_r = T
  )
  expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("seurat3 object clustify_lists-ing", {
  res <- clustify_lists(s_small3,
                        per_cell = F,
                        marker = pbmc4k_markers,
                        marker_inmatrix = F,
                        metric = "jaccard",
                        cluster_col = "RNA_snn_res.1",
                        seurat_out = T
  )
  res <- clustify_lists(s_small3,
                        per_cell = F,
                        marker = pbmc4k_markers,
                        marker_inmatrix = F,
                        metric = "jaccard",
                        cluster_col = "RNA_snn_res.1",
                        seurat_out = F
  )
  g <- plot_best_call(res,
                      use_seurat_meta(s_small3),
                      col = "RNA_snn_res.1",
                      plot_r = T
  )
  expect_true(ggplot2::is.ggplot(g[[1]]))
})
