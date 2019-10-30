context("compare_list")

test_that("warning if matrix is not binarized", {
  pbmc_mm <- matrixize_markers(pbmc_markers)
  pbmc_avg <- average_clusters(
    pbmc_matrix_small,
    pbmc_meta,
    cluster_col = "classified"
  )
  pbmc_avgb <- binarize_expr(pbmc_avg)
  gene_list_methods <- c("hyper")
  expect_warning(results <- lapply(
    gene_list_methods,
    function(x) {
      compare_lists(pbmc_avg,
        pbmc_mm,
        metric = x
      )
    }
  ))
})


test_that("run all gene list functions", {
  pbmc_mm <- matrixize_markers(pbmc_markers)
  pbmc_avg <- average_clusters(
    pbmc_matrix_small,
    pbmc_meta,
    cluster_col = "classified"
  )
  pbmc_avgb <- binarize_expr(pbmc_avg)
  gene_list_methods <- c("spearman", "hyper", "jaccard", "gsea")
  results <- lapply(
    gene_list_methods,
    function(x) {
      compare_lists(pbmc_avgb,
        pbmc_mm,
        metric = x
      )
    }
  )

  expect_equal(4, length(results))
})

test_that("gene list function options", {
  pbmc_mm <- matrixize_markers(pbmc_markers)
  pbmc_avg <- average_clusters(
    pbmc_matrix_small,
    pbmc_meta,
    cluster_col = "classified"
  )
  pbmc_avgb <- binarize_expr(pbmc_avg)
  expect_error(suppressWarnings(res <- compare_lists(pbmc_avgb,
    pbmc_mm,
    metric = "hyper",
    output_high = FALSE,
    n = 5
  )))
})

test_that("run all gene list functions in clustify_lists", {
  gene_list_methods <- c("spearman", "hyper", "jaccard", "gsea")
  results <- lapply(
    gene_list_methods,
    function(x) {
      clustify_lists(pbmc_matrix_small,
        per_cell = FALSE,
        metadata = pbmc_meta,
        cluster_col = "classified",
        marker = pbmc_markers,
        marker_inmatrix = FALSE,
        metric = x
      )
    }
  )

  expect_equal(4, length(results))
})

test_that("gsea outputs in cor matrix format", {
  res <- clustify_lists(pbmc_matrix_small,
    per_cell = FALSE,
    metadata = pbmc_meta,
    cluster_col = "classified",
    marker = pbmc_markers,
    marker_inmatrix = FALSE,
    metric = "gsea"
  )
  res2 <- cor_to_call(res)

  expect_equal(9, nrow(res2))
})

test_that("seurat object clustify_lists-ing", {
  res <- clustify_lists(s_small,
    per_cell = FALSE,
    marker = pbmc_markers,
    marker_inmatrix = FALSE,
    metric = "jaccard",
    cluster_col = "res.1",
    seurat_out = FALSE,
    dr = "tsne"
  )
  res <- clustify_lists(s_small,
    per_cell = FALSE,
    marker = pbmc_markers,
    marker_inmatrix = FALSE,
    metric = "jaccard",
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

test_that("clustify_lists inserts seurat metadata correctly", {
  res <- clustify_lists(s_small,
    per_cell = FALSE,
    marker = pbmc_markers,
    marker_inmatrix = FALSE,
    metric = "jaccard",
    cluster_col = "res.1",
    seurat_out = TRUE,
    dr = "tsne"
  )
  res2 <- clustify_lists(s_small,
    per_cell = TRUE,
    marker = pbmc_markers,
    marker_inmatrix = FALSE,
    metric = "jaccard",
    cluster_col = "res.1",
    seurat_out = TRUE,
    dr = "tsne"
  )
  expect_true(class(res) %in% c("matrix", "seurat"))
})

test_that("seurat3 object clustify_lists-ing", {
  res <- clustify_lists(s_small3,
    per_cell = FALSE,
    marker = pbmc_markers,
    marker_inmatrix = FALSE,
    metric = "jaccard",
    cluster_col = "RNA_snn_res.1",
    seurat_out = TRUE,
    dr = "tsne"
  )
  res <- clustify_lists(s_small3,
    per_cell = FALSE,
    marker = pbmc_markers,
    marker_inmatrix = FALSE,
    metric = "jaccard",
    cluster_col = "RNA_snn_res.1",
    seurat_out = FALSE,
    dr = "tsne"
  )
  g <- plot_best_call(res,
    seurat_meta(s_small3,
      dr = "tsne"
    ),
    cluster_col = "RNA_snn_res.1",
    plot_r = TRUE,
    x = "tSNE_1",
    y = "tSNE_2"
  )
  expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("clustify_lists inserts seurat3 metadata correctly", {
  res <- clustify_lists(s_small3,
    per_cell = FALSE,
    marker = pbmc_markers,
    marker_inmatrix = FALSE,
    metric = "jaccard",
    cluster_col = "RNA_snn_res.1",
    seurat_out = TRUE,
    dr = "tsne"
  )
  res2 <- clustify_lists(s_small3,
    per_cell = TRUE,
    marker = pbmc_markers,
    marker_inmatrix = FALSE,
    metric = "jaccard",
    cluster_col = "RNA_snn_res.1",
    seurat_out = TRUE,
    dr = "tsne"
  )
  expect_true(class(res) %in% c("matrix", "Seurat"))
})

test_that("run all gene list functions and then use consensus_call", {
  pbmc_mm <- matrixize_markers(pbmc_markers)
  pbmc_avg <- average_clusters(
    pbmc_matrix_small,
    pbmc_meta,
    cluster_col = "classified"
  )
  pbmc_avgb <- binarize_expr(pbmc_avg)
  gene_list_methods <- c("spearman", "hyper", "jaccard", "gsea")
  results <- lapply(
    gene_list_methods,
    function(x) {
      compare_lists(pbmc_avgb,
        pbmc_mm,
        metric = x
      )
    }
  )
  call_list <- lapply(
    results,
    cor_to_call_rank
  )
  calls <- call_consensus(call_list)
  expect_equal(4, length(results))
})

test_that("run all gene list functions in clustify_lists", {
  res <- clustify_lists(
    pbmc_matrix_small,
    cbmc_m,
    metadata = pbmc_meta,
    cluster_col = "classified",
    metric = "consensus"
  )
  expect_equal(9, nrow(res))
})

test_that("run all gene list functions in clustify_lists and seurat object", {
  res <- clustify_lists(
    s_small3,
    marker = cbmc_m,
    dr = "tsne",
    cluster_col = "RNA_snn_res.1",
    metric = "consensus",
    seurat_out = T
  )
  expect_true(is.data.frame(res) | "Seurat" %in% class(res))
})

test_that("run all gene list functions in clustify_lists and seurat object", {
  res <- clustify_lists(
    s_small,
    marker = cbmc_m,
    dr = "tsne",
    cluster_col = "res.1",
    metric = "consensus",
    seurat_out = T
  )
  expect_true(is.data.frame(res) | "seurat" %in% class(res))
})
