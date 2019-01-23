context("compare_list")

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
