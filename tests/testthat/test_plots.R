context("plot_cor")

test_that("plots can be generated", {
  res <- clustify(
    expr_mat = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )

  plts <- plot_cor(res, pbmc4k_meta,
           bulk_data_to_plot = colnames(res)[1:2],
           cluster_col = "cluster")

  expect_true(is.list(plts))
  expect_true(ggplot2::is.ggplot(plts[[1]]))
})
