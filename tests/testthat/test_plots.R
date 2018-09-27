context("plotting")

res <- clustify(
  expr_mat = pbmc4k_matrix,
  metadata = pbmc4k_meta,
  bulk_mat = pbmc_bulk_matrix,
  query_genes = pbmc4k_vargenes,
  cluster_col = "cluster"
)

test_that("plots can be generated", {
  plts <- plot_cor(res, pbmc4k_meta,
           bulk_data_to_plot = colnames(res)[1:2],
           cluster_col = "cluster")

  expect_true(is.list(plts))
  expect_true(ggplot2::is.ggplot(plts[[1]]))
})

test_that("plot_gene can handle strange and normal genenames", {
  genes <- c("RP11-196G18.24",
             "RP11-442N24__B.1",
             "ACTB")
  plts <- plot_gene(pbmc4k_matrix, pbmc4k_meta,
            genes = genes, cell_col = "rn")

  expect_true(is.list(plts))
  expect_true(all(sapply(plts, ggplot2::is.ggplot)))
})
