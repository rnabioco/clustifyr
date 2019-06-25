context("plotting")

res <- clustify(
  input = pbmc4k_matrix,
  metadata = pbmc4k_meta,
  ref_mat = pbmc_bulk_matrix,
  query_genes = pbmc4k_vargenes,
  cluster_col = "cluster"
)

res2 <- clustify(
  input = pbmc4k_matrix,
  metadata = pbmc4k_meta,
  ref_mat = pbmc_bulk_matrix,
  query_genes = pbmc4k_vargenes,
  cluster_col = "cluster",
  per_cell = TRUE
)

test_that("plots can be generated", {
  plts <- plot_best_call(res, pbmc4k_meta)
  expect_true(ggplot2::is.ggplot(plts))
})

test_that("plot_best_call warns about colnames", {
  pbmc4k_meta2 <- pbmc4k_meta
  pbmc4k_meta2$type <- 1
  plts <- plot_best_call(res, pbmc4k_meta2)
  expect_true(is.null(plts))
})

test_that("call plots can be generated", {
  plts <- plot_cor(res, pbmc4k_meta,
    ref_data_to_plot = colnames(res)[1:2],
    cluster_col = "cluster"
  )

  expect_error(plts <- plot_cor(res, pbmc4k_meta,
                                ref_data_to_plot = "nonsense",
                                cluster_col = "cluster"
  ))

  expect_true(is.list(plts))
  expect_true(ggplot2::is.ggplot(plts[[1]]))
})

test_that("plot_cor for all clusters by default", {
  plts <- plot_cor(res,
                   pbmc4k_meta,
                   cluster_col = "cluster"
  )

  plts2 <- plot_cor(res2,
                   pbmc4k_meta,
                   cluster_col = "cluster"
  )

  expect_true(length(plts) == 14)
})

test_that("plot_cor works with scale_legends option", {
  plts <- plot_cor(res,
                   pbmc4k_meta,
                   cluster_col = "cluster",
                   scale_legends = TRUE
  )

  plts2 <- plot_cor(res,
                   pbmc4k_meta,
                   cluster_col = "cluster",
                   scale_legends = c(0,1)
  )
  expect_true(length(plts) == 14)
})

test_that("plot_gene can handle strange and normal genenames", {
  genes <- c(
    "RP4-673M15.1",
    "CD24"
  )
  plts <- plot_gene(pbmc4k_matrix, pbmc4k_meta,
    genes = genes, cell_col = "rn"
  )

  expect_true(is.list(plts))
  expect_true(all(sapply(plts, ggplot2::is.ggplot)))
})

test_that("plot_gene automatically plots all cells", {
  genes <- c(
    "TNPO3"
  )
  expect_error(plts <- plot_gene(pbmc4k_matrix,
                                 tibble::column_to_rownames(pbmc4k_meta, "rn"),
                                 genes = genes,
                                 cell_col = "nonsense"
  ))

  plts <- plot_gene(pbmc4k_matrix,
                    tibble::column_to_rownames(pbmc4k_meta, "rn"),
                    genes = genes
  )

  expect_true(all(sapply(plts, ggplot2::is.ggplot)))
})

test_that("plot_best_call threshold works as intended, on per cell and collapsing", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    per_cell = TRUE
  )
  call1 <- plot_best_call(res,
    metadata = pbmc4k_meta,
    col = "rn",
    collapse_to_cluster = "cluster",
    threshold = 0.3
  )

  expect_true(ggplot2::is.ggplot(call1))
})

test_that("plot_gene checks for presence of gene name", {
  plot_gene(pbmc4k_matrix,
    pbmc4k_meta,
    c("CD24", "ZFP36L3"),
    cell_col = "rn",
    do.label = TRUE,
    do.legend = FALSE
  )
  expect_error(plot_gene(pbmc4k_matrix,
    pbmc4k_meta,
    c("ZFP36L3"),
    cell_col = "rn"
  ))
})

test_that("plot_cols returns a ggplot object", {
  g <- plot_cols(
    pbmc4k_meta,
    "cluster",
    "classified",
    "tSNE_1",
    pbmc4k_meta,
    "cluster",
    "tSNE_1"
  )
  expect_true(ggplot2::is.ggplot(g))
})

test_that("plot_cor_heatmap returns a ggplot object", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    per_cell = FALSE
  )
  g <- plot_cor_heatmap(res)
  expect_true(class(g) == "Heatmap")
})

test_that("plot_call works on defaults", {
  g <- plot_call(res,
                 pbmc4k_meta,
                 cluster_col = "cluster")

  expect_true(ggplot2::is.ggplot(g[[1]]))
})
