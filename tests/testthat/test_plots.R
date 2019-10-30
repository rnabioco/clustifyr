context("plotting")

res <- clustify(
  input = pbmc_matrix_small,
  metadata = pbmc_meta,
  ref_mat = pbmc_bulk_matrix,
  query_genes = pbmc_vargenes,
  cluster_col = "classified"
)

res2 <- clustify(
  input = pbmc_matrix_small,
  metadata = pbmc_meta,
  ref_mat = pbmc_bulk_matrix,
  query_genes = pbmc_vargenes,
  cluster_col = "classified",
  per_cell = TRUE
)

test_that("plots can be generated", {
  plts <- plot_best_call(
    res,
    pbmc_meta,
    cluster_col = "classified"
  )
  expect_true(ggplot2::is.ggplot(plts))
})

test_that("plot_best_call warns about colnames", {
  pbmc_meta2 <- pbmc_meta
  pbmc_meta2$type <- 1
  plts <- plot_best_call(res, pbmc_meta2)
  expect_true(is.null(plts))
})

test_that("call plots can be generated", {
  plts <- plot_cor(res, pbmc_meta,
    data_to_plot = colnames(res)[1:2],
    cluster_col = "classified"
  )

  expect_error(plts <- plot_cor(res, pbmc_meta,
    data_to_plot = "nonsense",
    cluster_col = "classified"
  ))

  expect_true(is.list(plts))
  expect_true(ggplot2::is.ggplot(plts[[1]]))
})

test_that("plot_cor for all clusters by default", {
  plts <- plot_cor(res,
    pbmc_meta,
    cluster_col = "classified",
    x = "UMAP_1",
    y = "UMAP_2"
  )

  plts2 <- plot_cor(res2,
    pbmc_meta %>% tibble::rownames_to_column("rn"),
    cluster_col = "rn",
    x = "UMAP_1",
    y = "UMAP_2"
  )

  expect_true(length(plts) == 14)
})

test_that("plot_cor works with scale_legends option", {
  plts <- plot_cor(res,
    pbmc_meta,
    cluster_col = "classified",
    scale_legends = TRUE
  )

  plts2 <- plot_cor(res,
    pbmc_meta,
    cluster_col = "classified",
    scale_legends = c(0, 1)
  )
  expect_true(length(plts) == 14)
})

test_that("plot_gene can handle strange and normal genenames", {
  genes <- c(
    "RP11-314N13.3",
    "ARF4"
  )
  plts <- plot_gene(
    pbmc_matrix_small,
    pbmc_meta %>% tibble::rownames_to_column("rn"),
    genes = genes,
    cell_col = "rn"
  )

  expect_true(is.list(plts))
  expect_true(all(sapply(plts, ggplot2::is.ggplot)))
})

test_that("plot_gene automatically plots all cells", {
  genes <- c(
    "ZYX"
  )
  expect_error(plts <- plot_gene(pbmc_matrix_small,
    tibble::column_to_rownames(pbmc_meta, "rn"),
    genes = genes,
    cell_col = "nonsense"
  ))

  plts <- plot_gene(
    pbmc_matrix_small,
    pbmc_meta,
    genes = genes
  )

  expect_true(all(sapply(plts, ggplot2::is.ggplot)))
})

test_that("plot_best_call threshold works as intended, on per cell and collapsing", {
  res <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    per_cell = TRUE
  )
  call1 <- plot_best_call(res,
    metadata = pbmc_meta,
    per_cell = TRUE,
    collapse_to_cluster = "classified",
    threshold = 0.3
  )

  expect_true(ggplot2::is.ggplot(call1))
})

test_that("plot_gene checks for presence of gene name", {
  plot_gene(pbmc_matrix_small,
    pbmc_meta %>% tibble::rownames_to_column("rn"),
    c("INIP", "ZFP36L3"),
    cell_col = "rn",
    do_label = TRUE,
    do_legend = FALSE,
    x = "UMAP_1",
    y = "UMAP_2"
  )
  expect_error(plot_gene(pbmc_matrix_small,
    pbmc_meta %>% tibble::rownames_to_column("rn"),
    c("ZFP36L3"),
    cell_col = "rn",
    x = "UMAP_1",
    y = "UMAP_2"
  ))
})

test_that("plot_cols returns a ggplot object", {
  g <- plot_cols(
    pbmc_meta,
    "seurat_clusters",
    "classified",
    "UMAP_1",
    pbmc_meta,
    "classified",
    "UMAP_1"
  )
  expect_true(ggplot2::is.ggplot(g))
})

test_that("plot_cor_heatmap returns a ggplot object", {
  res <- clustify(
    input = pbmc_matrix_small,
    metadata = pbmc_meta,
    ref_mat = pbmc_bulk_matrix,
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    per_cell = FALSE
  )
  g <- plot_cor_heatmap(res)
  expect_true(is(g, "Heatmap"))
})

test_that("plot_call works on defaults", {
  g <- plot_call(res,
    pbmc_meta,
    cluster_col = "classified"
  )

  expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("plot_tsne works with alpha_col", {
  pbmc_meta2 <- pbmc_meta
  pbmc_meta2$al <- 0
  pbmc_meta2$al[1] <- 1 # 1:nrow(pbmc_meta)/nrow(pbmc_meta)
  g <- plot_tsne(pbmc_meta2,
    feature = "classified",
    alpha_col = "al", do_legend = FALSE
  )
  expect_true(ggplot2::is.ggplot(g))
})
