context("run_gsea")

test_that("output is correctly formatted", {
  data("pbmc4k_vargenes")
  res <- run_gsea(pbmc4k_matrix,
    query_genes = pbmc4k_vargenes[1:100],
    n_perm = 10,
    cluster_ids = pbmc4k_meta$cluster
  )

  expect_equal(nrow(res), length(unique(pbmc4k_meta$cluster)))
  expect_true(all(res$pval >= 0 & res$pval <= 1))
})

test_that("run_gsea checks for matching number of clusters", {
  data("pbmc4k_vargenes")
  expect_error(res <- run_gsea(pbmc4k_matrix,
    query_genes = pbmc4k_vargenes[1:100],
    n_perm = 10,
    cluster_ids = pbmc4k_meta$cluster[1:3]
  ))
})

test_that("run_gsea warns slow runs", {
  data("pbmc4k_vargenes")
  res <- run_gsea(pbmc4k_matrix[, 1:3],
    query_genes = pbmc4k_vargenes[1:2],
    n_perm = 10001,
    per_cell = TRUE,
    cluster_ids = pbmc4k_meta$cluster
  )
  expect_equal(length(res), 3)
})

test_that("run_gsea warning suppression", {
  data("pbmc4k_vargenes")
  res <- run_gsea(pbmc4k_matrix[, 1:3],
    query_genes = pbmc4k_vargenes[1:2],
    n_perm = 1,
    per_cell = TRUE,
    cluster_ids = pbmc4k_meta$cluster,
    no_warnings = FALSE
  )
  expect_equal(length(res), 3)
})


test_that("calculate_pathway_gsea gives appropriate output", {
  gl <- list(
    "nmd" = c("UPF1", "SMG1", "UPF2"),
    "amd" = c("ZFP36", "ZFP36L1", "ZFP36L2")
  )
  pbmc4k_avg <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  res <- calculate_pathway_gsea(pbmc4k_avg, gl, scale = TRUE)

  expect_equal(nrow(res), length(unique(pbmc4k_meta$cluster)))
})

test_that("plot_pathway_gsea gives appropriate output", {
  gl <- list(
    "nmd" = c("UPF1", "SMG1", "UPF2"),
    "amd" = c("ZFP36", "ZFP36L1", "ZFP36L2")
  )
  pbmc4k_avg <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  g <- plot_pathway_gsea(pbmc4k_avg, gl, 5)
  expect_equal(length(g), 2)
})

test_that("plot_pathway_gsea gives output depending on returning option", {
  gl <- list(
    "nmd" = c("UPF1", "SMG1", "UPF2"),
    "amd" = c("ZFP36", "ZFP36L1", "ZFP36L2")
  )
  pbmc4k_avg <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  g <- plot_pathway_gsea(pbmc4k_avg, gl, 5, returning = "plot")
  g2 <- plot_pathway_gsea(pbmc4k_avg, gl, 5, returning = "res")
  expect_true(class(g) == "Heatmap" & class(g2) == "data.frame")
})
