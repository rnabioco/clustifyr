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

test_that("calculate_pathway_gsea gives appropriate output", {
  gl <- list(
    "nmd" = c("UPF1", "SMG1", "UPF2"),
    "amd" = c("ZFP36", "ZFP36L1", "ZFP36L2")
  )
  pbmc4k_avg <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  res <- calculate_pathway_gsea(pbmc4k_avg, gl, scale = T)

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
