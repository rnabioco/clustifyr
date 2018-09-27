context("run_gsea")

test_that("output is correctly formatted", {
  data("pbmc4k_vargenes")
  res <- run_gsea(pbmc4k_matrix,
                  query_genes = pbmc4k_vargenes[1:100],
                  n_perm = 10,
                  cluster_ids = pbmc4k_meta$cluster)

  expect_equal(nrow(res), length(unique(pbmc4k_meta$cluster)))
  expect_true(all(res$pval >= 0 & res$pval <= 1))
  })
