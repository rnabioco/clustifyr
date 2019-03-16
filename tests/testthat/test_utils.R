context("utils")

test_that("get_vargenes works for both matrix and dataframe form", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  var1 <- get_vargenes(pbmc4k_mm)
  var2 <- get_vargenes(pbmc4k_markers)

  expect_equal(var1[1], var2[1])
})

test_that("matrixize_markers with remove_rp option", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  pbmc4k_mm2 <- matrixize_markers(pbmc4k_markers, remove_rp = T)

  expect_true(nrow(pbmc4k_mm) != nrow(pbmc4k_mm2))
})

test_that("average_clusters works as intended", {
  pbmc4k_avg2 <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  expect_equal(nrow(pbmc4k_avg2), nrow(pbmc4k_avg))
})

test_that("average_clusters works when one cluster contains only 1 cell", {
  pbmc4k_meta2 <- pbmc4k_meta
  pbmc4k_meta2[1,"cluster"] <- 15
  pbmc4k_avg2 <- average_clusters(pbmc4k_matrix, pbmc4k_meta2)
  expect_equal(ncol(pbmc4k_avg2), ncol(pbmc4k_avg) + 1)
})

test_that("average_clusters_filter works on strings", {
  avg1 <- average_clusters_filter(pbmc4k_matrix, pbmc4k_meta,
                                  filter_on = "cluster",
                                  filter_method = "==",
                                  filter_value = "1")
  remove_background(pbmc4k_matrix, avg1, 1)
  expect_equal(class(avg1), "numeric")
})

test_that("cor_to_call threshold works as intended", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )
  call1 <- cor_to_call(res,
              metadata = pbmc4k_meta,
              col = "cluster",
              collapse_to_cluster = FALSE,
              threshold = 0.5)

  expect_true("r<0.5, unassigned" %in% call1$type)
})

test_that("cor_to_call threshold works as intended, on per cell and collapsing", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    per_cell = T
  )
  call1 <- cor_to_call(res,
                       metadata = pbmc4k_meta,
                       col = "rn",
                       collapse_to_cluster = "cluster",
                       threshold = 0.1)

  expect_true(!any(is.na(call1)))
})

test_that("assign_ident works with equal length vectors and just 1 ident", {
  m1 <- assign_ident(pbmc4k_meta,
               ident_col = "classified",
               clusters = c("1","2"),
               idents = c("whatever1", "whatever2"))
  m2 <- assign_ident(pbmc4k_meta,
                     ident_col = "classified",
                     clusters = c("1","2"),
                     idents = "whatever1")
  expect_true(nrow(m1) == nrow(m2))
})

test_that("cor_to_call_topn works as intended", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )
  call1 <- cor_to_call_topn(res,
                       metadata = pbmc4k_meta,
                       col = "cluster",
                       collapse_to_cluster = FALSE,
                       threshold = 0.5)

  expect_true(nrow(call1) == 2*nrow(res))
})

test_that("cor_to_call_topn works as intended on collapse to cluster option", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster",
    per_cell = T
  )
  call1 <- cor_to_call_topn(res,
                            metadata = pbmc4k_meta,
                            col = "rn",
                            collapse_to_cluster = "cluster",
                            threshold = 0)

  expect_true(nrow(call1) == 2*nrow(res))
})

test_that("gene_pct and gene_pct_markerm work as intended", {
  res <- gene_pct(pbmc4k_matrix,
                  cbmc_m$B,
                  pbmc4k_meta$cluster)

  res2 <- gene_pct_markerm(pbmc4k_matrix,
                           cbmc_m,
                           pbmc4k_meta,
                           cluster_col = "cluster")

  expect_true(nrow(res2) == 10)
})

test_that("clustify_nudge works with options and seruat2", {
  res <- clustify_nudge(input = s_small,
                        bulk_mat = cbmc_ref,
                        marker = cbmc_m,
                        cluster_col = "res.1",
                        threshold = 0.8,
                        seurat_out = F)
  expect_true(nrow(res) == 4)
})

test_that("clustify_nudge works with options", {
  res <- clustify_nudge(input = pbmc4k_matrix,
                        bulk_mat = cbmc_ref,
                        metadata = pbmc4k_meta,
                        marker = cbmc_m,
                        query_genes = pbmc4k_vargenes,
                        cluster_col = "cluster",
                        threshold = 0.8)
  expect_true(nrow(res) == 10)
})
