context("run_gsea")

# use capture.output to quiet progress bar from fgsea
shush <- function(...) {
  capture.output(..., file = nullfile())
}
test_that("output is correctly formatted", {
    data("pbmc_vargenes")
   
    shush(res <- run_gsea(
        pbmc_matrix_small,
        query_genes = pbmc_vargenes[1:100],
        n_perm = 10,
        cluster_ids = pbmc_meta$classified,
        no_warnings = TRUE
    ))

    expect_equal(nrow(res), length(unique(pbmc_meta$classified)))
    expect_true(all(res$pval >= 0 & res$pval <= 1))
})

test_that("run_gsea checks for matching number of clusters", {
    data("pbmc_vargenes")
    expect_error(
        res <- run_gsea(
            pbmc_matrix_small,
            query_genes = pbmc_vargenes[1:100],
            n_perm = 10,
            cluster_ids = pbmc_meta$classified[1:3],
            no_warnings = TRUE
        )
    )
})

test_that("run_gsea warns slow runs", {
    data("pbmc_vargenes")

    expect_warning(
      shush(res <- run_gsea(pbmc_matrix_small[, 1:3],
        query_genes = pbmc_vargenes[1:2],
        n_perm = 10001,
        per_cell = TRUE,
        cluster_ids = pbmc_meta$classified,
        no_warnings = TRUE
        )
      ) 
    )
})

test_that("run_gsea warning suppression", {
    data("pbmc_vargenes")
    expect_warning(
      shush(res <- run_gsea(
            pbmc_matrix_small[, 1:3],
            query_genes = pbmc_vargenes[1:2],
            n_perm = 1,
            per_cell = TRUE,
            cluster_ids = pbmc_meta$classified,
            no_warnings = FALSE
        )
      ) 
    )
})

test_that("calculate_pathway_gsea gives appropriate output", {
    gl <- list(
        "n" = c("PPBP", "LYZ", "S100A9"),
        "a" = c("IGLL5", "GNLY", "FTL")
    )
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    shush(res <- calculate_pathway_gsea(pbmc_avg, gl, scale = TRUE))

    expect_equal(nrow(res), length(unique(pbmc_meta$classified)))
})

test_that("plot_pathway_gsea gives appropriate output", {
    gl <- list(
        "n" = c("PPBP", "LYZ", "S100A9"),
        "a" = c("IGLL5", "GNLY", "FTL")
    )
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    shush(g <- plot_pathway_gsea(pbmc_avg, gl, 5))
    expect_equal(length(g), 2)
})

test_that("plot_pathway_gsea gives output depending on returning option", {
    gl <- list(
        "n" = c("PPBP", "LYZ", "S100A9"),
        "a" = c("IGLL5", "GNLY", "FTL")
    )
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    shush(g <- plot_pathway_gsea(pbmc_avg, gl, 5, returning = "plot"))
    shush(g2 <- plot_pathway_gsea(pbmc_avg, gl, 5, returning = "res"))
    
    expect_true(is(g, "Heatmap") & is.data.frame(g2))
})
