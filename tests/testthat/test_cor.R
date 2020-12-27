context("clustify")

test_that("output is correctly formatted", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        verbose = TRUE
    )
    n_clusters <- length(unique(pbmc_meta$classified))
    n_ref_samples <- ncol(cbmc_ref)

    expect_equal(ncol(res), n_ref_samples)
    expect_equal(n_clusters, nrow(res))
})

test_that("clustify takes accidental dataframe as well", {
    res <- clustify(
        input = as.data.frame(as.matrix(pbmc_matrix_small)),
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        verbose = TRUE
    )
    n_clusters <- length(unique(pbmc_meta$classified))
    n_ref_samples <- ncol(cbmc_ref)

    expect_equal(ncol(res), n_ref_samples)
    expect_equal(n_clusters, nrow(res))
})

test_that("clustify takes factor for metadata", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta$classified,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        verbose = TRUE
    )
    n_clusters <- length(unique(pbmc_meta$classified))
    n_ref_samples <- ncol(cbmc_ref)

    expect_equal(ncol(res), n_ref_samples)
    expect_equal(n_clusters, nrow(res))
})

test_that("run all correlation functions", {
    results <- lapply(
        clustifyr_methods,
        function(x) {
            clustify(
                input = pbmc_matrix_small,
                metadata = pbmc_meta,
                ref_mat = cbmc_ref,
                query_genes = pbmc_vargenes,
                cluster_col = "classified",
                compute_method = x
            )
        }
    )
    nrows <- lapply(results, nrow) %>% unlist()
    ncols <- lapply(results, ncol) %>% unlist()

    expect_equal(1, length(unique(nrows)))
    expect_equal(1, length(unique(ncols)))
})

test_that("test bad inputs", {
    expect_error(clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        compute_method = "foo"
    ))
})

test_that("test per cell", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cell_col = "rn",
        per_cell = TRUE
    )

    expect_equal(nrow(res), ncol(pbmc_matrix_small))
})

test_that("test permutation", {
    res1 <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )

    res_full <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        n_perm = 2,
        return_full = TRUE
    )

    expect_equal(res1, res_full$score)
    expect_equal(length(res_full), 2)
    expect_true(all(res_full$p_val >= 0 | res_full$p_val <= 0))
})

test_that("seurat object clustifying", {
    res <- clustify(s_small,
        cbmc_ref,
        cluster_col = "res.1",
        dr = "tsne"
    )

    res <- clustify(s_small,
        cbmc_ref,
        cluster_col = "res.1",
        seurat_out = FALSE,
        per_cell = TRUE,
        dr = "tsne"
    )

    res <- clustify(s_small,
        cbmc_ref,
        cluster_col = "res.1",
        seurat_out = FALSE,
        dr = "tsne"
    )
    g <- plot_best_call(
        res,
        seurat_meta(s_small,
            dr = "tsne"
        ),
        cluster_col = "res.1",
        plot_r = TRUE,
        x = "tSNE_1",
        y = "tSNE_2"
    )
    expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("clustify reinserts seurat metadata correctly", {
    res <- clustify(s_small,
        cbmc_ref,
        cluster_col = "res.1",
        seurat_out = TRUE,
        per_cell = TRUE,
        dr = "tsne"
    )
    res2 <- clustify(s_small,
        cbmc_ref,
        cluster_col = "res.1",
        seurat_out = TRUE,
        dr = "tsne"
    )
    if ("Seurat" %in% loadedNamespaces()) {
        expect_true(class(res) %in% c("seurat"))
    } else {
        expect_true(is.matrix(res))
    }
})

test_that("seurat3 object clustifying", {
    res <- clustify(s_small3,
        cbmc_ref,
        cluster_col = "RNA_snn_res.1",
        dr = "tsne"
    )
    res <- clustify(s_small3,
        cbmc_ref,
        cluster_col = "RNA_snn_res.1",
        seurat_out = FALSE,
        per_cell = TRUE,
        dr = "tsne"
    )
    res <- clustify(s_small3,
        cbmc_ref,
        cluster_col = "RNA_snn_res.1",
        seurat_out = FALSE,
        dr = "tsne"
    )
    g <- plot_best_call(res,
        seurat_meta(s_small3,
            dr = "tsne"
        ),
        cluster_col = "RNA_snn_res.1",
        plot_r = TRUE
    )
    expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("object with passing vector as metadata", {
    res <- clustify(s_small3,
                    cbmc_ref,
                    metadata = s_small3@meta.data$RNA_snn_res.1,
                    dr = "tsne"
    )
    res <- clustify_lists(
        s_small3,
        marker = cbmc_m,
        metadata = s_small3@meta.data$RNA_snn_res.1,
        dr = "tsne",
        metric = "posneg",
        seurat_out = FALSE
    )
    res <- clustify(s_small,
                    cbmc_ref,
                    metadata = s_small@meta.data$res.1,
                    dr = "tsne"
    )
    res <- clustify_lists(
        s_small,
        marker = cbmc_m,
        metadata = s_small@meta.data$res.1,
        dr = "tsne",
        metric = "posneg",
        seurat_out = FALSE
    )
    expect_true(is.matrix(res))
})

test_that("clustify reinserts seurat3 metadata correctly", {
    res <- clustify(s_small3,
        cbmc_ref,
        cluster_col = "RNA_snn_res.1",
        seurat_out = TRUE,
        per_cell = TRUE,
        dr = "tsne"
    )

    res2 <- clustify(s_small3,
        cbmc_ref,
        cluster_col = "RNA_snn_res.1",
        seurat_out = TRUE,
        dr = "tsne"
    )
    if ("Seurat" %in% loadedNamespaces()) {
        expect_true(class(res) %in% c("Seurat"))
    } else {
        expect_true(is.matrix(res))
    }
})

test_that("get_similarity handles NA entries", {
    pbmc_meta2 <- pbmc_meta
    pbmc_meta2[1, "classified"] <- NA
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta2,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    n_clusters <- length(unique(pbmc_meta$classified))
    n_ref_samples <- ncol(cbmc_ref)

    expect_equal(ncol(res), n_ref_samples)
    expect_equal(n_clusters + 1, nrow(res))
})

test_that("get_similarity handles NA entries", {
    pbmc_meta2 <- pbmc_meta
    pbmc_meta2[1, "classified"] <- NA
    res <- get_similarity(
        pbmc_matrix_small[intersect(rownames(pbmc_matrix_small), rownames(cbmc_ref)), ],
        ref_mat = cbmc_ref[intersect(rownames(pbmc_matrix_small), rownames(cbmc_ref)), ],
        pbmc_meta2$classified,
        compute_method = "spearman"
    )
    n_clusters <- length(unique(pbmc_meta$classified))
    n_ref_samples <- ncol(cbmc_ref)

    expect_equal(ncol(res), n_ref_samples)
    expect_equal(n_clusters + 1, nrow(res))
})

test_that("get_similarity can exclude 0s as missing data", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        per_cell = TRUE,
        rm0 = TRUE
    )

    expect_equal(ncol(pbmc_matrix_small), nrow(res))
})

test_that("permute_similarity runs per cell", {
    res <- permute_similarity(
        pbmc_matrix_small[c("PPBP", "LYZ", "S100A9"), c(7, 11)],
        cbmc_ref[c("PPBP", "LYZ", "S100A9"), 1:3],
        colnames(pbmc_matrix_small[c("PPBP", "LYZ", "S100A9"), c(7, 11)]),
        n_perm = 2,
        per_cell = TRUE,
        compute_method = "spearman"
    )
    expect_equal(length(res), 2)
})

test_that("error for unsupported method", {
    expect_error(
        res <- permute_similarity(
            pbmc_matrix_small[c("RBM28", "CCDC136", "TNPO3"), c(7, 11)],
            cbmc_ref[c("RBM28", "CCDC136", "TNPO3"), 1:3],
            pbmc_meta$rn[c(7, 11)],
            n_perm = 2,
            per_cell = TRUE,
            compute_method = "a"
        )
    )
})

test_that("cor throws readable error when mat has 0 rows", {
    expect_error(res <- clustify(
        input = pbmc_matrix_small[0, ],
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        verbose = TRUE
    ))
})

test_that("cor throws readable error when mat has wrong number of cols", {
    expect_error(res <- clustify(
        input = pbmc_matrix_small[, 1:2630],
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        verbose = TRUE
    ))
})

test_that("cor throws readable error when ref_mat has 0 rows", {
    expect_error(res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref[0, ],
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        verbose = TRUE
    ))
})

test_that("cor throws readable error when ref_mat has 0 cols", {
    expect_error(res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref[, 0],
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        verbose = TRUE
    ))
})

test_that("sparse matrix is accepted as input", {
    res <- clustify(
        input = s_small3@assays$RNA@counts,
        metadata = s_small3@meta.data,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "letter.idents",
        verbose = TRUE
    )

    expect_equal(2, nrow(res))
})

test_that("correct error message is displayed for nonexistent cluster_col", {
    expect_error(res <- clustify(
        input = s_small3@assays$RNA@counts,
        metadata = s_small3@meta.data,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "a",
        verbose = TRUE
    ))
})

test_that("input Seurat metadata columns are not changed (type, r, rn, etc). #259", {
    skip_if_not_installed("Seurat")
    tmp <- s_small3
    tmp@meta.data$type <- 0L
    tmp@meta.data$rn <- 0L
    tmp@meta.data$r <- 0L

    res <- clustify(
        input = tmp,
        ref_mat = cbmc_ref,
        cluster_col = "RNA_snn_res.1",
        dr = "tsne"
    )

    expect_true(all(c("type", "rn", "r") %in% colnames(res@meta.data)))
    expect_true(all(res@meta.data$type == 0L))
    expect_true(all(res@meta.data$rn == 0L))
    expect_true(all(res@meta.data$r == 0L))
})

test_that("clustify_lists works with pos_neg_select and Seurat3 object", {
    res <- clustify_lists(
        s_small3,
        marker = cbmc_m,
        cluster_col = "RNA_snn_res.1",
        dr = "tsne",
        metric = "posneg",
        seurat_out = FALSE
    )
    expect_true(nrow(res) == 3)
})

test_that("clustify_lists works with pos_neg_select, Seurat3 object, and lists of genes", {
    res <- clustify_lists(
        s_small3,
        marker = as.list(cbmc_m),
        marker_inmatrix = FALSE,
        cluster_col = "RNA_snn_res.1",
        dr = "tsne",
        metric = "posneg",
        seurat_out = FALSE
    )
    expect_true(nrow(res) == 3)
})

test_that("clustify_lists works with pos_neg_select, Seurat3 object, and matrix preprocessed by pos_neg_marker", {
    res <- clustify_lists(
        s_small3,
        marker = pos_neg_marker(as.list(cbmc_m)),
        marker_inmatrix = FALSE,
        cluster_col = "RNA_snn_res.1",
        dr = "tsne",
        metric = "posneg",
        seurat_out = FALSE
    )
    expect_true(nrow(res) == 3)
})

test_that("clustify_lists works with pct and Seurat3 object", {
    res <- clustify_lists(
        s_small3,
        marker = cbmc_m,
        cluster_col = "RNA_snn_res.1",
        dr = "tsne",
        metric = "pct",
        seurat_out = FALSE
    )
    expect_true(nrow(res) == 3)
})

test_that("clustify_lists gives correct error message upon unrecognized method", {
    expect_error(
        res <- clustify_lists(
            s_small3,
            marker = cbmc_m,
            cluster_col = "RNA_snn_res.1",
            dr = "tsne",
            metric = "ptc",
            seurat_out = FALSE
        )
    )
})

test_that("clustify takes factor for metadata", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta$classified,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        verbose = TRUE
    )

    res2 <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta$classified,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        verbose = TRUE,
        exclude_genes = c("CD27", "ZNF232", "ZYX")
    )

    expect_true(res[1, 1] != res2[1, 1])
})

test_that("sce object clustifying", {
    res <- clustify(sce_small,
        cbmc_ref,
        cluster_col = "cell_type1",
        obj_out = FALSE
    )
    expect_true(nrow(res) == 13)
})

test_that("sce object clustify_lists", {
    other <- c("TAF12", "SNHG3")
    delta <- c("PCSK1", "LEPR")
    panm <- data.frame(other, delta)

    res <- clustify_lists(
        sce_small,
        marker = panm,
        cluster_col = "cell_type1",
        obj_out = FALSE,
        n = 100,
        metric = "pct"
    )
    expect_true(nrow(res) == 13)
})

test_that("clustify filters low cell number clusters", {
    pbmc_meta2 <- pbmc_meta %>% mutate(classified = as.character(classified))
    pbmc_meta2[1, "classified"] <- 15
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta2$classified,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        dr = "tsne",
        low_threshold_cell = 2,
        seurat_out = FALSE
    )
    expect_true(nrow(res) == 9)
})

test_that("clustify_lists filters low cell number clusters", {
    pbmc_meta2 <- pbmc_meta %>% mutate(classified = as.character(classified))
    pbmc_meta2[1, "classified"] <- 15
    res <- clustify_lists(
        input = pbmc_matrix_small,
        metadata = pbmc_meta2$classified,
        marker = cbmc_m,
        dr = "tsne",
        low_threshold_cell = 2,
        seurat_out = FALSE
    )
    expect_true(nrow(res) == 9)
})

test_that("clustify n_genes options limits number of variable genes", {
    res <- clustify(s_small3,
                    cbmc_ref,
                    cluster_col = "RNA_snn_res.1",
                    dr = "tsne",
                    obj_out = FALSE
    )
    res2 <- clustify(s_small3,
                     cbmc_ref,
                     n_genes = 2,
                     cluster_col = "RNA_snn_res.1",
                     dr = "tsne",
                     obj_out = FALSE
    )
    expect_true(res[1,1] != res2[1,1])
})

test_that("clustify n_genes options ignored if too large", {
    res <- clustify(s_small3,
                    cbmc_ref,
                    cluster_col = "RNA_snn_res.1",
                    dr = "tsne",
                    obj_out = FALSE
    )
    res2 <- clustify(s_small3,
                     cbmc_ref,
                     n_genes = 20,
                     cluster_col = "RNA_snn_res.1",
                     dr = "tsne",
                     obj_out = FALSE
    )
    expect_true(res[1,1] == res2[1,1])
})

# test_that("clustify finds SCT variable genes if needed", {
#     s <- CreateSeuratObject(pbmc_matrix_small)
#     s <- SCTransform(s)
#     res <- clustify(s,
#                     cbmc_ref,
#                     per_cell = TRUE,
#                     verbose = TRUE,
#                     dr = "tsne",
#                     obj_out = FALSE
#     )
#     expect_true(res[1,1] == res2[1,1])
# })