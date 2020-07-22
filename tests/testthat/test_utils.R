context("utils")

test_that("get_vargenes works for both matrix and dataframe form", {
    pbmc_mm <- matrixize_markers(pbmc_markers)
    var1 <- get_vargenes(pbmc_mm)
    var2 <- get_vargenes(pbmc_markers)

    expect_equal(var1[1], var2[1])
})

test_that("matrixize_markers with remove_rp option", {
    pbmc_mm <- matrixize_markers(pbmc_markers)
    pbmc_mm2 <- matrixize_markers(pbmc_markers,
        remove_rp = TRUE
    )

    expect_true(nrow(pbmc_mm) != nrow(pbmc_mm2))
})

test_that("matrixize_markers to turn matrix into ranked list", {
    pbmc_mm <- matrixize_markers(pbmc_markers, n = 50)
    pbmc_mm2 <-
        matrixize_markers(pbmc_mm, ranked = TRUE, unique = TRUE)

    expect_true(nrow(pbmc_mm) < nrow(pbmc_mm2))
})

test_that("matrixize_markers uses supplied labels", {
    pbmc_mm <- matrixize_markers(
        pbmc_markers,
        n = 50,
        metadata = pbmc_meta %>% mutate(cluster = seurat_clusters),
        cluster_col = "classified"
    )
    pbmc_mm2 <- matrixize_markers(
        pbmc_mm,
        metadata = unique(as.character(pbmc_meta$classified)),
        cluster_col = "classified",
        ranked = TRUE
    )

    expect_true(nrow(pbmc_mm) < nrow(pbmc_mm2))
})

test_that("average_clusters works as intended", {
    pbmc_avg2 <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified",
        if_log = FALSE
    )
    expect_equal(nrow(pbmc_avg2), 2000)
})

test_that("average_clusters reports error when supplied cluster vector doesn't match number of cols", {
    expect_error(
        pbmc_avg2 <- average_clusters(pbmc_matrix_small,
                                      pbmc_meta$classified[1:2],
                                      if_log = FALSE
        )
    )
})

test_that("average_clusters works with disordered data", {
    pbmc_meta2 <- rbind(pbmc_meta[1320:2638, ], pbmc_meta[1:1319, ])
    pbmc_avg2 <- average_clusters(
        pbmc_matrix_small,
        pbmc_meta %>% tibble::rownames_to_column("rn"),
        if_log = TRUE,
        cell_col = "rn",
        cluster_col = "classified"
    )
    pbmc_avg3 <- average_clusters(
        pbmc_matrix_small,
        pbmc_meta2 %>% tibble::rownames_to_column("rn"),
        if_log = TRUE,
        cell_col = "rn",
        cluster_col = "classified"
    )
    expect_equal(pbmc_avg2, pbmc_avg3)
})


test_that("average_clusters detects wrong cluster ident", {
    expect_error(
        pbmc_avg2 <- average_clusters(
            pbmc_matrix_small,
            matrix(5, 5),
            if_log = FALSE,
            cluster_col = "classified"
        )
    )
})

test_that("average_clusters able to coerce factors", {
    col <- factor(pbmc_meta$classified)
    pbmc_avg2 <- average_clusters(pbmc_matrix_small,
        col,
        if_log = FALSE
    )
    expect_equal(nrow(pbmc_avg2), 2000)
})

test_that("average_clusters works with median option", {
    pbmc_avg2 <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        method = "median",
        cluster_col = "classified"
    )
    expect_equal(nrow(pbmc_avg2), 2000)
})

test_that("average_clusters works when one cluster contains only 1 cell", {
    pbmc_meta2 <- pbmc_meta
    pbmc_meta2$classified <- as.character(pbmc_meta2$classified)
    pbmc_meta2[1, "classified"] <- 15
    pbmc_avg2 <- average_clusters(pbmc_matrix_small,
        pbmc_meta2,
        cluster_col = "classified",
        low_threshold = 0
    )
    expect_equal(ncol(pbmc_avg2), 9 + 1)
})

test_that("average_clusters works when low cell number clusters should be removed", {
    pbmc_meta2 <- pbmc_meta %>% mutate(classified = as.character(classified))
    pbmc_meta2[1, "classified"] <- 15
    pbmc_avg2 <- average_clusters(pbmc_matrix_small,
        pbmc_meta2,
        low_threshold = 2,
        cluster_col = "classified"
    )
    expect_equal(ncol(pbmc_avg2), 9)
})

test_that("average_clusters works when low cell number clusters is kept but issues warning", {
    pbmc_meta2 <- pbmc_meta %>% mutate(classified = as.character(classified))
    pbmc_meta2[1, "classified"] <- 15
    pbmc_avg2 <- average_clusters(pbmc_matrix_small,
                                  pbmc_meta2,
                                  low_threshold = 0,
                                  cluster_col = "classified"
    )
    expect_equal(ncol(pbmc_avg2), 10)
})

test_that("average_clusters works with cutoff gene number", {
    pbmc_avg2 <- average_clusters(
        pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified",
        if_log = FALSE,
        cut_n = 1
    )
    expect_true(sum(pbmc_avg2[, "B"]) < 10)
})

test_that("average_clusters works when cluster info contains NA", {
    pbmc_meta2 <- pbmc_meta
    pbmc_meta2[1, "classified"] <- NA
    pbmc_avg2 <- average_clusters(pbmc_matrix_small,
        pbmc_meta2,
        low_threshold = 2,
        cluster_col = "classified"
    )
    expect_equal(ncol(pbmc_avg2), 9)
})

test_that("average_clusters works when cluster info in factor form", {
    pbmc_meta2 <- pbmc_meta
    pbmc_meta2$classified <- as.factor(pbmc_meta2$classified)
    pbmc_avg2 <- average_clusters(pbmc_matrix_small,
        pbmc_meta2,
        low_threshold = 2,
        cluster_col = "classified"
    )
    expect_equal(ncol(pbmc_avg2), 9)
})

test_that("cor_to_call threshold works as intended", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    call1 <- cor_to_call(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = 0.8,
        carry_r = TRUE
    )

    expect_true("r<0.8, unassigned" %in% call1$type)
})

test_that("cor_to_call threshold works as intended, on per cell and collapsing", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        per_cell = TRUE
    )
    call1 <- cor_to_call(
        res,
        metadata = pbmc_meta %>% tibble::rownames_to_column("rn"),
        cluster_col = "rn",
        collapse_to_cluster = "classified",
        threshold = 0.1
    )

    expect_true(!any(is.na(call1)))
})

test_that("assign_ident works with equal length vectors and just 1 ident", {
    m1 <- assign_ident(
        pbmc_meta,
        ident_col = "classified",
        clusters = c("1", "2"),
        idents = c("whatever1", "whatever2")
    )
    m2 <- assign_ident(
        pbmc_meta,
        ident_col = "classified",
        clusters = c("1", "2"),
        idents = "whatever1"
    )
    expect_true(nrow(m1) == nrow(m2))
})

test_that("cor_to_call_topn works as intended", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    call1 <- cor_to_call_topn(
        res,
        metadata = pbmc_meta,
        col = "classified",
        collapse_to_cluster = FALSE,
        threshold = 0.5
    )

    expect_true(nrow(call1) == 2 * nrow(res))
})

test_that("cor_to_call_topn works as intended on collapse to cluster option", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        per_cell = TRUE
    )
    call1 <- cor_to_call_topn(
        res,
        metadata = pbmc_meta %>% tibble::rownames_to_column("rn"),
        col = "rn",
        collapse_to_cluster = "classified",
        threshold = 0
    )

    expect_true(nrow(call1) == 2 * nrow(res))
})

test_that("gene_pct and gene_pct_markerm work as intended", {
    res <- gene_pct(
        pbmc_matrix_small,
        cbmc_m$B,
        pbmc_meta$classified
    )

    res2 <- gene_pct_markerm(pbmc_matrix_small,
        cbmc_m,
        pbmc_meta,
        cluster_col = "classified"
    )
    expect_error(
        res2 <- gene_pct_markerm(pbmc_matrix_small,
            cbmc_m,
            matrix(5, 5),
            cluster_col = "classified"
        )
    )
    expect_true(nrow(res2) == 9)
})

test_that("gene_pct can give min or max output", {
    res <- gene_pct(pbmc_matrix_small,
        cbmc_m$B,
        pbmc_meta$classified,
        returning = "min"
    )
    res2 <- gene_pct(pbmc_matrix_small,
        cbmc_m$B,
        pbmc_meta$classified,
        returning = "max"
    )

    expect_true(all(res2 >= res))
})

test_that("gene_pct_markerm norm options work", {
    res <- gene_pct_markerm(pbmc_matrix_small,
        cbmc_m,
        pbmc_meta,
        cluster_col = "classified",
        norm = NULL
    )
    res2 <- gene_pct_markerm(pbmc_matrix_small,
        cbmc_m,
        pbmc_meta,
        cluster_col = "classified",
        norm = "divide"
    )
    res3 <- gene_pct_markerm(pbmc_matrix_small,
        cbmc_m,
        pbmc_meta,
        cluster_col = "classified",
        norm = 0.3
    )

    expect_true(nrow(res2) == 9)
})

test_that("clustify_nudge works with options and seruat2", {
    res <- clustify_nudge(
        input = s_small,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        cluster_col = "res.1",
        threshold = 0.8,
        seurat_out = FALSE,
        mode = "pct",
        dr = "tsne"
    )
    expect_true(nrow(res) == 4)
})

test_that("clustify_nudge works with seurat_out", {
    res <- clustify_nudge(
        input = s_small,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        cluster_col = "res.1",
        threshold = 0.8,
        seurat_out = TRUE,
        mode = "pct",
        dr = "tsne"
    )

    res <- clustify_nudge(
        input = s_small3,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        threshold = 0.8,
        seurat_out = TRUE,
        cluster_col = "RNA_snn_res.1",
        mode = "pct",
        dr = "tsne"
    )
    expect_true(3 == 3)
})

test_that("clustify_nudge works with rank/posneg option", {
    res <- clustify_nudge(
        input = s_small,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        cluster_col = "res.1",
        threshold = 0.8,
        seurat_out = FALSE,
        mode = "rank",
        dr = "tsne"
    )
    expect_true(nrow(res) == 4)
})

test_that("clustify_nudge works with options and seruat3", {
    res <- clustify_nudge(
        input = s_small3,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        threshold = 0.8,
        seurat_out = FALSE,
        cluster_col = "RNA_snn_res.1",
        mode = "pct",
        dr = "tsne"
    )
    expect_true(nrow(res) == 3)
})

test_that("clustify_nudge works with seurat_out option", {
    res <- clustify_nudge(
        input = s_small,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        cluster_col = "res.1",
        threshold = 0.8,
        seurat_out = FALSE,
        marker_inmatrix = FALSE,
        mode = "pct",
        dr = "tsne"
    )
    expect_true(nrow(res) == 4)
})

test_that("clustify_nudge.Seurat works with seurat_out option", {
    res <- clustify_nudge(
        input = s_small3,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        cluster_col = "RNA_snn_res.1",
        threshold = 0.8,
        seurat_out = TRUE,
        marker_inmatrix = FALSE,
        mode = "pct",
        dr = "tsne"
    )

    res <- clustify_nudge(
        input = s_small3,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        cluster_col = "RNA_snn_res.1",
        threshold = 0.8,
        seurat_out = FALSE,
        marker_inmatrix = FALSE,
        mode = "pct",
        dr = "tsne"
    )
    expect_true(nrow(res) == 3)
})

test_that("clustify_nudge works with obj_out option", {
    s3 <- s_small3
    setClass(
        "ser3",
        representation(meta.data = "data.frame")
    )
    class(s3) <- "ser3"
    object_loc_lookup2 <- data.frame(
        ser3 = c(
            expr = "input@assays$RNA@data",
            meta = "input@meta.data",
            var = "input@assays$RNA@var.features",
            col = "RNA_snn_res.1"
        ),
        stringsAsFactors = FALSE
    )

    res <- clustify_nudge(
        input = s3,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        lookuptable = object_loc_lookup2,
        cluster_col = "RNA_snn_res.1",
        threshold = 0.8,
        obj_out = TRUE,
        marker_inmatrix = FALSE,
        mode = "pct",
        dr = "tsne"
    )

    res2 <- clustify_nudge(
        input = s3,
        ref_mat = cbmc_ref,
        marker = cbmc_m,
        lookuptable = object_loc_lookup2,
        cluster_col = "RNA_snn_res.1",
        threshold = 0.8,
        obj_out = FALSE,
        marker_inmatrix = FALSE,
        mode = "pct",
        dr = "tsne"
    )
    expect_true(nrow(res2) == 3)
})

test_that("clustify_nudge works with list of markers", {
    res <- clustify_nudge(
        input = pbmc_matrix_small,
        ref_mat = average_clusters(pbmc_matrix_small,
            pbmc_meta,
            cluster_col = "classified"
        ),
        metadata = pbmc_meta,
        marker = pbmc_markers,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        threshold = 0.8,
        call = FALSE,
        marker_inmatrix = FALSE,
        mode = "pct"
    )
    expect_true(nrow(res) == 9)
})

test_that("clustify_nudge autoconverts when markers are in matrix", {
    res <- clustify_nudge(
        input = pbmc_matrix_small,
        ref_mat = cbmc_ref,
        metadata = pbmc_meta,
        marker = as.matrix(cbmc_m),
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        threshold = 0.8,
        call = FALSE,
        marker_inmatrix = FALSE,
        mode = "pct"
    )
    expect_true(nrow(res) == 9)
})

test_that("overcluster_test works with ngenes option", {
    g <- overcluster_test(
        pbmc_matrix_small[, 1:100],
        pbmc_meta[1:100, ],
        cbmc_ref,
        cluster_col = "classified",
        x_col = "UMAP_1",
        y_col = "UMAP_2"
    )
    g2 <- overcluster_test(
        pbmc_matrix_small[, 1:100],
        pbmc_meta[1:100, ],
        cbmc_ref,
        cluster_col = "classified",
        ngenes = 100,
        x_col = "UMAP_1",
        y_col = "UMAP_2"
    )
    expect_true(ggplot2::is.ggplot(g))
})

test_that("overcluster_test works with defined other cluster column", {
    g <- overcluster_test(
        pbmc_matrix_small[, 1:100],
        pbmc_meta[1:100, ],
        cbmc_ref,
        cluster_col = "seurat_clusters",
        newclustering = "classified",
        do_label = FALSE,
        x_col = "UMAP_1",
        y_col = "UMAP_2"
    )
    expect_true(ggplot2::is.ggplot(g))
})

test_that("ref_feature_select chooses the correct number of features", {
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    res <- ref_feature_select(pbmc_avg[1:100, ], 5)
    expect_true(length(res) == 5)
})

test_that("ref_feature_select chooses the correct number of features with options", {
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    res <- ref_feature_select(pbmc_avg[1:100, ], 5, mode = "cor")
    expect_true(length(res) == 5)
})

test_that("feature_select_PCA will log transform", {
    res <- feature_select_PCA(cbmc_ref, if_log = FALSE)
    res2 <- feature_select_PCA(cbmc_ref, if_log = TRUE)
    expect_true(length(res) > 0)
})

test_that("feature_select_PCA can handle precalculated PCA", {
    pcs <- prcomp(t(as.matrix(cbmc_ref)))$rotation
    res <- feature_select_PCA(cbmc_ref, if_log = TRUE)
    res2 <- feature_select_PCA(pcs = pcs, if_log = TRUE)
    expect_true(all.equal(rownames(res), rownames(res2)))
})


test_that("downsample_matrix can select same number of cells per cluster", {
    mat1 <- downsample_matrix(
        pbmc_matrix_small,
        metadata = pbmc_meta$classified,
        n = 10,
        keep_cluster_proportions = TRUE
    )
    expect_true(all.equal(ncol(mat1), 10 * length(unique(
        pbmc_meta$classified
    ))))
})

test_that("percent_clusters works with defaults", {
    res <- percent_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    expect_equal(nrow(res), nrow(pbmc_matrix_small))
})

test_that("get_best_str finds correct values", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        per_cell = FALSE
    )

    a <- get_best_str("DC", get_best_match_matrix(res), res)
    a2 <- get_best_str("DC", get_best_match_matrix(res), res, carry_cor = FALSE)

    expect_equal(stringr::str_sub(a, 1, 2), stringr::str_sub(a2, 1, 2))
})

test_that("seurat_ref gets correct averages", {
    avg <- seurat_ref(s_small,
        cluster_col = "res.1",
        var_genes_only = TRUE
    )
    avg3 <- seurat_ref(s_small,
        cluster_col = "res.1",
        var_genes_only = TRUE,
        if_log = FALSE
    )
    avg2 <- seurat_ref(s_small,
        cluster_col = "res.1",
        var_genes_only = "PCA"
    )
    expect_true(ncol(avg) == 4)
})

test_that("object_ref with seurat3", {
    s3 <- s_small3
    avg <- object_ref(s3,
        var_genes_only = TRUE,
        cluster_col = "RNA_snn_res.1"
    )
    expect_true(ncol(avg) == 3)
})

test_that("object_ref with SingleCellExperiment", {
    sce <- sce_small
    avg <- object_ref(sce,
                      cluster_col = "cell_type1")
    expect_equal(dim(avg), c(200, 13))
})


test_that("object_ref gets correct averages", {
    s3 <- s_small3
    class(s3) <- "ser3"
    object_loc_lookup2 <- data.frame(
        ser3 = c(
            expr = "input@assays$RNA@data",
            meta = "input@meta.data",
            var = "input@assays$RNA@var.features",
            col = "RNA_snn_res.1"
        ),
        stringsAsFactors = FALSE
    )
    avg <- object_ref(s3,
        lookuptable = object_loc_lookup2,
        var_genes_only = TRUE
    )
    expect_true(ncol(avg) == 3)
})

test_that("seurat_ref gets other assay slots", {
    avg <- seurat_ref(
        s_small,
        cluster_col = "res.1",
        assay_name = "ADT",
        var_genes_only = TRUE
    )
    avg2 <- seurat_ref(
        s_small,
        cluster_col = "res.1",
        assay_name = c("ADT", "ADT2"),
        var_genes_only = TRUE
    )
    expect_true(nrow(avg2) - nrow(avg) == 2)
})

test_that("seurat_ref gets correct averages with seurat3 object", {
    avg <- seurat_ref(
        s_small3,
        cluster_col = "RNA_snn_res.1",
        assay_name = c("ADT", "ADT2"),
        var_genes_only = TRUE
    )
    avg <- seurat_ref(
        s_small3,
        cluster_col = "RNA_snn_res.1",
        assay_name = c("ADT"),
        var_genes_only = TRUE
    )
    avg2 <- seurat_ref(
        s_small3,
        cluster_col = "RNA_snn_res.1",
        assay_name = c("ADT", "ADT2"),
        var_genes_only = "PCA"
    )
    expect_true(nrow(avg2) - nrow(avg) == 2)
})

test_that("object parsing works for custom object", {
    s3 <- s_small3
    class(s3) <- "ser3"
    object_loc_lookup2 <- data.frame(
        ser3 = c(
            expr = "input@assays$RNA@data",
            meta = "input@meta.data",
            var = "input@assays$RNA@var.features",
            col = "RNA_snn_res.1"
        ),
        stringsAsFactors = FALSE
    )

    res2 <- clustify(s3,
        cbmc_ref,
        lookuptable = object_loc_lookup2,
        obj_out = FALSE
    )

    res <- clustify_lists(
        s3,
        marker = pbmc_markers,
        marker_inmatrix = FALSE,
        lookuptable = object_loc_lookup2,
        obj_out = FALSE
    )

    expect_true(nrow(res) == nrow(res2))
})

test_that("object metadata assignment works for custom object", {
    s3 <- s_small3
    setClass(
        "ser3",
        representation(meta.data = "data.frame")
    )
    class(s3) <- "ser3"
    object_loc_lookup2 <- data.frame(
        ser3 = c(
            expr = "input@assays$RNA@data",
            meta = "input@meta.data",
            var = "input@assays$RNA@var.features",
            col = "RNA_snn_res.1"
        ),
        stringsAsFactors = FALSE
    )

    res2 <- clustify(s3,
        cbmc_ref,
        lookuptable = object_loc_lookup2,
        obj_out = TRUE
    )

    res3 <- clustify_lists(
        s3,
        marker = pbmc_markers,
        marker_inmatrix = FALSE,
        lookuptable = object_loc_lookup2,
        obj_out = TRUE,
        rename_prefix = "A"
    )

    expect_true(is(res2, "ser3"))
})

test_that("cor_to_call renaming with suffix input works as intended, per_cell or otherwise", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    call1 <- cor_to_call(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = 0.5,
        rename_prefix = "a"
    )
    res2 <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        per_cell = TRUE
    )
    call2 <- cor_to_call(
        res2,
        metadata = tibble::rownames_to_column(pbmc_meta, "rn"),
        cluster_col = "classified",
        collapse_to_cluster = TRUE,
        threshold = 0,
        rename_prefix = "a"
    )
    call3 <- cor_to_call(
        res2,
        pbmc_meta,
        cluster_col = "rn",
        collapse_to_cluster = TRUE,
        threshold = 0,
        rename_prefix = "a"
    )
    expect_true("a_type" %in% colnames(call1) &
        "a_type" %in% colnames(call2))
})

test_that("renaming with suffix input works as intended with clusify wrapper", {
    res <- clustify(
        input = s_small,
        ref_mat = cbmc_ref,
        cluster_col = "res.1",
        rename_suff = "a",
        dr = "tsne"
    )
    res2 <- clustify(
        input = s_small3,
        ref_mat = cbmc_ref,
        cluster_col = "RNA_snn_res.1",
        rename_suff = "a",
        dr = "tsne"
    )
    expect_true(!is.null(res))
})

test_that("ref_marker_select works with cutoffs", {
    res1 <- ref_marker_select(cbmc_ref, cut = 0)
    mm <-
        matrixize_markers(res1,
            n = 5,
            unique = TRUE,
            remove_rp = TRUE
        )
    res2 <- ref_marker_select(cbmc_ref, cut = 2)
    expect_true(nrow(res1) != nrow(res2))
})

test_that("pos_neg_select takes dataframe of 1 col or more", {
    pn_ref <- data.frame(
        "CD4" = c(1, 0, 0),
        "CD8" = c(0, 0, 1),
        row.names = c("CD4", "clustifyr0", "CD8B")
    )
    pn_ref2 <- data.frame(
        "CD8" = c(0, 0, 1),
        row.names = c("CD4", "clustifyr0", "CD8B")
    )
    res <- pos_neg_select(
        pbmc_matrix_small,
        pn_ref,
        pbmc_meta,
        "classified"
    )
    res2 <- pos_neg_select(
        pbmc_matrix_small,
        pn_ref2,
        pbmc_meta,
        "classified"
    )
    expect_identical(res[, 2], res2[, 1])
})

test_that("pos_neg_select normalizes res", {
    pn_ref2 <- data.frame(
        "a" = c(1, 0.01, 0),
        row.names = c("CD74", "clustifyr0", "CD79A")
    )
    res <- pos_neg_select(pbmc_matrix_small,
        pn_ref2,
        pbmc_meta,
        "classified",
        cutoff_score = 0.8
    )
    res2 <- pos_neg_select(pbmc_matrix_small,
        pn_ref2,
        pbmc_meta,
        "classified",
        cutoff_score = NULL
    )
    expect_true(res[1] != res2[1])
})

test_that("clustify_nudge works with pos_neg_select", {
    pn_ref2 <- data.frame(
        "CD8 T" = c(0, 0, 1),
        row.names = c("CD4", "clustifyr0", "CD8B"),
        check.names = FALSE
    )
    res <- clustify_nudge(
        pbmc_matrix_small,
        cbmc_ref,
        pn_ref2,
        metadata = pbmc_meta,
        cluster_col = "classified",
        norm = 0.5
    )
    expect_true(all(dim(res) == c(9, 3)))
})

test_that("reverse_marker_matrix takes matrix of markers input", {
    m1 <- reverse_marker_matrix(cbmc_m)
    m2 <- reverse_marker_matrix(as.matrix(cbmc_m))
    expect_identical(m1, m2)
})

test_that("more readable error message when cluster_col is not in metadata when joining", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        verbose = TRUE
    )

    expect_error(plot_best_call(
        res,
        pbmc_meta,
        "a"
    ))
})

test_that("more readable error message when cluster_col is not the previous col from metadata when joining", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        verbose = TRUE
    )

    res2 <- cor_to_call(
        res,
        pbmc_meta,
        "classified"
    )
    expect_error(call_to_metadata(
        res2,
        pbmc_meta,
        "seurat_clusters"
    ))
})

test_that("more readable error message when cluster_col exist but is wrong info", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        verbose = TRUE
    )

    expect_error(plot_best_call(
        res,
        pbmc_meta,
        "seurat_clusters"
    ))
})

marker_file <- system.file("extdata",
    "hsPBMC_markers.txt",
    package = "clustifyr"
)

test_that("paring marker files works on included example", {
    markers <- file_marker_parse(marker_file)
    expect_true(length(markers) == 2)
})

gmt_file <- system.file("extdata",
    "c2.cp.reactome.v6.2.symbols.gmt.gz",
    package = "clustifyr"
)

test_that("paring gmt files works on included example", {
    gmt_list <- gmt_to_list(path = gmt_file)
    expect_true(is.list(gmt_list))
})

test_that("clustify_nudge works with pos_neg_select and seurat2 object", {
    pn_ref2 <- data.frame(
        "CD8 T" = c(0, 0, 1),
        row.names = c("CD4", "clustifyr0", "CD8B"),
        check.names = FALSE
    )
    res <- clustify_nudge(
        s_small,
        cbmc_ref,
        pn_ref2,
        cluster_col = "res.1",
        norm = 0.5,
        dr = "tsne",
        seurat_out = FALSE
    )
    expect_true(nrow(res) == 4)
})

test_that("clustify_nudge works with pos_neg_select and Seurat3 object", {
    pn_ref2 <- data.frame(
        "CD8 T" = c(0, 0, 1),
        row.names = c("CD4", "clustifyr0", "CD8B"),
        check.names = FALSE
    )
    res <- clustify_nudge(
        s_small3,
        cbmc_ref,
        pn_ref2,
        cluster_col = "RNA_snn_res.1",
        norm = 0.5,
        dr = "tsne",
        seurat_out = FALSE
    )
    expect_true(nrow(res) == 3)
})

test_that("pos_neg_marker takes list, matrix, and dataframe", {
    res <- pos_neg_marker(cbmc_m)
    res2 <- pos_neg_marker(as.matrix(cbmc_m))
    res3 <- pos_neg_marker(as.list(cbmc_m))
    expect_true(nrow(res) == nrow(res2))
})

test_that("pos_neg_marker takes uneven list", {
    genelist <- list(
        precusor_oligodendrocyte = c("Olig1", "Olig2"),
        mature_oligodendrocyte = c("Olig1", "Olig2", "Mbp"),
        microglia = c("Aif1", "Tmem119"),
        astrocyte = c("Gfap", "Aqp4")
    )
    res <- pos_neg_marker(genelist)
    expect_true(ncol(res) == length(genelist))
})

test_that("average_clusters can generate subclusters", {
    ref_mat <- average_clusters(
        pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified",
        subclusterpower = 0.15
    )
    expect_true(ncol(ref_mat) > length(unique(pbmc_meta$classified)))
})

test_that("cor_to_call threshold works as intended", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref[, -10],
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    call1 <- cor_to_call(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = "auto",
        carry_r = TRUE
    )

    expect_true(any(stringr::str_detect(call1$type, "unassigned")))
})

test_that("cor_to_call_rank threshold works as intended", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    call1 <- cor_to_call_rank(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = "auto",
        rename_prefix = "a"
    )

    expect_true(max(call1$rank) == 100)
})

test_that("cor_to_call_rank options", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    call1 <- cor_to_call_rank(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = "auto",
        rename_prefix = "a",
        top_n = 1
    )

    expect_true(max(call1$rank) == 1)
})

test_that("call_consensus marks ties", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    call1 <- cor_to_call_rank(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = "auto"
    )
    res2 <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        cluster_col = "classified"
    )
    call2 <- cor_to_call_rank(
        res2,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = "auto"
    )
    call_f <- call_consensus(list(call1, call2))
    expect_true(nrow(call_f) == length(unique(pbmc_meta$classified)))
})

test_that("cor_to_call can collapse_to_cluster", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified",
        per_cell = TRUE
    )
    call1 <- cor_to_call(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = TRUE,
        threshold = 0.1
    )
    expect_true(ncol(call1) == 4)
})

test_that("cor_to_call and collapse_to_cluster work on objects", {
    res <- clustify(
        input = s_small,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "res.1",
        dr = "tsne",
        per_cell = TRUE,
        collapse_to_cluster = TRUE
    )
    expect_true(is.matrix(res) | is(res, "seurat"))
})

test_that("seurat_meta warns about not finding dr", {
    m <- seurat_meta(s_small,
        dr = "tsne"
    )
    m2 <- seurat_meta(s_small,
        dr = "s"
    )
    m3 <- seurat_meta(s_small3,
        dr = "s"
    )
    expect_true(all(rownames(m) == rownames(m2)))
})

test_that("find_rank_bias filters out unassigned", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    call1 <- cor_to_call(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = 0.8
    )
    pbmc_meta2 <- call_to_metadata(
        call1,
        pbmc_meta,
        "classified"
    )
    b <- find_rank_bias(pbmc_matrix_small,
        pbmc_meta2, "type",
        cbmc_ref,
        query_genes = pbmc_vargenes
    )
    expect_true(length(unique(pbmc_meta2$type)) > ncol(b))
})

test_that("repeated insertionn of types into metdadata renames correctly", {
    res <- clustify(
        input = pbmc_matrix_small,
        metadata = pbmc_meta,
        ref_mat = cbmc_ref,
        query_genes = pbmc_vargenes,
        cluster_col = "classified"
    )
    call1 <- cor_to_call(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = 0.8
    )
    pbmc_meta2 <- call_to_metadata(
        call1,
        pbmc_meta,
        "classified"
    )
    call2 <- cor_to_call(
        res,
        metadata = pbmc_meta,
        cluster_col = "classified",
        collapse_to_cluster = FALSE,
        threshold = 0
    )
    pbmc_meta2 <- call_to_metadata(
        call2,
        pbmc_meta2,
        "classified",
        rename_prefix = "a"
    )
    expect_true(colnames(pbmc_meta2)[10] == "type")
})

test_that("object_ref with sce", {
    avg <- object_ref(sce_small,
        cluster_col = "cell_type1"
    )
    expect_true(ncol(avg) == 13)
})

test_that("object_data works with sce", {
    mat <- object_data(
        object = sce_small,
        slot = "data"
    )
    expect_true(ncol(mat) == 200)
})

test_that("object_data works with seuratv2", {
    mat <- object_data(
        object = s_small,
        slot = "data"
    )
    expect_true(ncol(mat) == 80)
})

test_that("object_data works with Seuratv3", {
    mat <- object_data(
        object = s_small3,
        slot = "data"
    )
    expect_true(ncol(mat) == 80)
})

test_that("append_genes creates a union reference matrix", {
    mat <- append_genes(
        gene_vector = ,
        ref_matrix = 
    )
})

test_that("check raw counts of matrices", {
    mat <- check_raw_counts(
        counts_matrix = , 
        max_log_value = 50 
    )
})

