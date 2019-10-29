s <- readRDS("segerstolpe.rds")
sce_small <- s[1:200,1:200]

attr(attr(sce_small@int_elementMetadata, "class"), "package") <- NULL
attr(attr(sce_small@int_colData, "class"), "package") <- NULL
attr(attr(sce_small@elementMetadata, "class"), "package") <- NULL
attr(attr(sce_small@reducedDims, "class"), "package") <- NULL
attr(attr(sce_small@rowRanges, "class"), "package") <- NULL
attr(attr(sce_small@colData, "class"), "package") <- NULL
attr(attr(sce_small, "class"), "package") <- NULL
attr(attr(sce_small@assays, "class"), "package") <- NULL
attr(attr(sce_small@rowRanges@partitioning, "class"), "package") <- NULL
attr(attr(sce_small@rowRanges@elementMetadata, "class"), "package") <- NULL
attr(attr(sce_small@rowRanges@unlistData, "class"), "package") <- NULL
attr(attr(sce_small@rowRanges@unlistData@elementMetadata, "class"), "package") <- NULL
attr(attr(sce_small@rowRanges@unlistData@seqnames, "class"), "package") <- NULL
attr(attr(sce_small@rowRanges@unlistData@ranges, "class"), "package") <- NULL
attr(attr(sce_small@rowRanges@unlistData@strand, "class"), "package") <- NULL
attr(attr(sce_small@rowRanges@unlistData@seqinfo, "class"), "package") <- NULL
attr(attr(sce_small@assays@.xData[["data"]], "class"), "package") <- NULL
attr(attr(sce_small@assays$`.->data`, "class"), "package") <- NULL
attr(attr(sce_small@assays$data, "class"), "package") <- NULL
sce_small@assays$.objectPackage <- NULL

usethis::use_data(sce_small, compress = "xz", overwrite = TRUE)
