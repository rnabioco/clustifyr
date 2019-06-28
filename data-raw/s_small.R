library(usethis)

# need seurat v2
s_small <- Seurat::pbmc_small
attr(attr(s_small, "class"), "package") <- NULL
attr(attr(attr(s, "dr")[["pca"]], "class"), "package") <- NULL
attr(attr(attr(s, "dr")[["tsne"]], "class"), "package") <- NULL
attr(attr(attr(s, "spatial"), "class"), "package") <- NULL
attr(attr(attr(s@dr$pca, "jackstraw"), "class"), "package") <- NULL

usethis::use_data(s_small, compress = "xz", overwrite = TRUE)