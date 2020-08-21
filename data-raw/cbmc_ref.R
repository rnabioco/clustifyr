library(usethis)
library(Seurat)
library(tidyverse)
library(clustifyr)

# following seurat tutorial from https://satijalab.org/seurat/v3.0/multimodal_vignette.html#identify-differentially-expressed-proteins-between-clusters
cbmc.rna <- as.sparse(read.csv(
    file = "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",",
    header = TRUE, row.names = 1
))
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)
cbmc.adt <- as.sparse(read.csv(
    file = "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",",
    header = TRUE, row.names = 1
))
cbmc.adt <- cbmc.adt[setdiff(rownames(x = cbmc.adt), c("CCR5", "CCR7", "CD10")), ]
cbmc <- CreateSeuratObject(counts = cbmc.rna)
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- FindNeighbors(cbmc, dims = 1:25)
cbmc <- FindClusters(cbmc, resolution = 0.8)
cbmc <- RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")
new.cluster.ids <- c(
    "Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono", "Mouse", "B",
    "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk",
    "Mouse", "DC", "pDCs"
)
names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)
cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)
cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
cbmc <- ScaleData(cbmc, assay = "ADT")
cbmc <- subset(cbmc, idents = c("Multiplets", "Mouse"), invert = TRUE)
DefaultAssay(cbmc) <- "ADT"
cbmc <- RunPCA(cbmc,
    features = rownames(cbmc), reduction.name = "pca_adt", reduction.key = "pca_adt_",
    verbose = FALSE
)
adt.data <- GetAssayData(cbmc, slot = "data")
adt.dist <- dist(t(adt.data))
cbmc[["rnaClusterID"]] <- Idents(cbmc)
cbmc[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
cbmc[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
cbmc <- FindClusters(cbmc, resolution = 0.2, graph.name = "adt_snn")
new.cluster.ids <- c(
    "CD4 T", "CD14+ Mono", "NK", "B", "CD8 T", "NK", "CD34+", "T/Mono doublets",
    "CD16+ Mono", "pDCs", "B"
)
names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)
cbmc[["citeID"]] <- Idents(cbmc)

m <- cbmc@meta.data %>%
    rownames_to_column("rn") %>%
    mutate(ID = ifelse(citeID != "CD8 T" & citeID != "CD4 T", as.character(rnaClusterID), as.character(citeID))) %>%
    mutate(ID = ifelse((rnaClusterID == "CD4 T" & citeID != "CD4 T") | (rnaClusterID == "CD8 T" & citeID != "CD8 T"), "Unknown", as.character(ID))) %>%
    column_to_rownames("rn")
cbmc@meta.data <- m

cbmc_refm <- use_seurat_comp(
    cbmc,
    cluster_col = "ID",
    var_genes_only = FALSE,
    assay_name = NULL
)

cbmc_ref <- as.data.frame(cbmc_refm) %>%
    as_tibble(rownames = "gene") %>%
    select(-Unknown, -`T/Mono doublets`) %>%
    filter(str_sub(gene, 1, 5) != "MOUSE" & str_sub(gene, 1, 5) != "ERCC_") %>%
    column_to_rownames("gene") %>%
    as.matrix()

cbmc_ref <- cbmc_ref[Matrix::rowSums(cbmc_ref) != 0, ]

var_genes <- ref_feature_select(cbmc_ref, n = 2000)
cbmc_ref <- cbmc_ref[var_genes, ]

usethis::use_data(cbmc_ref, compress = "xz", overwrite = TRUE)
