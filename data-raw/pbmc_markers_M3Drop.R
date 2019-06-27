library(Seurat)
library(tidyverse)
library(M3Drop)
library(usethis)

# follow seurat tutorial from https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
pbmc.data <- Read10X(data.dir = "/Users/rf/Downloads/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

tm <- expm1(as.matrix(pbmc@assays$RNA@data))
Normalized_data <- M3DropCleanData(tm,
                                   labels = rownames(pbmc@assays$RNA@data),
                                   is.counts = FALSE)
fits <- M3DropDropoutModels(Normalized_data$data)
pbmc_markers_M3Drop <- M3DropDifferentialExpression(Normalized_data$data,
                                         mt_method = "fdr",
                                         mt_threshold = 0.05)
usethis::use_data(pbmc_markers_M3Drop, compress = "xz", overwrite = TRUE)