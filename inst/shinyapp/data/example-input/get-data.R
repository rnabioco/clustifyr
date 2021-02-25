library(Seurat)
library(tidyverse)
library(here)

# data from lung injury single cell RNA-seq dataset (GSE113049)

mdata <- read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113049/suppl/GSE113049_cell_metadata.tsv.gz")
mat <- read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113049/suppl/GSE113049_count_matrix.tsv.gz",
                skip = 1,
                col_names = c("gene", mdata$cell))

mat <- column_to_rownames(mat, "gene") %>%
  as.matrix()

# downsample to 25 cells per class
set.seed(42)
cell_ids <- mdata %>%
  group_by(sample_type, cell_type) %>%
  sample_frac(0.05) %>%
  pull(cell)

mat <- mat[, cell_ids]

mdata <- mdata %>% filter(cell %in% cell_ids) %>%
  column_to_rownames("cell") %>%
  .[cell_ids, ]

var_genes <- CreateSeuratObject(mat, meta.data = mdata) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  VariableFeatures()

mat <- mat[var_genes, ]

mat <- as.data.frame(mat) %>%
  rownames_to_column("gene")
mdata <- mdata %>% rownames_to_column("cell")

write_csv(mat, here("data", "lung-data", "matrix.csv"))
write_csv(mdata, here("data", "lung-data", "meta-data.csv"))
