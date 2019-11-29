library(usethis)
library(SingleCellExperiment)

s <- readRDS(url("https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/segerstolpe.rds"))
sce_small <- s[1:200, 1:200]

use_data(sce_small, compress = "xz", overwrite = TRUE)
