library(tidyverse)
# read 10x gene.tsv file for all genes
human_genes_10x <- read_tsv("genes.tsv", col_names = FALSE) %>% 
  pull(X2)

usethis::use_data(human_genes_10x, compress = "xz", overwrite = TRUE)
