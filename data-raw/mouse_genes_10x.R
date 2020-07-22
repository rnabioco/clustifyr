library(tidyverse)
# read 10x gene.tsv file for all genes
mouse_genes_10x <- read_tsv("genes.tsv", col_names = FALSE) %>% 
  pull(X2)

usethis::use_data(mouse_genes_10x, compress = "xz", overwrite = TRUE)
