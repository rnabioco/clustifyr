library(tidyverse)

d1 <- c(
    paste0("[", "ref_MCA", "](", "https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_MCA.rda", ")"),
    "Mouse Cell Atlas",
    dim(clustifyrdata::ref_MCA)[2],
    dim(clustifyrdata::ref_MCA)[1],
    "mouse",
    paste0("[from](", "https://www.cell.com/cell/fulltext/S0092-8674(18)30116-8", ")")
)

d2 <- c(
    paste0("[", "ref_tabula_muris_drop", "](", "https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_tabula_muris_drop.rda", ")"),
    "Tabula Muris (10X)",
    dim(clustifyrdata::ref_tabula_muris_drop)[2],
    dim(clustifyrdata::ref_tabula_muris_drop)[1],
    "mouse",
    paste0("[from](", "https://www.nature.com/articles/s41586-018-0590-4", ")")
)

d3 <- c(
    paste0("[", "ref_tabula_muris_facs", "](", "https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_tabula_muris_facs.rda", ")"),
    "Tabula Muris (SmartSeq2)",
    dim(clustifyrdata::ref_tabula_muris_facs)[2],
    dim(clustifyrdata::ref_tabula_muris_facs)[1],
    "mouse",
    paste0("[from](", "https://www.nature.com/articles/s41586-018-0590-4", ")")
)

d4 <- c(
    paste0("[", "ref_mouse.rnaseq", "](", "https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_mouse.rnaseq.rda", ")"),
    "Mouse RNA-seq from 28 cell types",
    dim(clustifyrdata::ref_mouse.rnaseq)[2],
    dim(clustifyrdata::ref_mouse.rnaseq)[1],
    "mouse",
    paste0("[from](", "https://genome.cshlp.org/content/early/2019/03/11/gr.240093.118", ")")
)

d5 <- c(
    paste0("[", "ref_moca_main", "](", "https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_moca_main.rda", ")"),
    "Mouse Organogenesis Cell Atlas (main cell types)",
    dim(clustifyrdata::ref_moca_main)[2],
    dim(clustifyrdata::ref_moca_main)[1],
    "mouse",
    paste0("[from](", "https://www.nature.com/articles/s41586-019-0969-x", ")")
)

d6 <- c(
    paste0("[", "ref_hema_microarray", "](", "https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_hema_microarray.rda", ")"),
    "Human hematopoietic cell microarray",
    dim(clustifyrdata::ref_hema_microarray)[2],
    dim(clustifyrdata::ref_hema_microarray)[1],
    "human",
    paste0("[from](", "https://www.cell.com/fulltext/S0092-8674(11)00005-5", ")")
)

d7 <- c(
    paste0("[", "ref_cortex_dev", "](", "https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_cortex_dev.rda", ")"),
    "Human cortex development scRNA-seq",
    dim(clustifyrdata::ref_cortex_dev)[2],
    dim(clustifyrdata::ref_cortex_dev)[1],
    "human",
    paste0("[from](", "https://science.sciencemag.org/content/358/6368/1318.long", ")")
)

d8 <- c(
    paste0("[", "ref_pan_indrop", "](", "https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_pan_indrop.rda", ")"),
    "Human pancreatic cell scRNA-seq (inDrop)",
    dim(clustifyrdata::ref_pan_indrop)[2],
    dim(clustifyrdata::ref_pan_indrop)[1],
    "human",
    paste0("[from](", "https://www.cell.com/fulltext/S2405-4712(16)30266-6", ")")
)

d9 <- c(
    paste0("[", "ref_pan_smartseq2", "](", "https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_pan_smartseq2.rda", ")"),
    "Human pancreatic cell scRNA-seq (SmartSeq2)",
    dim(clustifyrdata::ref_pan_smartseq2)[2],
    dim(clustifyrdata::ref_pan_smartseq2)[1],
    "human",
    paste0("[from](", "https://www.sciencedirect.com/science/article/pii/S1550413116304363", ")")
)

d <- data.frame(d1, d2, d3, d4, d5, d6, d7, d8, d9) %>% t()
colnames(d) <- c("name", "desc", "ntypes", "ngenes", "org", "from_pub")
downrefs <- d %>% as.tibble()
usethis::use_data(downrefs, compress = "xz", overwrite = TRUE)
