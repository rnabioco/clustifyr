library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(recount)

dl_recount <- function(sra_id){
  download_study(sra_id)
  load(file.path(sra_id, "rse_gene.Rdata"))
  # no longer need to downloaded data
  unlink(sra_id, recursive = TRUE)
  
  rse <- scale_counts(rse_gene)
  read_counts <- assay(rse, "counts")
  gene_ids <- rownames(read_counts)
  # get gene symbols, which are stored in rowData
  id2symbol <- data_frame(ids = rowData(rse_gene)$gene_id,
                          symbols = rowData(rse_gene)$symbol@listData) %>%
    mutate(symbols = map_chr(symbols, ~.x[1]))
  
  # clean up metadata into a dataframe
  mdata <- colData(rse)
  mdata_cols <- lapply(mdata$characteristics,
                       function(x){str_match(x, "^([^:]+):")[, 2]}) %>%
    unique() %>%
    unlist()
  
  mdata <- data_frame(run =  mdata$run,
                      all_data = as.list(mdata$characteristics)) %>%
    mutate(out = purrr::map_chr(all_data,
                                ~str_c(.x, collapse = "::"))) %>%
    tidyr::separate(out,
                    sep = "::",
                    into = mdata_cols) %>%
    select(-all_data) %>%
    mutate_at(.vars = vars(-matches("run")),
              .funs = function(x) str_match(x, ": (.+)")[, 2])
  
  # convert ids to symbols
  row_ids_to_symbols <- left_join(data_frame(ids = gene_ids),
                                  id2symbol, by = "ids")
  
  if(length(gene_ids) != nrow(row_ids_to_symbols)) {
    warning("gene id mapping to symbols produce more or less ids")
  }
  
  row_ids_to_symbols <- filter(row_ids_to_symbols, !is.na(symbols))
  
  out_df <- read_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(., row_ids_to_symbols,
              by = c("gene_id" = "ids")) %>%
    dplyr::select(-gene_id) %>%
    dplyr::select(symbols, everything()) %>%
    filter(!is.na(symbols))
  
  out_matrix <- tidyr::gather(out_df, library, expr, -symbols) %>%
    group_by(symbols, library) %>%
    summarize(expr = sum(expr)) %>%
    tidyr::spread(library, expr) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("symbols") %>%
    as.matrix()
  
  list(read_counts = out_matrix,
       meta_data = mdata)
}

# download
pbmc_data <- dl_recount("SRP051688")
# filter
good_libs <- filter(pbmc_data$meta_data, str_detect(time, "0"))
pbmc_data <- pbmc_data$read_counts[, good_libs$run]
# rename
new_ids <- left_join(data_frame(run = colnames(pbmc_data)), 
                     good_libs, by = "run") %>% 
  group_by(`cell type`) %>% 
  mutate(cell_id = stringr::str_c(`cell type`, " rep ", row_number())) %>% 
  pull(cell_id)

colnames(pbmc_data) <- new_ids
pbmc_bulk_matrix <- pbmc_data
usethis::use_data(pbmc_bulk_matrix, compress = "xz", overwrite = TRUE)
