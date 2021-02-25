library(shiny)
library(shinyjs)
library(shinyWidgets)
library(waiter)
library(dplyr)
library(readr)
library(tools)
library(clustifyr)
library(clustifyrdatahub)
library(rsconnect)
library(ExperimentHub)
library(Seurat)
library(shinydashboard)
#library(dashboardthemes)
library(tidyverse)
library(data.table)
library(R.utils)
library(DT)
library(GEOquery)
library(pheatmap)
library(googlesheets4)
library(openxlsx)

options(shiny.maxRequestSize = 1500 * 1024^2)
options(repos = BiocManager::repositories())
options(datatable.fread.datatable = FALSE)
options(shiny.reactlog = TRUE)
options(
  DT.options = list(
    dom = "tp",
    paging = TRUE,
    pageLength = 6,
    scrollX = TRUE,
    server = TRUE
  )
)

is_local <- Sys.getenv('SHINY_PORT') == ""

# google sheet
gs4_auth(cache = "sheet", email = TRUE, use_oob = TRUE)
sheetid <- "https://docs.google.com/spreadsheets/d/107qXuwo568wmPikaNDIe1supkGOYakwnRhiaTaL7U8c"

# setup experimenthub
eh <- ExperimentHub()
refs <- query(eh, "clustifyrdatahub")
refs_meta <- mcols(refs) %>% as.data.frame()
ref_dict <- refs$ah_id %>% setNames(refs$title)

# get clicked column
js <- c(
  "table.on('click', 'td', function(){",
  "  var cell = table.cell(this);",
  "  var colindex = cell.index().column;",
  "  var colname = table.column(colindex).header().innerText;",
  "  Shiny.setInputValue('column_clicked', colname);",
  "});"
)

# get active tab
js2 <- '
$(document).ready(function(){
  $("a[data-toggle=tab]").on("show.bs.tab", function(e){
    Shiny.setInputValue("activeTab", $(this).attr("data-value"));
  });
});
'

# GEO functions
make_button <- function(tbl){
  function(i){
    sprintf(
      paste0('<button id="button_%s_%d', '_', format(Sys.time(), "%H_%M_%S"), '" type="button" onclick="%s">Preview</button>'),
      tbl, i, "Shiny.setInputValue('button', this.id);")
  }
}

get_tar <- function(link) {
  paste0(
    str_remove(link, "/GSE[0-9]+_RAW.tar"),
    "/filelist.txt"
  )
  # str_c("https://ftp.ncbi.nlm.nih.gov/geo/series/", str_sub(id, 1, str_length(id) - 3), "nnn/", id, "/suppl/filelist.txt")
}

get_file_size <- function(url) {
  response <- httr::HEAD(url)
  size <- tryCatch(httr::headers(response)[["Content-Length"]] %>% as.numeric(),
                   error = "error_get")
  if (is.null(size)) {
    return("error_get")
  }
  utils:::format.object_size(size, "auto")
}

list_geo <- function(id) {
  message("fetching info for all files available...")
  # look for files
  out <- tryCatch(suppressMessages(GEOquery::getGEOSuppFiles(id,
                                                             makeDirectory = FALSE,
                                                             fetch_files = FALSE))$fname,
                  error = function(e) {
                    "error_get"
                  })
  # make links
  if (is.null(out)) {
    return("error_get")
  } 
  
  if (out == "error_get") {
    return("error_get")
  }

  out <- data.frame(file = out) %>%
    mutate(link = str_c("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE",
                        str_extract(file, "[0-9]{3}"),
                        "nnn/",
                        id,
                        "/suppl/",
                        file))
  out
}

prep_email <- function(id) {
  out <- tryCatch(
    suppressMessages(GEOquery::getGEO(
      GEO = id,
      filename = NULL,
      GSElimits = NULL, GSEMatrix = FALSE,
      AnnotGPL = FALSE, getGPL = FALSE,
      parseCharacteristics = FALSE
    )),
    error = function(e) {
      "error_get"
    })
  if (class(out) != "GSE") {
    return(out)
  } else {
    name <- paste0("Dr ", out@header$contact_name %>% str_remove(".+,"))
    if (id == "GSE113049") {
      email <- "placeholder@forexample.com"
    } else {
      email <- out@header$contact_email
    }

    link <- paste0("mailto:",
                   email,
                   "?subject=additional info request for ",
                   id,
                   "&body=Dear ",
                   name,
                   ",%0D%0A%0D%0ACan you please provide additional metadata information for the single cell dataset deposited on GEO, ",
                   id,
                   ".",
                   "%0D%0A%0D%0AThank you so much")
    return(link)
  }
}

preview_link <- function(link, n_row = 5, n_col = 50, verbose = TRUE) {
  # make sure link works
  message(link)
  if (!str_starts(str_to_lower(link), "http")) {
    return(NA)
  }

  # stream in a few lines only
  message("read")
  url1 <- url(link)
  if (str_ends(link, "\\.gz")) {
    temp <- readLines(gzcon(url1), n = n_row)
  } else {
    temp <- readLines(url1, n = n_row)
  }
  close(url1)
  readable <- map(temp, function(x) {all(charToRaw(x[1]) <= as.raw(127))}) %>%
    unlist() %>%
    all()
  if (!readable) {
    return(NULL)
  }

  # parsing, using fread auto
  temp_df <- tryCatch(data.table::fread(text = temp),#, header = TRUE, fill = TRUE),
                      error = function() {"parsing failed"})

  return(temp_df)
}

# load rdata to file name
load_rdata <- function(file) {
  env <- new.env()
  nm <- load(file, envir = env)[1]
  env[[nm]]
}

# Plot correlation heatmap
plot_hmap <- function (cor_mat,
                       col = clustifyr:::not_pretty_palette,
                       legend_title = NULL,
                       ...) {
  pheatmap::pheatmap(cor_mat,
                     color = colorRampPalette(col)(100),
                     ...)
}

# pull in someta
someta <- readRDS(url("https://github.com/rnabioco/someta/raw/master/inst/extdata/current_geo.rds"))
someta <- someta[ , ] %>% select(id, organism = org, usable, files = suppfiles, tar_files = tarfiles, geo, pubmed) %>%
  mutate(files = map_chr(files, function(x) paste0(x, collapse = "; "))) %>% 
  mutate(tar_files = map_chr(tar_files, function(x) paste0(x, collapse = "; "))) %>% 
  mutate(files = ifelse(str_length(files) > 0, files, "none")) %>% 
  mutate(tar_files = ifelse(str_length(tar_files) > 0 & tar_files != "error_parse", tar_files, "none")) %>% 
  mutate(usable = factor(usable, levels = c("yes", "no"))) %>% 
  mutate(summary = map_chr(geo, function(x) {
    g <- tryCatch(x$summary,
                  error = function(e) return(NULL))
    if (is.null(g)) {
      g <- "none"
    }
    if (length(g) > 0) {
      g <- str_c(g, collapse = " ")
    }
    g
  })) %>% 
  mutate(pubmed_id = map_chr(pubmed, function(x) {
    g <- tryCatch(x$pmid[1],
             error = function(e) return(NULL))
    if (is.null(g)) {
      g <- "none"
    }
    g
    })) %>% 
  mutate(pubmed_title = map_chr(pubmed, function(x) {
    g <- tryCatch(x$title[1],
                  error = function(e) return(NULL))
    if (is.null(g)) {
      g <- "none"
    }
    g
  }))



