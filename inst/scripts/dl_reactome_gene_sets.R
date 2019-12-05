library(here)

dl_url <- "http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/c2.cp.reactome.v6.2.symbols.gmt"

proj_dir <- here()

download.file(dl_url,
              file.path(proj_dir,
                        "inst/extdata/c2.cp.reactome.v6.2.symbols.gmt"))
