library(here)

dl_url <- "http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/c2.cp.reactome.v6.2.symbols.gmt"

proj_dir <- here()

download.file(dl_url,
              file.path(proj_dir,
                        "inst/extdata/c2.cp.reactome.v6.2.symbols.gmt"))

R.utils::gzip(file.path(proj_dir,
                        "inst/extdata/c2.cp.reactome.v6.2.symbols.gmt"))
# hsPBMC_markers.txt taken from garnett website

dl_url <- "https://cole-trapnell-lab.github.io/garnett/marker_files/hsPBMC_markers.txt"

download.file(dl_url,
              file.path(proj_dir,
                        "inst/extdata/c2.cp.reactome.v6.2.symbols.gmt"))
