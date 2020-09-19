
#' Build reference atlases from external UCSC cellbrowsers
#'
#' @param cb_url URL of cellbrowser dataset (e.g. http://cells.ucsc.edu/?ds=cortex-dev).
#' Note that the URL must contain the ds=dataset-name suffix.
#' @param cluster_col annotation field for summarizing gene expression (e.g. clustering,
#' cell-type name, samples, etc.)
#' @param ... additional args passed to average_clusters
#'
#' @importFrom httr http_error parse_url build_url
#'
#' @examples
#'
#'\donttest{
#' get_ext_reference(cb_url = "http://cells.ucsc.edu/?ds=kidney-atlas%2FFetal_Immune",
#'                  cluster_col = "celltype")
#' }
#' @export
get_ucsc_reference <- function(cb_url,
                               cluster_col,
                               ...){

  url <- httr::parse_url(cb_url)
  base_url <- url
  ds <- url$query$ds
  base_url$query <- ""

  mdata_url <- httr::modify_url(base_url,
                                path = file.path(ds, "meta.tsv"))
  if(!httr::http_error(mdata_url)){
    mdata <- readr::read_tsv(mdata_url)
  } else {
    stop("unable to find metadata at url: ", mdata_url)
  }

  mat_url <- httr::modify_url(base_url,
                              path = file.path(ds, "exprMatrix.tsv.gz"))
  if(!httr::http_error(mat_url)){
    mat <- readr::read_tsv(mat_url)
  } else {
    stop("unable to find matrix at url: ", mat_url)
  }

  mat <- tibble::column_to_rownames(mat, "gene")
  mat <- as.matrix(mat)
  average_clusters(mat, mdata, cluster_col = cluster_col, ...)
}

