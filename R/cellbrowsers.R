
#' Build reference atlases from external UCSC cellbrowsers
#'
#' @param cb_url URL of cellbrowser dataset (e.g. http://cells.ucsc.edu/?ds=cortex-dev).
#' Note that the URL must contain the ds=dataset-name suffix.
#' @param cluster_col annotation field for summarizing gene expression (e.g. clustering,
#' cell-type name, samples, etc.)
#' @param ... additional args passed to average_clusters
#'
#' @importFrom httr http_error parse_url build_url
#' @return reference matrix
#' @examples
#' \dontrun{
#' 
#' # many datasets hosted by UCSC have UMI counts in the expression matrix
#' # set if_log = FALSE if the expression matrix has not been natural log transformed
#' 
#' get_ucsc_reference(cb_url = "https://cells.ucsc.edu/?ds=evocell+mus-musculus+marrow",
#'                    cluster_col = "Clusters", if_log = FALSE)
#'                    
#' get_ucsc_reference(cb_url = "http://cells.ucsc.edu/?ds=muscle-cell-atlas",
#'                    cluster_col = "cell_annotation",
#'                    if_log = FALSE)
#' }
#' @export
get_ucsc_reference <- function(cb_url,
                               cluster_col,
                               ...){
  if(!requireNamespace("R.utils", quietly = TRUE)) {
    stop("This function requires the R.utils package, please install\n",
         "install.packages('R.utils')")
  }
  
  if(!requireNamespace("data.table", quietly = TRUE)) {
    stop("This function requires the data.table package, please install\n",
         "install.packages('data.table')")
  }
  
  url <- httr::parse_url(cb_url)
  base_url <- url
  ds <- url$query$ds
  
  # ds can include sub-datasets with syntax, "dataset+subdataset+and-so-on"
  # files are hosted at urls: dataset/subdataset/andsoon/..."
  ds_split <- strsplit(ds, "+", fixed = TRUE)[[1]]
  ds <- paste0(ds_split, collapse = "/")
  base_url$query <- ""

  mdata_url <- httr::modify_url(base_url,
                                path = file.path(ds, "meta.tsv"))
  if(!httr::http_error(mdata_url)){
    mdata <- data.table::fread(mdata_url, data.table = FALSE, sep = "\t")
  } else {
    stop("unable to find metadata at url: ", mdata_url)
  }
  
  mat_url <- httr::modify_url(base_url,
                              path = file.path(ds, "exprMatrix.tsv.gz"))
  if(!httr::http_error(mat_url)){
    mat <- data.table::fread(mat_url, data.table = FALSE, sep = "\t")
  } else {
    stop("unable to find matrix at url: ", mat_url)
  }

  rownames(mat) <- mat[, 1]
  mat[, 1] <- NULL
  mat <- as.matrix(mat)
  
  mm <- max(mat)
  
  if(mm > 50) {
    dots <- list(...)
    if(!"if_log" %in% names(dots) || dots$if_log) {
      warning("the data matrix has a maximum value of ",   mm, "\n",
              "the data are likely not log transformed,\n",
              "please set the if_log argument for average clusters accordingly")
    }
  }
  
  average_clusters(mat, mdata, cluster_col = cluster_col, ...)
}

