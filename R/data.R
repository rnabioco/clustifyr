#' Matrix of single-cell RNA-seq PBMCs.
#'
#' Count matrix of 3k pbmcs from Seurat3 tutorial, with only var.features
#'
#' @format A sparseMatrix with genes as rows and cells as columns.
#'
#' \describe{
#' }
#' R
#' @source (https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)
"pbmc_matrix_small"

#' Meta-data for single-cell RNA-seq PBMCs.
#'
#' Metadata, including umap, of 3k pbmcs from Seurat3 tutorial
#'
#' @format A data frame with 5 variables:
#' \describe{
#' }
#'
#' @source `[pbmc_matrix]` processed by Seurat
"pbmc_meta"

#' Marker genes identified by Seurat from single-cell RNA-seq PBMCs.
#'
#' Dataframe of markers from Seurat FindAllMarkers function
#'
#' @format A data frame with 7 variables:
#' \describe{
#' }
#'
#' @source `[pbmc_matrix]` processed by Seurat
"pbmc_markers"

#' Marker genes identified by M3Drop from single-cell RNA-seq PBMCs.
#'
#' Selected features of 3k pbmcs from Seurat3 tutorial
#'
#' @format A data frame with 3 variables:
#' \describe{
#' }
#'
#' @source `[pbmc_matrix]` processed by `[M3Drop]`
"pbmc_markers_M3Drop"

#' Variable genes identified by Seurat from single-cell RNA-seq PBMCs.
#'
#' Top 2000 variable genes from 3k pbmcs from Seurat3 tutorial
#'
#' @format A vector:
#' \describe{
#' }
#'
#' @source `[pbmc_matrix]` processed by Seurat
"pbmc_vargenes"

#' Bulk RNA-Seq data from sorted populations isolated from PBMCs.
#'
#'
#' @format A read count matrix with 14 variables:
#' \describe{
#' }
#'
#' @source (http://duffel.rail.bio/recount/v2/SRP051688/rse_gene.Rdata)
"pbmc_bulk_matrix"

#' Small clustered Seurat2 object
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  `[pbmc_small]` processed by seurat
"s_small"

#' Small clustered Seurat3 object
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  `[pbmc_small]` processed by Seurat
"s_small3"

#' reference matrix from seurat citeseq CBMC tutorial
#'
#' @format expression matrix
#' \describe{
#' }
#'
#' @source (https://satijalab.org/seurat/v3.0/multimodal_vignette.html#identify-differentially-expressed-proteins-between-clusters)
"cbmc_ref"

#' reference marker matrix from seurat citeseq CBMC tutorial
#'
#' @format matrix of markers
#' \describe{
#' }
#'
#' @source (https://satijalab.org/seurat/v3.0/multimodal_vignette.html#identify-differentially-expressed-proteins-between-clusters)
"cbmc_m"

#' lookup table for single cell object structures
#'
#' @format dataframe
#' \describe{
#' }
#'
#' @source  various packages
"object_loc_lookup"

#' Small clustered SCE object
#'
#' @format sce object
#' \describe{
#' }
#'
#' @source subsampled from pancreas data
"sce_small"
