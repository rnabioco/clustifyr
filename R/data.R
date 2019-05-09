#' Matrix of single-cell RNA-seq PBMCs.
#'
#' A sample of 300 PBMC cells, downsampled from a bigger experiment, 30 cells per cluster
#'
#' @format A sparseMatrix with genes as rows and cells as columns.
#'
#' \describe{
#' }
#' R
#' @source (https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k)
"pbmc4k_matrix"

#' Meta-data for single-cell RNA-seq PBMCs.
#'
#' A sample of 300 PBMC cells, downsampled from a bigger experiment, 30 cells per cluster
#'
#' @format A data frame with 5 variables:
#' \describe{
#' }
#'
#' @source `[pbmc4k_matrix]` processed by `[seurat]`
"pbmc4k_meta"

#' Marker genes identified by Seurat from single-cell RNA-seq PBMCs.
#'
#' A sample of 4000 cells from a bigger experiment
#'
#' @format A data frame with 7 variables:
#' \describe{
#' }
#'
#' @source `[pbmc4k_matrix]` processed by `[seurat]`
"pbmc4k_markers"

#' Marker genes identified by M3Drop from single-cell RNA-seq PBMCs.
#'
#' A sample of 4000 cells from a bigger experiment
#'
#' @format A data frame with 3 variables:
#' \describe{
#' }
#'
#' @source `[pbmc4k_matrix]` processed by `[M3Drop]`
"pbmc4k_markers_M3Drop"

#' variable genes identified by Seurat from single-cell RNA-seq PBMCs.
#'
#' A sample of 4000 cells from a bigger experiment
#'
#' @format A vector:
#' \describe{
#' }
#'
#' @source `[pbmc4k_matrix]` processed by `[seurat]`
"pbmc4k_vargenes"

#' Bulk RNA-Seq data from sorted populations isolated from PBMCs.
#'
#'
#' @format A read count matrix with 14 variables:
#' \describe{
#' }
#'
#' @source (http://duffel.rail.bio/recount/v2/SRP051688/rse_gene.Rdata)
"pbmc_bulk_matrix"

#' Variable genes defined by Principal components analysis of pbmc4k data
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  `[pbmc4k_matrix]` processed by `[seurat]`
"pbmc_pca"

#' Small clustered Seurat2 object
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  `[pbmc_small]` processed by `[seurat]`
"s_small"

#' Small clustered Seurat3 object
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  `[pbmc_small]` processed by `[Seurat]`
"s_small3"

#' Small clustered FSCE object
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  `[fsce_small]` processed by `[Scrunchy]`
"fsce_small"

#' reference matrix from seurat citeseq CBMC
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  `[seurat]` tutorial, averaged per cluster
"cbmc_ref"

#' reference marker matrix from seurat citeseq CBMC
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  `[seurat]` tutorial, matrixized top markers
"cbmc_m"

#' lookup table for single cell object structures
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  various packages
"object_loc_lookup"
