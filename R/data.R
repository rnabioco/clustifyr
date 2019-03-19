#' Matrix of single-cell RNA-seq PBMCs.
#'
#' A sample of 4000 PBMC cells from a bigger experiment
#'
#' @format A sparseMatrix with genes as rows and cells as columns.
#'
#' \describe{
#' }
#'R
#' @source (https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k)
"pbmc4k_matrix"

#' Matrix of average gene expression per cluster from single-cell RNA-seq PBMCs.
#'
#' A sample of 4000 PBMC cells from a bigger experiment
#'
#' @format A sparseMatrix with genes as rows and clusters as columns.
#'
#' \describe{
#' }
#'
#' @source `[pbmc4k_matrix]` averaged per cluster
"pbmc4k_avg"

#' Meta-data for single-cell RNA-seq PBMCs.
#'
#' A sample of 4000 cells from a bigger experiment
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

#' Matrix of single-cell RNA-seq PBMCs from 5'end kit.
#'
#' A sample of 8000 cells from a bigger experiment
#'
#' @format A sparseMatrix with genes as rows and cells as columns.
#' \describe{
#' }
#'
#' @source (https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_pbmc_5gex)
"pbmc5_matrix"

#' Meta-data for single-cell RNA-seq PBMCs from 5'end kit.
#'
#' A sample of 8000 cells from a bigger experiment
#'
#' @format A data frame with 4 variables:
#' \describe{
#' }
#'
#' @source `[pbmc5_matrix]` processed by `[seurat]`
"pbmc5_meta"

#' Marker genes identified by Seurat from single-cell RNA-seq PBMCs from 5'end kit.
#'
#' A sample of 8000 cells from a bigger experiment
#'
#' @format A data frame with 7 variables:
#' \describe{
#' }
#'
#' @source `[pbmc5_matrix]` processed by `[seurat]`
"pbmc5_markers"

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

#' Small clustered Seurat object
#'
#' @format XXX
#' \describe{
#' }
#'
#' @source  `[pbmc_small]` processed by `[seurat]`
"s_small"

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
