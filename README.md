
# clustifyr

<!-- badges: start -->

[![R-CMD-check-bioc](https://github.com/rnabioco/clustifyr/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/rnabioco/clustifyr/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/rnabioco/clustifyr/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/rnabioco/clustifyr?branch=devel)
[![platforms](https://bioconductor.org/shields/availability/release/clustifyr.svg)](https://bioconductor.org/packages/release/bioc/html/clustifyr.html)
[![bioc](https://bioconductor.org/shields/years-in-bioc/clustifyr.svg)](https://bioconductor.org/packages/release/bioc/html/clustifyr.html)
[![\#downloads](https://img.shields.io/badge/%23%20downloads-8045-brightgreen)](https://bioconductor.org/packages/stats/bioc/clustifyr/clustifyr_stats.tab)
<!-- badges: end -->

clustifyr classifies cells and clusters in single-cell RNA sequencing
experiments using reference bulk RNA-seq data sets, sorted microarray
expression data, single-cell gene signatures, or lists of marker genes.

## Installation

Install the Bioconductor version with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clustifyr")
```

Install the development version with:

``` r
BiocManager::install("rnabioco/clustifyr")
```

## Example usage

In this example we use the following built-in input data:

- an expression matrix of single cell RNA-seq data (`pbmc_matrix_small`)
- a metadata data.frame (`pbmc_meta`), with cluster information stored
  (`"classified"`)
- a vector of variable genes (`pbmc_vargenes`)
- a matrix of mean normalized scRNA-seq UMI counts by cell type
  (`cbmc_ref`)

We then calculate correlation coefficients and plot them on a
pre-calculated projection (stored in `pbmc_meta`).

``` r
library(clustifyr)

# calculate correlation
res <- clustify(
  input = pbmc_matrix_small,
  metadata = pbmc_meta$classified,
  ref_mat = cbmc_ref,
  query_genes = pbmc_vargenes
)

# print assignments
cor_to_call(res)
#> # A tibble: 9 × 3
#> # Groups:   cluster [9]
#>   cluster      type           r
#>   <chr>        <chr>      <dbl>
#> 1 B            B          0.909
#> 2 CD14+ Mono   CD14+ Mono 0.915
#> 3 FCGR3A+ Mono CD16+ Mono 0.929
#> 4 Memory CD4 T CD4 T      0.861
#> 5 Naive CD4 T  CD4 T      0.889
#> 6 DC           DC         0.849
#> 7 Platelet     Mk         0.732
#> 8 CD8 T        NK         0.826
#> 9 NK           NK         0.894

# plot assignments on a projection
plot_best_call(
  cor_mat = res,
  metadata = pbmc_meta,
  cluster_col = "classified"
)
```

![](man/figures/readme_example-1.png)<!-- -->

`clustify()` can take a clustered `SingleCellExperiment` or `seurat`
object (both v2 and v3) and assign identities.

``` r
# for SingleCellExperiment
clustify(
  input = sce_small,          # an SCE object
  ref_mat = cbmc_ref,         # matrix of RNA-seq expression data for each cell type
  cluster_col = "cell_type1", # name of column in meta.data containing cell clusters
  obj_out = TRUE              # output SCE object with cell type inserted as "type" column
) 
#> class: SingleCellExperiment 
#> dim: 200 200 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(200): SGIP1 AZIN2 ... TAF12 SNHG3
#> rowData names(10): feature_symbol is_feature_control ... total_counts
#>   log10_total_counts
#> colnames(200): AZ_A1 AZ_A10 ... HP1502401_E18 HP1502401_E19
#> colData names(35): cell_quality cell_type1 ... type r
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):

library(Seurat)
# for Seurat3/4
clustify(
  input = s_small3,
  cluster_col = "RNA_snn_res.1",
  ref_mat = cbmc_ref,
  seurat_out = TRUE
)
#> An object of class Seurat 
#> 230 features across 80 samples within 1 assay 
#> Active assay: RNA (230 features, 20 variable features)
#>  2 dimensional reductions calculated: pca, tsne

# New output option, directly as a vector (in the order of the metadata), which can then be inserted into metadata dataframes and other workflows
clustify(
  input = s_small3,
  cluster_col = "RNA_snn_res.1",
  ref_mat = cbmc_ref,
  vec_out = TRUE
)
#>  [1] "Mk"         "Mk"         "Mk"         "Mk"         "Mk"        
#>  [6] "Mk"         "Mk"         "Mk"         "Mk"         "Mk"        
#> [11] "B"          "B"          "B"          "B"          "B"         
#> [16] "B"          "B"          "B"          "B"          "B"         
#> [21] "CD16+ Mono" "CD16+ Mono" "CD16+ Mono" "CD16+ Mono" "CD16+ Mono"
#> [26] "CD16+ Mono" "CD16+ Mono" "CD16+ Mono" "CD16+ Mono" "CD16+ Mono"
#> [31] "Mk"         "B"          "Mk"         "Mk"         "Mk"        
#> [36] "Mk"         "Mk"         "Mk"         "Mk"         "Mk"        
#> [41] "Mk"         "B"          "Mk"         "Mk"         "B"         
#> [46] "B"          "Mk"         "Mk"         "Mk"         "Mk"        
#> [51] "CD16+ Mono" "CD16+ Mono" "B"          "CD16+ Mono" "CD16+ Mono"
#> [56] "CD16+ Mono" "CD16+ Mono" "CD16+ Mono" "CD16+ Mono" "Mk"        
#> [61] "B"          "CD16+ Mono" "B"          "CD16+ Mono" "B"         
#> [66] "CD16+ Mono" "CD16+ Mono" "CD16+ Mono" "CD16+ Mono" "B"         
#> [71] "Mk"         "Mk"         "Mk"         "Mk"         "Mk"        
#> [76] "Mk"         "Mk"         "Mk"         "Mk"         "CD16+ Mono"
```

New reference matrix can be made directly from `SingleCellExperiment`
and `Seurat` objects as well. Other scRNAseq experiment object types are
supported as well.

``` r
# make reference from SingleCellExperiment objects
sce_ref <- object_ref(
  input = sce_small,               # SCE object
  cluster_col = "cell_type1"       # name of column in colData containing cell identities
)
#> The following clusters have less than 10 cells for this analysis: co-expression, ductal, endothelial, epsilon, MHC class II, PSC. Classification is likely inaccurate.

# make reference from seurat objects
s_ref <- seurat_ref(
  seurat_object = s_small3,
  cluster_col = "RNA_snn_res.1"
)

head(s_ref)
#>                 0        1        2
#> MS4A1    0.000000 1.126047 5.151065
#> CD79B    2.469341 2.920407 5.031316
#> CD79A    0.000000 2.535151 5.375681
#> HLA-DRA  3.640368 6.008446 7.055386
#> TCL1A    0.000000 1.495867 4.963367
#> HLA-DQB1 1.603068 3.836290 5.137422
```

`clustify_lists()` handles identity assignment of matrix or
`SingleCellExperiment` and `seurat` objects based on marker gene lists.

``` r
clustify_lists(
  input = pbmc_matrix_small,
  metadata = pbmc_meta,
  cluster_col = "classified",
  marker = pbmc_markers,
  marker_inmatrix = FALSE
)
#>                      0        1        2         3         4        5        6
#> Naive CD4 T  1.5639055 20.19469 31.77095  8.664074 23.844992 19.06931 19.06931
#> Memory CD4 T 1.5639055 20.19469 31.77095 10.568007 23.844992 17.97875 19.06931
#> CD14+ Mono   0.9575077 14.70716 76.21353 17.899569 11.687739 49.86699 16.83210
#> B            0.6564777 12.70976 31.77095 26.422929 13.536295 20.19469 13.53630
#> CD8 T        1.0785353 17.97875 31.82210 12.584823 31.822099 22.71234 40.45383
#> FCGR3A+ Mono 0.6564777 13.63321 72.43684 17.899569  9.726346 56.48245 14.61025
#> NK           0.6564777 14.61025 31.82210  7.757206 31.822099 22.71234 45.05072
#> DC           0.6564777 15.80598 63.34978 19.069308 13.758144 40.56298 17.97875
#> Platelet     0.5428889 13.34769 59.94938 14.215244 15.158755 46.92861 19.49246
#>                      7          8
#> Naive CD4 T   6.165348  0.6055118
#> Memory CD4 T  6.165348  0.9575077
#> CD14+ Mono   25.181595  1.0785353
#> B            17.899569  0.1401901
#> CD8 T         7.882145  0.3309153
#> FCGR3A+ Mono 21.409177  0.3309153
#> NK            5.358651  0.3309153
#> DC           45.101877  0.1401901
#> Platelet     19.492465 59.9493793

clustify_lists(
  input = s_small3,
  marker = pbmc_markers,
  marker_inmatrix = FALSE,
  cluster_col = "RNA_snn_res.1",
  seurat_out = TRUE
)
#> An object of class Seurat 
#> 230 features across 80 samples within 1 assay 
#> Active assay: RNA (230 features, 20 variable features)
#>  2 dimensional reductions calculated: pca, tsne
```

## Additional resources

- [Script](https://github.com/rnabioco/clustifyrdata/blob/master/inst/run_clustifyr.R)
  for benchmarking, compatible with
  [`scRNAseq_Benchmark`](https://github.com/tabdelaal/scRNAseq_Benchmark)

- Additional reference data (including tabula muris, immgen, etc) are
  available in a supplemental package
  [`clustifyrdatahub`](https://github.com/rnabioco/clustifyrdatahub).
  Also see
  [list](https://rnabioco.github.io/clustifyrdata/articles/download_refs.html)
  for individual downloads.

- See the
  [FAQ](https://github.com/rnabioco/clustifyr/wiki/Frequently-asked-questions)
  for more details.
