
# clustifyr <img src="man/figures/logo.png" align="right">

[![Build
Status](https://travis-ci.org/rnabioco/clustifyr.svg?branch=master)](https://travis-ci.org/rnabioco/clustifyr)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/rnabioco/clustifyr?branch=master&svg=true)](https://ci.appveyor.com/project/rnabioco/clustifyr)
[![Coverage
status](https://codecov.io/gh/rnabioco/clustifyr/branch/master/graph/badge.svg)](https://codecov.io/github/rnabioco/clustifyr?branch=master)

clustifyr classifies cells and clusters in single-cell RNA sequencing
experiments using reference bulk RNA-seq data sets, sorted microarray
expression data, single-cell gene signatures, or lists of marker genes.

Single cell transcriptomes are difficult to annotate without extensive
knowledge of the underlying biology of the system in question. Even with
this knowledge, accurate identification can be challenging due to the
lack of detectable expression of common marker genes. clustifyr solves
this problem by automatically annotating single cells or clusters using
scRNA-seq, bulk RNA-seq data, microarray or marker gene lists.
Additional functions allow for exploratory analysis of similarities
between single cell RNA-seq datasets and reference data.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("rnabioco/clustifyr")
```

## Example usage

In this example we use the following built-in input data:

  - an expression matrix of single cell RNA-seq data
    (`pbmc_matrix_small`)
  - a metadata data.frame (`pbmc_meta`), with cluster information stored
    (`"classified"`)
  - a vector of variable genes (`pbmc_vargenes`)
  - a matrix of mean normalized scRNA-seq UMI counts by cell type
    (`cbmc_ref`):

We then calculate correlation coefficients and plot them on a
pre-calculated tSNE projection (stored in `pbmc_meta`).

``` r
library(clustifyr)

# calculate correlation
res <- clustify(
  input = pbmc_matrix_small,
  metadata = pbmc_meta$classified,
  ref_mat = cbmc_ref,
  query_genes = pbmc_vargenes
)
res2 <- cor_to_call(res)
res2
#> # A tibble: 9 x 3
#> # Groups:   cluster [9]
#>   cluster      type           r
#>   <chr>        <chr>      <dbl>
#> 1 B            B          0.891
#> 2 CD14+ Mono   CD14+ Mono 0.908
#> 3 FCGR3A+ Mono CD16+ Mono 0.906
#> 4 Memory CD4 T CD4 T      0.895
#> 5 Naive CD4 T  CD4 T      0.916
#> 6 DC           DC         0.833
#> 7 Platelet     Mk         0.630
#> 8 CD8 T        NK         0.866
#> 9 NK           NK         0.896

# plotting
plot_best_call(
  cor_mat = res, 
  metadata = pbmc_meta, 
  cluster_col = "classified"
)
```

![](man/figures/readme_example-1.png)<!-- -->

Alternatively, `clustify` can take a clustered `seurat` object (both v2
and v3) and assign identities. Required data can be looked up by default
or specified. New reference matrix can be made directly from `seurat`
object as well. Other scRNAseq experiment object types are supported as
well.

  - a seuratV2 object (`s_small`)
  - or a SeuratV3 object (`s_small3`)
  - a matrix of mean normalized scRNA-seq UMI counts by cell type
    (`cbmc_ref`):

<!-- end list -->

``` r
library(Seurat)

# for seurat2
res <- clustify(
  input = s_small,
  cluster_col = "res.1",
  ref_mat = cbmc_ref,
  seurat_out = TRUE
)

# for Seurat3
res2 <- clustify(
  input = s_small3,
  cluster_col = "RNA_snn_res.1",
  ref_mat = cbmc_ref,
  seurat_out = TRUE
)
res2
#> An object of class Seurat 
#> 230 features across 80 samples within 1 assay 
#> Active assay: RNA (230 features)
#>  2 dimensional reductions calculated: pca, tsne

# make reference from seurat objects
new_ref_matrix <- seurat_ref(
  seurat_object = s_small,
  cluster_col = "res.1"
)

head(new_ref_matrix)
#>                 0        1        2        3
#> MS4A1    4.517255 3.204766 0.000000 0.000000
#> CD79B    4.504191 3.549095 2.580662 0.000000
#> CD79A    4.457349 4.199849 0.000000 0.000000
#> HLA-DRA  6.211779 6.430463 3.659590 4.169965
#> TCL1A    4.394310 2.837922 0.000000 0.000000
#> HLA-DQB1 4.380289 4.325293 0.000000 1.666167
```

Similarly, `clustify_lists` can also handle identity assignment of
matrix or `seurat` object based on marker gene lists.

``` r
res <- clustify_lists(
  input = pbmc_matrix_small,
  metadata = pbmc_meta,
  cluster_col = "classified",
  marker = pbmc_markers,
  marker_inmatrix = FALSE
)

res <- clustify_lists(
  input = s_small,
  marker = pbmc_markers,
  marker_inmatrix = FALSE,
  cluster_col = "res.1",
  seurat_out = TRUE
)
#> cannot find dr info
```

## Additional reference data

More reference data (including tabula muris, immgen, etc) are available
at <https://github.com/rnabioco/clustifyrdata>.

Also see list for individual downloads at
<https://rnabioco.github.io/clustifyr/articles/download_refs.html>
