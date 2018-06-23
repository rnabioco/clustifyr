[![Build Status](https://travis-ci.org/NCBI-Hackathons/clustifyR.svg?branch=master)](https://travis-ci.org/NCBI-Hackathons/clustifyR)

<p align="center">
  <img src="/inst/logo/logo_transparent.png">
</p>

### ClustifyR is an [R](https://www.r-project.org/) package that classifies cells and clusters in single-cell RNA sequencing experiments using reference bulk RNA-seq data sets, gene signatures or marker genes. 

Single cell transcriptomes are difficult to annotate without extensive knowledge of the underlying biology of the system in question. Even with this knowledge, accurate identification can be challenging due to the lack of detectable expression of common marker genes defined by bulk RNA-seq, flow cytometry, etc. `ClustifyR` solves this problem by providing functions to automatically annotate single cells or clusters using bulk RNA-seq data or marker gene lists (ranked or unranked). Additional functions allow for exploratory analysis of similarities between single cell RNA-seq datasets and reference data. Put another way:

**C**lustifyR **L**everages **U**ser **S**upplied **T**ranscripts to **I**dentify **F**eatures in **Y**our sc**R**NA-seq

## Installation
Installation from github in R is a two step process:

### Step 1:
```
# Install devtools
> install.packages("devtools")
```

### Step 2:
```
# Install classifyR from github
> devtools::install_github("NCBI-Hackathons/clustifyR")
```

## Usage

### Super-quickstart with sample data (included in package):

Generate a correlation matrix from a matrix of single cell RNA-seq data (`pbmc4k_matrix`), a metadata table describing the single cell data (`pbmc4k_meta`), a list of variable genes in the single cell data (`pbmc4k_vargenes`), and a matrix of bulk RNA-seq read counts (`pbmc_bulk_matrix`):

```
# run correlation (pearson by default)
> res <- run_cor(expr_mat = pbmc4k_matrix,
                 metadata = pbmc4k_meta,
                 bulk_mat = pbmc_bulk_matrix,
                 query_gene_list = pbmc4k_vargenes,
                 compute_method = corr_coef)
```

Plot the correlation coefficients on a pre-calculated tSNE projection (stored in `pbmc4k_meta`):

```
# plot correlation coefficients on tSNE for each identity class
> plot_cor(res,
           pbmc4k_meta,
           colnames(res)[c(1, 5)],
           cluster_col = "classified")
```

### For more infomation, see the [detailed vignettes and documentation](https://ncbi-hackathons.github.io/clustifyR/).
