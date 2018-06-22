# clustifyR

[![Travis build status](https://travis-ci.org/NCBI-Hackathons/clustifyR.svg?branch=master)](https://travis-ci.org/NCBI-Hackathons/clustifyR)

### Classify cell types in single-cell RNA sequencing expeirments using
reference bulk RNA-seq data sets. See detailed vignettes and documentation
[here](https://ncbi-hackathons.github.io/clustifyR/).


## Workflow

1. Load a single-cell RNA-seq data set (`SparseMatrix` or `data.frame`). Normalize (?) with `scran`.

1. Load the reference data type (i.e., bulk RNA-seq data, a gene set (ordered or unordered))

1. Compute correlations between single-cell and reference data types. Classify on per-cell or per-cluster basis.

1. Output per-cell or per-cluster, significance of the assignment.

1. Visualize per-cell, per-cluster correlations in e.g. a t-SNE.

## Previous work

* `scmap`: [[code](https://github.com/hemberg-lab/scmap)] [[paper](https://www.nature.com/articles/nmeth.4644)]
  
* `scfind` [[code](https://github.com/hemberg-lab/scfind)]

* `SingleR` [[code](https://github.com/dviraran/SingleR)] [[paper](https://doi.org/10.1101/284604)]

* `CellAtlasSearch` [[paper](https://doi.org/10.1093/nar/gky421)] (no code available)
