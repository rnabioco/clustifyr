[![Build Status](https://travis-ci.org/NCBI-Hackathons/clustifyR.svg?branch=master)](https://travis-ci.org/NCBI-Hackathons/clustifyR)

<p align="center">
  <img src="/inst/logo/logo_transparent.png">
</p>

### ClustifyR classifies cell and clusters in single-cell RNA sequencing experiments using reference bulk RNA-seq data sets, gene signatures or marker genes. 

Single cell transcriptomes are difficult to annotate without extensive knowledge of the underlying biology of the system in question. Even with this knowledge, accurate identification can be challenging due to the lack of detectable expression of common marker genes defined by bulk RNA-seq, flow cytometry, etc. `ClustifyR` solves this problem by providing functions to automatically annotate single cells or clusters using bulk RNA-seq data or marker gene lists (ranked or unranked). Additional functions allow for exploratory analysis of similarities between single cell RNA-seq datasets and reference data. Put another way:

**C**lustifyR **L**everages **U**ser **S**upplied **T**ranscripts to **I**dentify **F**eatures in **Y**our sc**R**NA-seq

## Installation
Installation in R from the github repo is a two step process:

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

## Usage:

### See detailed vignettes and documentation [here](https://ncbi-hackathons.github.io/clustifyR/).

