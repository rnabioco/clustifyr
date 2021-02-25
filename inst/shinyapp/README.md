# Clustifyr Web App <img src="logo.png" align="right">
The purpose of the Clustifyr Web App is to enable quick benchmarking and clustering of cells in single-cell RNA sequencing experiments through an interactive web interface and application. Users can directly upload their matrix and metadata files generated from single-cell RNA sequencing data and produce useful clustering data tables and visualizations. Using a built in library of references compiled from reference bulk RNA-seq experiments, microarray expression data, single-cell gene signatures, and marker gene lists, users can also conduct analysis of similarities between their novel data and reference data. An additional purpose of the Clustifyr Web App is to make sure that investigators have accurate metadata as part of their NCBI Gene Expression Omnibus (GEO) records. If no metadata is present from a GEO record, the user can directly auto-email the investigator to supply the sufficient metadata.

## Clustifyr Background
clustifyr classifies cells and clusters in single-cell RNA sequencing
experiments using reference bulk RNA-seq data sets, sorted microarray
expression data, single-cell gene signatures, or lists of marker genes.

Single cell transcriptomes are difficult to annotate without knowledge
of the underlying biology. Even with this knowledge, accurate
identification can be challenging due to the lack of detectable
expression of common marker genes. clustifyr solves this problem by
automatically annotating single cells or clusters of cells using
single-cell RNA-seq, bulk RNA-seq data, microarray, or marker gene
lists. Additional functions enable exploratory analysis of similarities
between single cell RNA-seq datasets and reference data.

## Atlas and Clustifyr Data Hub
The reference bulk RNA-seq data sets, sorted microarray
expression data, single-cell gene signatures, and lists of marker genes
are located in the accompanying [Clustifyr Data Hub](https://github.com/rnabioco/clustifyrdatahub).
Descriptions of each data set are present in the table below. 

## Available references include

| Title                    | Species      | Description                                      | RDataPath                                     | BiocVersion | Genome | SourceType | SourceUrl                                                                                            |
| :----------------------- | :----------- | :----------------------------------------------- | :-------------------------------------------- | ----------: | :----- | :--------- | :--------------------------------------------------------------------------------------------------- |
| ref\_MCA                 | Mus musculus | Mouse Cell Atlas                                 | clustifyrdatahub/ref\_MCA.rda                 |        3.12 | mm10   | Zip        | <https://ndownloader.figshare.com/files/10756795>                                                    |
| ref\_tabula\_muris\_drop | Mus musculus | Tabula Muris (10X)                               | clustifyrdatahub/ref\_tabula\_muris\_drop.rda |        3.12 | mm10   | Zip        | <https://ndownloader.figshare.com/articles/5821263>                                                  |
| ref\_tabula\_muris\_facs | Mus musculus | Tabula Muris (SmartSeq2)                         | clustifyrdatahub/ref\_tabula\_muris\_facs.rda |        3.12 | mm10   | Zip        | <https://ndownloader.figshare.com/articles/5821263>                                                  |
| ref\_mouse.rnaseq        | Mus musculus | Mouse RNA-seq from 28 cell types                 | clustifyrdatahub/ref\_mouse.rnaseq.rda        |        3.12 | mm10   | RDA        | <https://github.com/dviraran/SingleR/tree/master/data>                                               |
| ref\_moca\_main          | Mus musculus | Mouse Organogenesis Cell Atlas (main cell types) | clustifyrdatahub/ref\_moca\_main.rda          |        3.12 | mm10   | RDA        | <https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads>                    |
| ref\_immgen              | Mus musculus | Mouse sorted immune cells                        | clustifyrdatahub/ref\_immgen.rda              |        3.12 | mm10   | RDA        | <https://github.com/dviraran/SingleR/tree/master/data>                                               |
| ref\_hema\_microarray    | Homo sapiens | Human hematopoietic cell microarray              | clustifyrdatahub/ref\_hema\_microarray.rda    |        3.12 | hg38   | TXT        | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24759/matrix/GSE24759_series_matrix.txt.gz>     |
| ref\_cortex\_dev         | Homo sapiens | Human cortex development scRNA-seq               | clustifyrdatahub/ref\_cortex\_dev.rda         |        3.12 | hg38   | TSV        | <https://cells.ucsc.edu/cortex-dev/exprMatrix.tsv.gz>                                                |
| ref\_pan\_indrop         | Homo sapiens | Human pancreatic cell scRNA-seq (inDrop)         | clustifyrdatahub/ref\_pan\_indrop.rda         |        3.12 | hg38   | RDA        | <https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/baron-human.rds>                   |
| ref\_pan\_smartseq2      | Homo sapiens | Human pancreatic cell scRNA-seq (SmartSeq2)      | clustifyrdatahub/ref\_pan\_smartseq2.rda      |        3.12 | hg38   | RDA        | <https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/segerstolpe.rds>                   |
| ref\_mouse\_atlas        | Mus musculus | Mouse Atlas scRNA-seq from 321 cell types        | clustifyrdatahub/ref\_mouse\_atlas.rda        |        3.12 | mm10   | RDA        | <https://github.com/rnabioco/scRNA-seq-Cell-Ref-Matrix/blob/master/atlas/musMusculus/MouseAtlas.rda> |
