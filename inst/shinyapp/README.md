# Clustifyr Shiny App <img src="logo.png" align="right">

The purpose of this app is to enable quick classification of single-cell RNA sequencing data through an interactive web interface. Users can directly upload their matrix and metadata files or Seurat/SCE objects generated from single-cell RNA sequencing analyses and produce useful cell identity inference and visualization,using a built-in library of references (clustifyrdatahub) compiled from reference bulk RNA-seq experiments, microarray expression data, and single-cell gene signatures. An additional purpose of the app is to enable quick browsing, preview, and reference building directly from NCBI Gene Expression Omnibus (GEO) records. Data reuse for this application and many other reanalysis/extension purposes require accurate metadata, which is frustratingly rare. We call on data repositories, journals, and investigators to work together towards ensuring proper cell-level annotation deposition. Please see [someta](https://github.com/rnabioco/someta) for further discussions.

## Workflow

1. Upload expression, either raw counts or normalized, matrix. Or Seurat/SCE object. You can also retrieve GEO data through accession #.
2. Upload cell-level metadata in text formats. Or Seurat/SCE object. You can also retrieve GEO data through accession #.
3. Choose (in dropdown menu or just click on the preview) the column in metadata that represent clustering information.
4. Choose or upload reference dataset, a matrix containing average expression of each cell type. Mouse MCA is the default.
5. Go to clustify step and look at results: correlation matrix, called cell-types, and heatmap. Results can be downloaded as xlsx.

## Clustifyr Background

Single cell transcriptomes are difficult to annotate without knowledge
of the underlying biology. Even with this knowledge, accurate
identification can be challenging due to the lack of detectable
expression of common marker genes. [clustifyr](https://github.com/rnabioco/clustifyr) aims to alleviate this problem by
automatically annotating single cells or clusters of cells using
single-cell RNA-seq, bulk RNA-seq data, microarray, or marker gene
lists. Additional functions enable exploratory analysis of similarities
between single cell RNA-seq datasets and reference data.

## Clustifyr Data Hub

Reference cell type gene signatures are located in the accompanying [Clustifyr Data Hub](https://github.com/rnabioco/clustifyrdatahub).
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
