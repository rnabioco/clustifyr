library(usethis)

object_loc_lookup$SingleCellExperiment <- c(
  expr = "input@assays$data$logcounts",
  meta = "as.data.frame(input@colData)",
  var = NULL,
  col = "cell_type1"
)

object_loc_lookup$URD <- c(
  expr = "input@logupx.data",
  meta = "input@meta",
  var = "input@var.genes",
  col = "cluster"
)

object_loc_lookup$FunctionalSingleCellExperiment <- c(
  expr = "input@ExperimentList$rnaseq@assays$data$logcounts",
  meta = "input@ExperimentList$rnaseq@colData",
  var = NULL,
  col = "leiden_cluster"
)

object_loc_lookup$Seurat <- c(
  expr = "input@assays$RNA@data",
  meta = "input@meta.data",
  var = "input@assays$RNA@var.features",
  col = "RNA_snn_res.1"
)

object_loc_lookup$CellDataSet <- c(
  expr = "do.call(function(x) {row.names(x) <- input@featureData@data$gene_short_name; return(x)}, list(input@assayData$exprs))",
  meta = "as.data.frame(input@phenoData@data)",
  var = "as.character(input@featureData@data$gene_short_name[input@featureData@data$use_for_ordering == T])",
  col = "Main_Cluster"
)

usethis::use_data(object_loc_lookup, compress = "xz", overwrite = TRUE)
