run_clustifyr <- function(DataPath,
                          LabelsPath,
                          CV_RDataPath,
                          OutputDir,
                          GeneOrderPath = NULL,
                          NumGenes = NULL,
                          threshold = 0.75) {
  "
  run clustifyr
  Wrapper script to run clustifyr on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.

  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection,
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  threshold : whether calls are cut at min correlation coefficients
  "

  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder <- read.csv(GeneOrderPath)
  }

  #############################################################################
  #                              clustifyr                                    #
  #############################################################################
  library(clustifyr)
  library(M3Drop)
  True_Labels_clustifyr <- list()
  Pred_Labels_clustifyr <- list()
  Total_Time_clustifyr <- list()
  Data = t(as.matrix(Data))

  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      start_time <- Sys.time()
      ref = average_clusters(
        Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]],
        Labels[Train_Idx[[i]]]
      )
      meta <- data.frame("cluster" = Labels[Test_Idx[[i]]])
      rownames(meta) <- Test_Idx[[i]]
      res <- clustify(
        input = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]],
        ref = ref,
        metadata = Labels[Test_Idx[[i]]]
      )
      res2 <- cor_to_call(res, threshold = threshold, cluster_col = "cluster")
      meta2 <- call_to_metadata(res2, meta, cluster_col = "cluster")
      end_time <- Sys.time()
    }
    else{
      start_time <- Sys.time()
      Normalized_data <- M3DropCleanData(Data[,Train_Idx[[i]]],
                                         labels = rownames(Data[,Train_Idx[[i]]]),
                                         is.counts = FALSE,
                                         suppress.plot = TRUE
      )
      fits <- M3DropDropoutModels(Normalized_data$data, suppress.plot = TRUE)
      markers_M3Drop <- M3DropFeatureSelection(Normalized_data$data,
                                               mt_method = "fdr",
                                               mt_threshold = 0.01,
                                               suppress.plot = TRUE
      )

      ref = average_clusters(
        Data[,Train_Idx[[i]]],
        Labels[Train_Idx[[i]]]
      )

      meta <- data.frame("cluster" = Labels[Test_Idx[[i]]])
      rownames(meta) <- Test_Idx[[i]]

      res <- clustify(
        input = Data[,Test_Idx[[i]]],
        ref = ref,
        metadata = Labels[Test_Idx[[i]]],
        query_genes = markers_M3Drop$Gene
      )

      res2 <- cor_to_call(
        res,
        threshold = threshold, #0.75,
        cluster_col = "cluster"
      )
      meta2 <- call_to_metadata(
        res2,
        meta,
        cluster_col = "cluster"
      )
      end_time <- Sys.time()
    }
    Total_Time_clustifyr[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_clustifyr[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_clustifyr[i] <- list(meta2$type)
  }
  True_Labels_clustifyr <- as.vector(unlist(True_Labels_clustifyr))
  Pred_Labels_clustifyr <- as.vector(unlist(Pred_Labels_clustifyr))
  Total_Time_clustifyr <- as.vector(unlist(Total_Time_clustifyr))

  setwd(OutputDir)

  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    write.csv(True_Labels_clustifyr,paste('clustifyr_',NumGenes,'_True_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Pred_Labels_clustifyr,paste('clustifyr_',NumGenes,'_Pred_Labels.csv', sep = ''),row.names = FALSE)
    write.csv(Total_Time_clustifyr,paste('clustifyr_',NumGenes,'_Total_Time.csv', sep = ''),row.names = FALSE)
  }
  else{
    write.csv(True_Labels_clustifyr,'clustifyr_True_Labels.csv',row.names = FALSE)
    write.csv(Pred_Labels_clustifyr,'clustifyr_Pred_Labels.csv',row.names = FALSE)
    write.csv(Total_Time_clustifyr,'clustifyr_Total_Time.csv',row.names = FALSE)
  }
}
