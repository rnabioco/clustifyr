appendGenes <- function(geneVector, GSEMatrix)
{
  rownamesGSEMatrix <- rownames(GSEMatrix) #Get rownames from GSEMatrix (new GSE file)

  rowCountHumanGenes <- nrow(geneVector) #Calculate number of rows from list of full human genes
  rowCountNewGSEFile <- nrow(GSEMatrix) #Calculate number of rows of GSE matrix

  missing_rows <- setdiff(geneVector, rownamesGSEMatrix) #Use setdiff function to figure out rows which are different/missing from GSE matrix

  zeroExpressionMatrix <- matrix(0, nrow = length(missing_rows), ncol = ncol(GSEMatrix)) #Create a placeholder matrix with zeroes and missing_rows length as row length

  rownames(zeroExpressionMatrix) <- missing_rows #Assign row names
  colnames(zeroExpressionMatrix) <- colnames(GSEMatrix) #Assign column names

  fullMatrix <- rbind(GSEMatrix, zeroExpressionMatrix) #Bind GSEMatrix and zeroExpressionMatrix together

  #Reorder matrix
  fullMatrix <- fullMatrix[geneVector, ] #Reorder fullMatrix to preserve gene order
  return(fullMatrix) #Return fullMatrix
}

checkRawCounts <- function(GSEMatrix, max_log_value = 50)
{
  if(!is.matrix(GSEMatrix))
  {
    GSEMatrix <- as.matrix(GSEMatrix)
  }
  if (is.integer(GSEMatrix))
  {
    return("raw counts")
  }
  else if (is.double(GSEMatrix))
  {
    if (all(GSEMatrix == floor(GSEMatrix)))
    {
      return("raw counts")
    }
    if(max(GSEMatrix) > max_log_value)
    {
      return("normalized")
    }
    else if (min(GSEMatrix) < 0)
    {
      stop("negative values detected, likely scaled data")
    }
    else
    {
      return("log-normalized")
    }
  }
  else
  {
    stop("unknown matrix format: ", typeof(GSEMatrix))
  }
}


#' Function to combine records into single atlas
#'
#' @param matrix_fns character vector of paths to study matrices stored as .rds files.
#' If a named character vector, then the name will be added as a suffix to the cell type
#' name in the final matrix. If it is not named, then the filename will be used (without .rds)
#' @param genes_fn text file with a single column containing genes and the ordering desired
#' in the output matrix
#' @param matrix_objs Checks to see whether .rds files will be read or R objects in a local environment 
#' @param output_fn output filename for .rds file. If NULL the matrix will be returned instead of
#' saving
build_atlas <- function(matrix_fns = NULL,
                        genes_fn,
                        matrix_objs = NULL,
                        output_fn = NULL)
{
  genesVector <- read_lines(genes_fn)
  if(is.null(matrix_obs) && !is.null(matrix_fns))
    {
      ref_mats <- lapply(matrix_fns, readRDS)
      if(is.null(names(matrix_fns)))
        {
          names(ref_mats) <- basename(ref_matrices_fns) %>% str_remove(".rds$")
        } 
      else 
      {
        names(ref_mats) <- names(matrix_fns)
      }
    } 
  else if(!is.null(matrix_obs)) 
    {
    ref_mats <- matrix_obs
      if(is.null(names(matrix_obs)))
        {
          names(ref_mats) <- basename(ref_matrices_fns) %>% str_remove(".rds$")
        } 
      else 
        {
          names(ref_mats) <- names(matrix_obs)
        }
    }
    new_mats <- list()
    for(i in seq_along(ref_mats))
      {
        # standardize genes in matrix
        mat <- appendGenes(geneVector = genesVector,
                           GSEMatrix = as.matrix(ref_mats[[i]]))
        # get study name
        mat_name <- names(ref_mats)[i]
        
        # append study name to cell type names
        new_cols <- paste0(colnames(mat),
                           " (",
                           mat_name,
                           ")")
        colnames(mat) <- new_cols
        
        # assign to list
        new_mats[[i]] <- mat
      }
      
    # cbind a list of matrices
    atlas <- do.call(cbind, new_mats)
      
      if(!is.null(output_fn))
        {
          saveRDS(atlas, output_fn)
        } 
      else 
        {
          return(atlas)
        }
}


