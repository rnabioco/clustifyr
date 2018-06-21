# test correlation function
# remember to execute functions at the end of the document
data("pbmc4k_avg");
data("pbmc_bulk_matrix");
gene_constraints <- list(rownames(pbmc4k_avg), rownames(pbmc_bulk_matrix));
pbmc4k_avg <- select_gene_subset(pbmc4k_avg, gene_constraints);
pbmc_bulk_matrix <- select_gene_subset(pbmc_bulk_matrix, gene_constraints);
compute_similarity(pbmc4k_avg[,"1"], pbmc_bulk_matrix[,"primary human monocytes rep 1"], kl_divergence)

# pick and sort genes
# expr_matrix: expression matrix with row names as the gene names (short name for now)
# gene_constraints: a list of vectors, where each vector is a candidate list of selected genes. can be HVG, row names of scRNAseq data and/or bulk data.
select_gene_subset <- function(expr_matrix, gene_constraints) {
  gene_subset <- gene_constraints[[1]];
  for (i in 2:length(gene_constraints)) {
    gene_subset <- intersect(gene_subset, gene_constraints[[i]]);
  }
  return(expr_matrix[sort(gene_subset),]);
}

# Compute the similarity score between two vectors using a customized scoring function
# vec1, vec2: input expression vectors, can come from single cell or bulk RNA seq data.
# the length of vec1 and vec2 must match, and corresponding elements must refer to the same gene.
# vec1 and vec2 should be AFTER pre-processing, feature selection and dimension reduction
# compute_method and ...: function to compute the similarity score, and relevant parameters for compute_method
compute_similarity <- function(vec1, vec2, compute_method, ...) {
  # examine whether two vectors are of the same size
  if (!is.numeric(vec1) || !is.numeric(vec2) || length(vec1) != length(vec2)) {
    stop("compute_similarity: two input vectors are not numeric or of different sizes.");
  }

  # return the similarity score, must be
  return(compute_method(vec1, vec2, ...));
}

## candidate scoring methods for similarity.
## score must be between -1 and 1.
## 1. completely positively correlated
## -1. completely negatively correlated
corr_coef <- function(vec1, vec2, method="pearson") {
  return(switch(method,
                pearson=cor(vec1, vec2, method="pearson"),
                spearman=cor(vec1, vec2, method="spearman"),
                cosine=sum(vec1*vec2)/sqrt(sum(vec1^2)*sum(vec2^2))
  ));
}
# use package entropy to compute kl divergence
kl_divergence <- function(vec1, vec2, if_logcounts=FALSE, total_reads=1000, max_KL=1) {
  library(entropy); # load entropy library
  if (if_logcounts) {
    vec1 <- exp(vec1); vec2 <- exp(vec2);
  }
  count1 <- round(vec1*total_reads/sum(vec1)); count2 <- round(vec2*total_reads/sum(vec2));
  est_KL <- KL.shrink(count1, count2, unit="log2")
  return((max_KL-est_KL)/max_KL*2-1);
}

