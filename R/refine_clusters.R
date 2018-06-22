# math model
#   P(x belongs to cluster i)
# = P(x is similar to bulk i)^\lambda  # must similar to a pre-defined cell type
#   * P(x is close to centroid of cluster i)^(1-\lambda)  # must make sure points within a cluster are "closed" to each other

# Algorithm description
# 1. Initiate a cluster using k-means, with k = known cell-types = # bulk data
# 2. For each cluster i and each cell (expr vector x), compute P(x, i). Assign x to the cluster which maximizes P(x,i)
# 3. Create a new cluster and update the centroid of each cluster
# 4. Repeat Step 2-3 until convergence (!!need to prove)

# implementation issues
# 1. log-scale to avoid issue of machine precision
# 2. P(x is similar to bulk i): make use of the similarity score (\in [-1, 1]).
#    Use exponential scoring function P(similarity) = C exp(\epsilon * similarity_score), where C=\epsilon/(exp(\epsilon) - exp(-\epsilon)).
#TODO: Consider changing similarity score to [0, 1] to simplify the expression of C. Probably can ignore because we are comparing relative probabilities
# 3. P(x is closed to centroid of cluster i): make use of Euclidean distance between x and centroid \mu_i
#    Use Euclidean distance because of better change to have convergence (k-means using Euclidean distance is guaranteed to converge to the unique point)
#    Assume something like random walk (i.e. diffusion) -> Gaussian distribution
#    P = 1/sqrt(2*\pi*\sigma^2) exp(-|x-\mu_i|^2/2/\sigma^2)
#TODO: In principle it's a multivariate Gaussian distribution, so \sigma^2 should be replaced by the covariance matrix \Sigma
#    If we assume that all cells/centroids have the same set of genes (which should be!), and every gene has the same amount of error (after normalizing genes!), this expression should be ok.
#    Otherwise, consider model \Sigma explicitly. As a first-order approximation, can assume \Sigma is a diagonal matrix,
#    where each diagonal element represents the stochasticity of each gene (which can be estimate from looking at cluster-wise/data-wise variance)
#    In that case, P \propto exp(-0.5*\sum_j (x_j-\mu_{i,j})^2/\sigma_j^2)/\prod_j \sigma_j, where j iterates over every gene. Note that \sigma_j should at least constant per cluster, in which case \sigma_j -> \sigma_{i,j}

# Simplification using logarithm scale
# ln P(x similar to bulk i)^\lambda = \lambda \epsilon *similarity_score + \lambda ln C   (EQ 1)
# ln P(x closed to centroid i)^(1-\lambda) = -(1-\lambda) 0.5*\sum_j (x_j-\mu_{i,j})^2/\sigma_{i,j}^2 - (1-\lambda)\sum_j ln \sigma_{i,j} - 0.5*(1-\lambda)*n ln(2*\pi)  (EQ 2)
# overall probability
# ln P = (EQ 1) + (EQ 2) = \lambda \epsilon *similarity_score - (1-\lambda) 0.5*\sum_j (x_j-\mu_{i,j})^2/\sigma_{i,j}^2 - (1-\lambda)\sum_j ln \sigma_{i,j} + some constant
# similarity score can be computed by compute_similarity function
# \mu_{i,j} and \sigma_{i,j} can be computed per cluster

# Numerical consideration
# What happened if a cluster runs out of data points? Is it possible that everything will be converge to one cluster, possibly due to any normalization issue? Need a constant penalty factor if that is the case (i.e. P -> P/(P+some constant))

# Testing consideration
# Can we have authentic data? Running permutation definitely works.

