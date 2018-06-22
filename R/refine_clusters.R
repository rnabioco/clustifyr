# math model
#   P(x belongs to cluster i)
# = P(x is similar to bulk i)^\lambda  # must similar to a pre-defined cell type
#   * P(x is close to centroid of cluster i)^(1-\lambda)  # must make sure points within a cluster are "closed" to each other
# Algorithm description
# 1. Initiate a cluster using k-means, with k = known cell-types = # bulk data
# 2. For each cluster i and each cell (expr vector x), compute P(x, i). Assign x to the cluster which maximizes P(x,i)
# 3. Create a new cluster and update the centroid of each cluster
# 4. Repeat Step 2-3 until convergence (!!need to prove)

