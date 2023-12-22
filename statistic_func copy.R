library(dplyr)
#library(mvtnorm)
#library(nnet)
#library('MASS')
#library('LambertW')

# GENERAL FUNCTIONS
f_lp_distance <- function(point1, point2, p) {
  return(sum(abs(point1 - point2)^p)^(1/p))
}
f_enumerate <- function(X){
  cbind(1:length(X), X)
}

f_logistic <- function(x){
  1/(1+exp(-x))
}

f_logit<- function(p){
  log(p/(1-p))
}

my_maha2 = function(x, mu, sig){
  t(x-mu)%*%solve(sig)%*%(x-mu)
}

f_frob_pnorm <- function(A){
  sum = 0 
  for (i in 1:(nrow(A))){
    for (j in 1:(ncol(A))){
      sum = sum+ (A[i,][j])^(2)
    }
  }
  sqrt(sum)
}

f_lpnorm <- function(vec, p){
  (sum(abs(vec)^p))^(1/p)
}

# HYPOTHESIS TESTING
#install.packages("mvtnorm")
#library(mvtnorm) # Multivariate Normal Function





# STATISTICAL TEST FUNCTIONS
roy_stat <- function(sigma_SSE, sigma_SSM){
  max(eigen(solve(sigma_SSE) %*% sigma_SSM)$values)
}

pillai_stat <- function(sigma_SSE, sigma_SSM){
  sum(diag(solve(sigma_SSE + sigma_SSM) %*% sigma_SSM))
}

wilk_stat <- function(n, g, d, sigma_SSE, sigma_SST){
  (n-g -0.5*(d+2))*log(det(sigma_SST)/det(sigma_SSE))
}

m_wilk_stat <- function(n, g, d, sigma_SSE, sigma_SST){
  (n-g -d +1)*log(det(sigma_SST)/det(sigma_SSE))/(d*(g-1))
}

hotelling_stat <- function(sigma_SSE, sigma_SSM){
  sum(diag(solve(sigma_SSE) %*% sigma_SSM))
}

# CLUSTER INITIALIZATION STRATEGY

f_initialBBox <- function(K, X){
  min_ls <- apply(X,2,min)
  max_ls <- apply(X,2,max)
  centroids <- matrix(nrow = K, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    centroids[, i] <- runif(K, min = min_ls[i], max = max_ls[i])
  }
  
  return(centroids)
}

f_initialForgy <- function(K, X){
  X[sample(nrow(X), K, replace = FALSE), ]
}

f_initialRandomPartition <- function(K, X){
  # Randomly assign ids 
  n = nrow(X)
  indices = sample(1:K, n, replace = TRUE)
  newX <- rbind(X, indices)
  
  # Define new centroids
  centroids = c()
  for (i in 1:K){
    filtered_obs <- newX[newX$indices == i, ]
    centroids <- cbind(centroids, colMeans(filtered_obs))
  }
  return(centroids)
}

# CLUSTERING METHODS

# Hierarchical 

f_hierarchical_cluster <- function(X, threhold, p){
  n = nrow(X)
  clusters <- list()
  cluster_labels <- rep(NA, n)
  for (i in 1:n){
    if (is.na(cluster_labels[i])) { # Not Clustered Yet
      # Current point becomes the centroid of a new cluster
      centroid <- X[i,]
      # Initialize an empty cluster
      new_cluster <- c(i)
      cur_id <-  length(clusters) + 1
      cluster_labels[i] <- cur_id
      for (j in i+1:n){
        if (is.na(cluster_labels[j])) { # Not Clustered Yet
          distance <- f_lp_distance(centroid, X[j,], p)
          if(distance <= threshold){
            # Assign obs cluster id
            new_cluster <- c(new_cluster, j) 
            cluster_labels[j] <- cur_id
          }
        }
      }
      #Add set of indices corresponding to cluster
      clusters[[length(clusters) + 1]] <- new_cluster
    }
  }
  # Clean up the clusters list to remove empty entries
  clusters <- clusters[!sapply(clusters, is.null)]
  
  # Output the number of clusters and their sizes
  num_clusters <- length(clusters)
  #cat("Number of clusters:", num_clusters, "\n")
  cluster_sizes <- sapply(clusters, length)
  #cat("Sizes of each cluster:", cluster_sizes, "\n")
  
  return(clusters, num_clusters, cluster_sizes)
}

# Lloyd 
# Assume K is number of clusters, X is dataset, C is initial centroids, and p is the norm
f_lloyd_one_iteration <- function(K, X, C, p) {
  centroids <- C
  clusters <- rep(0, nrow(X))
  
  # Cluster assignment step
  new_clusters <- apply(X, 1, function(x) {
    which.min(sapply(1:K, function(j) sum(abs(x - centroids[j,])^p)^(1/p)))
  })
  
  clusters <- new_clusters
  for (k in 1:K) {
    if (any(clusters == k)) {
      centroids[k, ] <- colMeans(X[clusters == k, ])
    }
  }
  
  list(centroids = centroids, clusters = clusters)
}

# Assume K is number of clusters, X is dataset, C is initial centroids, and p is the norm
f_lloyd <- function(K, X, C, p){
  centroids <- C
  clusters <- rep(0, nrow(X))
  while(TRUE){
    new_clusters <- apply(X, 1, function(x) {
      which.min(sapply(1:K, function(j) sum(abs(x - centroids[j,])^p)^(1/p)))
    })
    if (identical(new_clusters, clusters)) {
      break
    }
    clusters <- new_clusters
    # Update centroids
    centroids <- sapply(unique(clusters), function(c) {
      colMeans(X[clusters == c, ])
    }, simplify = 'data.frame') 
  }
  
  list(centroids = centroids, clusters = clusters)
}

# Assume K is number of clusters, X is dataset, C is initial centroids, and p is the norm
f_macqueen <- function(K, X, C, p){
  # high level: introduce a point, centroid, another point, centroid
  # We want to assume that the dataset is ordered (for random, randomize). 
  n <- nrow(X)
  clusters <- rep(NA, n)
  centroids <- C
  switch_ls <- rep(0, K)
  for (i in 1:n){
    obs <- X[i,]
    c_id <- which.min(sapply(1:K, function(j) f_lp_distance(obs, centroids[j], p)))
    clusters[i] = c_id
    if(switch_ls[c_id] == 0){
      switch_ls[c_id] = 1
      centroids[c_id] = obs
    }else{ #non-zero
      centroids[c_id] = (as.double(centroids[c_id])*switch_ls[c_id] + obs)/(switch_ls[c_id] + 1)
      switch_ls[c_id] = switch_ls[c_id] + 1
    }
  }
  return(list(centroids = centroids, clusters = clusters))
}

# Fuzzy Clustering


# CLASSIFICATION
# Logistic Regression - Classification
# Assume X is the dataset of the v1s and features, and Y are the labels
f_logistic_reg_model <- function(X, Y) {
  # Create a formula dynamically
  feature_data <- X[,2:ncol(X)]
  predictor_vars <- colnames(X)[2:ncol(X)] # Omit the v1s
  formula_str <- paste("Y ~", paste(predictor_vars, collapse = " + "))
  logistic_formula <- as.formula(formula_str)
  # Fit the logistic model
  logistic_model <- glm(logistic_formula, data = as.data.frame(cbind(Y, feature_data)), family = 'binomial')
  # Predict probabilities
  beta_hat <- logistic_model$coefficients
  p_hat <- 1/(1 + exp(-X%*%beta_hat))
  # Calculate statistics
  ssr <- sum((Y - p_hat)^2)
  rmse <- sqrt(ssr)
  mean_Y <- mean(Y)
  # Return a list of statistics
  return(list(ssr = ssr, rmse = rmse, mean_Y = mean_Y))
}

# K nearest neighbors (KNN) - classification
# Assume that K is the number of neighbors, X is data, Y are the labels.
f_knn_classifier <- function(K, x_train, Y_train, x_test, p){
  predictions <- c()
  agg_data <- cbind(x_train, Y_train)
  
  if (!is.null(nrow(x_test))){
    for (i in 1:nrow(x_test)){
      distances <- apply(x_train, 1, function(x) f_lp_distance(x, x_test[i,], 2))
      neighbors <- order(distances)[1:K]
      neighbor_labels <- agg_data[neighbors, ncol(agg_data)]
      predictions <- c(predictions,  names(sort(table(neighbor_labels), decreasing = TRUE))[1])#max(agg_data[neighbors, ncol(agg_data)]) )
    }
  }else{ #NULL

    distances <- apply(x_train, 1, function(x) f_lp_distance(x, x_test, 2))

    enum_dist <- f_enumerate(distances)

    neighbors <- order(distances)[1:K]

    neighbor_labels <- agg_data[neighbors, ncol(agg_data)]
    print(neighbor_labels)
    predictions <- c(predictions,  names(sort(table(neighbor_labels), decreasing = TRUE))[1])#max(agg_data[neighbors, ncol(agg_data)]) )
    print(names(sort(table(neighbor_labels), decreasing = TRUE))[1])
  }
  list(predictions, neighbor_labels)
}


f_SSC <- function(centroid, X){
  
  centroid_mat <- matrix(rep(centroid, nrow(X)), 
                         nrow = nrow(X), 
                         ncol = ncol(X), 
                         byrow = TRUE)
  sum((X - centroid_mat)^2)
}

f_gini_impurity <- function(bag){
  prob_ls <- table(bag)/length(bag)
  return(1- sum(prob_ls^2))
}

f_entropy <- function(bag, b = exp(1)){
  prob_ls <- table(bag)/length(bag)
  return(-1*sum(prob_ls*log(prob_ls, base = b)))
}

# Assume that X is the dataset with relevant features, Y is the ground truth, 
# newobs is the new dataset on which we run the decision tree
f_exhaustive_decision_tree <- function(X, Y, newobs){ 
  # Let's just do it layerwise by features in the following order
  # Features
  #out <- list()
  #n <- nrow(X)
  #for (i in 1:n){
  #  layer_feature <- X[, i]
  # unique
  # which
  #}
  return(NA)
}

f_quick_exhaustive_decision_tree <- function(X, Y, newObs){
  matches <- apply(X, 1, function(row) all(row == newObs))
  # Subset Y based on matches
  matched_Y <- Y[matches]
  return(table(matched_Y)/length(matched_Y))
}

f_random_tree <-  function(X){
  # Bootstrap
  
  return(NA)
}
