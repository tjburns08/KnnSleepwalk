# Code: a helper function for sleepwalk for KNN visualization
# Author: Tyler J Burns
# Date: November 16, 2022

#' @title K-nearest neighbors Sleepwalk
#' @description Takes a data matrix and a 2-D embedding as input. It produces
#' a 'KNN matrix' and places that along with the aforementioned inputs into the
#' sleepwalk function.
#' @param mat A data matrix with data points as rows and features as columns.
#' the matrix must already be filtered by the markers you care about.
#' @param embedding The 2-D embedding of the data matrix, in matrix format with
#' data points as rows and two columns.
#' @param k The number of nearest neighbors to be visualized
KnnSleepwalk <- function(mat, embedding, k = 100, output_file = NULL, point_size = 1.5, dimr_names = c("t-SNE", "UMAP")) {
  message('Building distance matrix')
  dist_mat <- dist(mat) %>% as.matrix()
  message("Finding k-nearest neighbors")
  nn_mat <- lapply(seq(nrow(dist_mat)), function(i) {
    curr <- dist_mat[i,]
    max_dist <- sort(curr, decreasing = FALSE)[k] # This is the K, decreasing set to false
    curr <- ifelse(curr <= max_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)
  sleepwalk::sleepwalk(embeddings = embedding, distances = nn_mat, saveToFile = output_file, pointSize = point_size, titles = dimr_names)
}

#' @title K-farthest neighbors Sleepwalk
#' @description Takes a data matrix and a 2-D embedding as input. It produces
#' a 'KFN matrix' and places that along with the aforementioned inputs into the
#' sleepwalk function.
#' @param mat A data matrix with data points as rows and features as columns.
#' the matrix must already be filtered by the markers you care about.
#' @param embedding The 2-D embedding of the data matrix, in matrix format with
#' data points as rows and two columns.
#' @param k The number of farthest neighbors to be visualized
KfnSleepwalk <- function(mat, embedding, k = 100, output_file = NULL, dimr_names = c("t-SNE", "UMAP")) {
  message('Building distance matrix')
  dist_mat <- dist(mat) %>% as.matrix()
  message("Finding k-farthest neighbors")
  nn_mat <- lapply(seq(nrow(dist_mat)), function(i) {
    curr <- dist_mat[i,]
    min_dist <- sort(curr, decreasing = TRUE)[k] # This is the K, decreasing set to false
    curr <- ifelse(curr >= min_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)
  sleepwalk::sleepwalk(embeddings = embedding, distances = nn_mat, saveToFile = output_file, titles = dimr_names)
}

#' @title K-nearest neighbors Sleepwalk Direct
#' @description Takes a data matrix and a 2-D embedding as input. It produces
#' a 'KNN matrix' and places that along with the aforementioned inputs into the
#' sleepwalk function. Unlike the default KNN sleepwalk function, this one makes
#' two KNN matrices, like UMAP space vs high-D space for comparisons across the
#' same embedding.
#' @param mat1 A data matrix with data points as rows and features as columns.
#' the matrix must already be filtered by the markers you care about.
#' @param mat2 A second data matrix in the format as above
#' @param embedding The 2-D embedding of the data matrix, in matrix format with
#' data points as rows and two columns.
#' @param k The number of nearest neighbors to be visualized
KnnSleepwalkDirect <- function(mat1, mat2, embedding, k = 100, output_file = NULL, point_size = 1.5, dimr_names = c("KNN UMAP space", "KNN high-dim space")) {
  message('Building distance matrix')

  # First distance matrix
  dist_mat1 <- dist(mat1) %>% as.matrix()
  dist_mat2 <- dist(mat2) %>% as.matrix()

  message("Finding k-nearest neighbors for mat1")
  nn_mat1 <- lapply(seq(nrow(dist_mat1)), function(i) {
    curr <- dist_mat1[i,]
    max_dist <- sort(curr, decreasing = FALSE)[k] # This is the K, decreasing set to false
    curr <- ifelse(curr <= max_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)

  # Second distance matrix
  dist_mat2 <- dist(mat2) %>% as.matrix()
  message("Finding k-nearest neighbors for mat2")
  nn_mat2 <- lapply(seq(nrow(dist_mat2)), function(i) {
    curr <- dist_mat2[i,]
    max_dist <- sort(curr, decreasing = FALSE)[k] # This is the K, decreasing set to false
    curr <- ifelse(curr <= max_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)

  sleepwalk::sleepwalk(embeddings = embedding, compare = "distances", distances = list(nn_mat1, nn_mat2), saveToFile = output_file, pointSize = point_size, titles = dimr_names)
}

#' @title K-farthest neighbors Sleepwalk Direct
#' @description Takes a data matrix and a 2-D embedding as input. It produces
#' a 'KFN matrix' and places that along with the aforementioned inputs into the
#' sleepwalk function. Unlike the default KFN sleepwalk function, this one makes
#' two KFN matrices, like UMAP space vs high-D space for comparisons across the
#' same embedding.
#' @param mat1 A data matrix with data points as rows and features as columns.
#' the matrix must already be filtered by the markers you care about.
#' @param mat2 A second data matrix in the format as above
#' @param embedding The 2-D embedding of the data matrix, in matrix format with
#' data points as rows and two columns.
#' @param k The number of nearest neighbors to be visualized
KfnSleepwalkDirect <- function(mat1, mat2, embedding, k = 100, output_file = NULL, point_size = 1.5, dimr_names = c("KFN UMAP space", "KFN high-dim space")) {
  message('Building distance matrix')

  # First distance matrix
  dist_mat1 <- dist(mat1) %>% as.matrix()
  dist_mat2 <- dist(mat2) %>% as.matrix()

  message("Finding k-nearest neighbors for mat1")
  nn_mat1 <- lapply(seq(nrow(dist_mat1)), function(i) {
    curr <- dist_mat1[i,]
    min_dist <- sort(curr, decreasing = TRUE)[k] # This is the K, decreasing set to false
    curr <- ifelse(curr >= min_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)

  # Second distance matrix
  dist_mat2 <- dist(mat2) %>% as.matrix()
  message("Finding k-nearest neighbors for mat2")
  nn_mat2 <- lapply(seq(nrow(dist_mat2)), function(i) {
    curr <- dist_mat2[i,]
    min_dist <- sort(curr, decreasing = TRUE)[k] # This is the K, decreasing set to false
    curr <- ifelse(curr >= min_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)

  sleepwalk::sleepwalk(embeddings = embedding, compare = "distances", distances = list(nn_mat1, nn_mat2), saveToFile = output_file, pointSize = point_size, titles = dimr_names)
}

