# Code: KNN sleepwalk functions
# Author: Tyler J Burns
# Date: May 7, 2024

#' @importFrom magrittr %>%
#' @importFrom stats dist
#' @import sleepwalk
NULL

#' @title K-nearest neighbors Sleepwalk
#' @description Takes a data matrix and a 2-D embedding as input. It produces
#' a 'KNN matrix' and places that along with the aforementioned inputs into the
#' sleepwalk function. Unlike the default KNN sleepwalk function, this one makes
#' two KNN matrices, like UMAP space vs high-D space for comparisons across the
#' same embedding.
#' @param embedding The 2-D embedding of the data matrix, in matrix format with
#' data points as rows and two columns.
#' @param orig_data A data matrix with data points as rows and features as columns.
#' the matrix must already be filtered by the markers you care about.
#' @param k The number of nearest neighbors to be visualized
#' @param output_file The file to save your sleepwalk html page to
#' @param point_size How big you want the point on the plots to be
#' @param plot_names What you want the comparison plots to be named
#' @export
KnnSleepwalk <- function(embedding, orig_data, k = 100, output_file = NULL, point_size = 1.5, plot_names = c("KNN embedding space", "KNN high-dim space")) {
  message('Building distance matrix')

  # Make distance matrices
  dist_mat1 <- dist(embedding) %>% as.matrix()
  dist_mat2 <- dist(orig_data) %>% as.matrix()

  # KNN from the first distance matrix
  message("Finding k-nearest neighbors for the embedding")
  nn_mat1 <- lapply(seq(nrow(dist_mat1)), function(i) {
    curr <- dist_mat1[i,]
    max_dist <- sort(curr, decreasing = FALSE)[k] # This is the K, decreasing set to false
    curr <- ifelse(curr <= max_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)

  # KNN from the second distance matrix
  message("Finding k-nearest neighbors for original data")
  nn_mat2 <- lapply(seq(nrow(dist_mat2)), function(i) {
    curr <- dist_mat2[i,]
    max_dist <- sort(curr, decreasing = FALSE)[k] # This is the K, decreasing set to false
    curr <- ifelse(curr <= max_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)

  sleepwalk::sleepwalk(embeddings = embedding,
                       compare = "distances",
                       distances = list(nn_mat1, nn_mat2),
                       saveToFile = output_file,
                       pointSize = point_size,
                       titles = plot_names)
}

#' @title K-farthest neighbors Sleepwalk Direct
#' @description Takes a data matrix and a 2-D embedding as input. It produces
#' a 'KFN matrix' and places that along with the aforementioned inputs into the
#' sleepwalk function. Unlike the default KFN sleepwalk function, this one makes
#' two KFN matrices, like UMAP space vs high-D space for comparisons across the
#' same embedding.
#' @param embedding The 2-D embedding of the data matrix, in matrix format with
#' data points as rows and two columns.
#' @param orig_data A data matrix with data points as rows and features as columns.
#' the matrix must already be filtered by the markers you care about.
#' @param k The number of nearest neighbors to be visualized
#' @param output_file The file to save your sleepwalk html page to
#' @param point_size How big you want the point on the plots to be
#' @param plot_names What you want the comparison plots to be named
#' @export
KfnSleepwalk <- function(embedding, orig_data, embedding, k = 100, output_file = NULL, point_size = 1.5, plot_names = c("KFN embedding space", "KFN high-dim space")) {
  message('Building distance matrix')

  # First distance matrix
  dist_mat1 <- dist(embedding) %>% as.matrix()
  dist_mat2 <- dist(orig_data) %>% as.matrix()

  message("Finding k-nearest neighbors for the embedding")
  nn_mat1 <- lapply(seq(nrow(dist_mat1)), function(i) {
    curr <- dist_mat1[i,]
    min_dist <- sort(curr, decreasing = TRUE)[k] # This is the K, decreasing set to true
    curr <- ifelse(curr >= min_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)

  # Second distance matrix
  message("Finding k-nearest neighbors for the original data")
  nn_mat2 <- lapply(seq(nrow(dist_mat2)), function(i) {
    curr <- dist_mat2[i,]
    min_dist <- sort(curr, decreasing = TRUE)[k] # This is the K, decreasing set to false
    curr <- ifelse(curr >= min_dist, curr, 1000) # A large number
    return(curr)
  }) %>% do.call(rbind, .)

  sleepwalk::sleepwalk(embeddings = embedding,
                       compare = "distances",
                       distances = list(nn_mat1, nn_mat2),
                       saveToFile = output_file,
                       pointSize = point_size,
                       titles = plot_names)
}

