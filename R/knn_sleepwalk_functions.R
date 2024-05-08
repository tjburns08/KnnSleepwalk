# Code: KNN sleepwalk functions
# Author: Tyler J Burns
# Date: May 7, 2024

#' @importFrom magrittr %>%
#' @importFrom stats dist
#' @import sleepwalk
NULL

#' @title Make nearest neighbor matrix
#' @description A helper function for KnnSleepwalk which takes a matrix as
#' input, computes the distance matrix, and then changes the values such that
#' the output of the KnnSleepwalk will show the K-nearest or fartherst neighbors
#' as black, with the rest of the cells being a lighter color. This is what we
#' call the nearest neighbor matrix.
#' @param mat A data matrix, presumably from a single-cell dataset
#' @param kfn Whether you want to look at the K-nearest or K-farthest neighbors
#' @param k The number of nearest neighbors
#' @return The afrorementioned nearest neighbor matrix
#' @keywords internal
MakeNnMatrix <- function(mat, kfn = FALSE, k = 100) {
  dist_mat <- dist(mat) %>% as.matrix()
  nn_mat <- lapply(seq(nrow(dist_mat)), function(i) {
    curr <- dist_mat[i,]
    if(kfn) {
      # KFN case
      working_dist <- sort(curr, decreasing = TRUE)[k]
      curr <- ifelse(curr <= working_dist, 1000, 0)
    } else {
      # KNN case
      working_dist <- sort(curr, decreasing = FALSE)[k]
      curr <- ifelse(curr >= working_dist, 1000, 0)
    }
    return(curr)
  }) %>% do.call(rbind, .)
  return(nn_mat)
}

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
#' @param kfn Whether you want to look at the k-farthest neighbors. If false (default) then you look at the K-nearest nighbors.
#' @export
KnnSleepwalk <- function(embedding, orig_data, k = 100, output_file = NULL, point_size = 1.5, plot_names = c("KNN embedding space", "KNN high-dim space"), kfn = FALSE) {
  message('Building distance matrix')

  # KNN from the first distance matrix
  message("Finding k-nearest neighbors for the embedding")
  nn_mat1 <- MakeNnMatrix(embedding, kfn = kfn, k = k)

  # KNN from the second distance matrix
  message("Finding k-nearest neighbors for original data")
  nn_mat2 <- MakeNnMatrix(orig_data, kfn = kfn, k = k)

  sleepwalk::sleepwalk(embeddings = embedding,
                       compare = "distances",
                       distances = list(nn_mat1, nn_mat2),
                       saveToFile = output_file,
                       pointSize = point_size,
                       titles = plot_names)
}



