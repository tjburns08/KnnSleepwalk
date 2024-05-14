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
#' @param metric The distance metric to be used. Available options: euclidean,
#' manhattan, and cosine
#' @return The aforementioned nearest neighbor matrix
#' @keywords internal
MakeNnMatrix <- function(mat, kfn = FALSE, k = 100, metric = "euclidean") {
  # Create the distance matrix
  if(metric != "cosine") {
    dist_mat <- dist(mat, method = metric) %>% as.matrix()
  } else {
    dist_mat <- lsa::cosine(t(mat))
    dist_mat <- 1 - dist_mat # To make it cosine distance rather than similarity

    # We have to re-name the rows and columns
    rownames(dist_mat) <- colnames(dist_mat) <- 1:nrow(dist_mat)
  }

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
#' @param kfn Whether you want to look at the k-farthest neighbors. If false
#' (default) then you look at the K-nearest neighbors.
#' @param metric The distance metric you're going to use for the neighbor
#' finding computation for the original high-dimensional data. For the embedding,
#' Euclidean is used. Default is also set to euclidean, with manhattan and cosine
#' as other options.
#' @export
KnnSleepwalk <- function(embedding,
                         orig_data,
                         k = 100,
                         output_file = NULL,
                         point_size = 1.5,
                         plot_names = c("KNN embedding space", "KNN high-dim space"),
                         kfn = FALSE,
                         metric = "euclidean") {

  message('Building distance matrix')

  # KNN from the first distance matrix
  message("Finding k-nearest neighbors for the embedding")
  nn_mat1 <- MakeNnMatrix(embedding, kfn = kfn, k = k, metric = "euclidean")

  # KNN from the second distance matrix
  message("Finding k-nearest neighbors for original data")
  nn_mat2 <- MakeNnMatrix(orig_data, kfn = kfn, k = k, metric = metric)

  sleepwalk::sleepwalk(embeddings = embedding,
                       compare = "distances",
                       distances = list(nn_mat1, nn_mat2),
                       saveToFile = output_file,
                       pointSize = point_size,
                       titles = plot_names)
}

#' @title Biaxial Sleepwalk
#' @description A generalization of KNN sleepwalk, whereby the K-nearest
#' neighbors of the first biaxial plot are calculated, and the corresponding
#' cells are highlighted on any of the other biaxials that are on the map
#' @param root_biax The biaxial plot from which the KNN are calculated
#' @param biax_list A list that starts with root_biax and continues with the
#' remaining biaxials that you care about. These will be visualized.
#' @param k The number of nearest neighbors
#' @param output_file The file that the output gets saved to. If set to NULL,
#' then no save takes place
#' @param point_size The size of the points on the plots.
#' @param plot_names These will place names on top of the plots. This vector of
#' strings must equal the total number of plots being visualized.
#' @export
BiaxialSleepwalk <- function(root_biax,
                             biax_list,
                             k = 100,
                             output_file = NULL,
                             point_size = 1.5,
                             plot_names = c()) {

  dist_mat <- dist(root_biax) %>% as.matrix()
  nn_mat <- MakeNnMatrix(mat = dist_mat, k = k, kfn = FALSE)

  sleepwalk::sleepwalk(embeddings = biax_list,
                       distances = nn_mat,
                       pointSize = point_size,
                       saveToFile = output_file,
                       titles = plot_names)
}
