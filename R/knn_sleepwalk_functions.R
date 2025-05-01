# Code: KNN sleepwalk functions
# Author: Tyler J Burns
# Date: May 7, 2024

# Edited by David Novak on 1 May, 2025
# -> option to compare preservation of true local structures between embeddings
#    (new `SleepwalkSP` function); this requires exporting `MakeNnMatrix`
# -> option to use RANN for faster k-NN search, so downsampling is not needed
#    for large data (`rann` argument in `MakeNnMatrix`, `KnnSleepwalk`,
#    `SleepwalkSP`)
# -> option to color points by neighbor rank (`ranks` argument in 
#    `MakeNnMatrix`, `KnnSleepwalk`, `SleepwalkSP`)
# -> helper function to plot neighborhood preservation of a single point by an
#    embedding (`PlotNeighborhood`)
# -> helper function to uniformly sample representative points from an
#    embedding on a grid, to plot their neighborhood preservation and compare
#    that to another embedding (`GridSamplePoints`)
# -> `BiaxialSleepwalk` documentation has minimal example code now
# (All new functions have code examples in their documentation.)

#' @importFrom magrittr %>%
#' @import lsa
#' @import sleepwalk
NULL

#' @title Make nearest (or farthest) neighbor matrix
#' @description A helper function for KnnSleepwalk which takes a matrix as
#' input, and generates a nearest neighbor matrix, so that Sleepwalk shows the
#' K-nearest or farthest neighbors as black (or using a gradient by neighbor
#' rank), with the rest of the cells being a lighter color.
#' @param mat A data matrix, presumably from a single-cell dataset
#' @param kfn Whether you want to look at the K-nearest or K-farthest neighbors
#' @param rann Whether to use the faster RANN for neighbor search (only for
#' nearest neighbors and Euclidean distances)
#' @param ranks Whether to retain neighbourhood ranks for a colorscale
#' @param k The number of nearest neighbors
#' @param metric The distance metric to be used. Available options: euclidean,
#' manhattan, and cosine
#' @param ... Any other keyword arguments to pass onto `RANN::nn2` if used
#' @return The aforementioned nearest neighbor matrix
#' @details The default setting computes a distance matrix for `mat`, which is 
#' infeasible for large data. For accelerated nearest-neighbors with Euclidean 
#' distances, set `rann=TRUE`. For further speed-up, also specify a positive 
#' value for the keyword argument `eps` passed onto `RANN::nn2`. This uses an
#' approximate algorithm. Even if `eps` is set to 0, some discrepancies can 
#' occur between the default and `ran=TRUE` setting due to inconsistent 
#' treatment of tied neighbor ranks.
#' @keywords internal
MakeNnMatrix <- function(mat, kfn = FALSE, rann = FALSE, ranks = FALSE, k = 100,
                         metric = "euclidean", ...) {
  metric <- match.arg(metric, choices = c("euclidean", "cosine", "manhattan"))
  if (rann&&(kfn||metric!="euclidean")) {
    stop("RANN is only compatible with nearest-neighbor search using Euclidean
         distance")
  }
  if (rann) { # approximate k-NN with Euclidean distances
    require(RANN)
    
    nn <- RANN::nn2(mat, k = k)$`nn.idx`
    n <- nrow(mat)
    
    nn_mat <- t(apply(
      X = nn,
      MARGIN = 1,
      FUN =
        if (ranks) {
          function(neighb) {
            res <- rep(k+1, times = n)
            res[neighb] <- seq_len(k)
            res
          }
        } else {
          function(neighb) as.numeric(!seq_len(n)%in%neighb)*1000
        }
    ))
    
  } else { # exact k/f-NN via distance matrix
    
    # Create the distance matrix
    if(metric != "cosine") {
      dist_mat <- dist(mat, method = metric) %>% as.matrix()
    } else {
      dist_mat <- lsa::cosine(t(mat))
      dist_mat <- 1 - dist_mat # To make it cosine distance rather than
                               # similarity
    }
    nn_mat <- lapply(seq(nrow(dist_mat)), function(i) {
      curr <- dist_mat[i,]
      working_dist <- sort(curr, decreasing = kfn)[k]
      mask <- curr <= working_dist
      if (ranks) {
        res <- rep(k+1, times = length(curr))
        res[mask] <- rank(curr[mask], ties.method = 'min')
        curr <- res
      } else {
        curr <- as.numeric(!mask)*1000
      }
      return(curr)
    }) %>% do.call(rbind, .)
  }
  
  attributes(nn_mat)$ranks <- ranks
  rownames(nn_mat) <- colnames(nn_mat) <- 1:nrow(nn_mat)
  nn_mat[lower.tri(nn_mat)] = t(nn_mat)[lower.tri(nn_mat)]
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
#' @param k The number of the nearest neighbors to be visualized
#' @param ranks Whether to visualize neighbor ranks using a color gradient
#' @param rann Whether to use faster neighbor search (only for nearest neighbors
#' and not used with cosine/manhattan distances)
#' @param output_file The file to save your sleepwalk html page to
#' @param point_size How big you want the point on the plots to be
#' @param plot_names What you want the comparison plots to be named. This vector
#' of strings needs to match the number of plots (eg. 2 if we're dealing with
#' the default of 2 plots here).
#' @param kfn Whether you want to look at the k-farthest neighbors. If false
#' (default) then you look at the K-nearest neighbors.
#' @param metric The distance metric you're going to use for the neighbor
#' finding computation for the original high-dimensional data. For the embedding,
#' Euclidean is used. Default is also set to euclidean, with manhattan and cosine
#' as other options.
#' @param ... Any other keyword arguments to pass onto `RANN::nn2` if used
#' @export
KnnSleepwalk <- function(embedding,
                         orig_data,
                         k = 100,
                         ranks = FALSE,
                         rann = FALSE,
                         output_file = NULL,
                         point_size = 1.5,
                         plot_names = c("KNN embedding space", "KNN high-dim space"),
                         kfn = FALSE,
                         metric = "euclidean",
                         ...) {
  
  message('Building distance matrix')
  
  # KNN from the first distance matrix
  message("Finding k-",
          if (kfn) "farthest" else "nearest",
          " neighbors for the embedding")
  nn_mat1 <- MakeNnMatrix(embedding, kfn = kfn, k = k, metric = "euclidean",
                          ranks = ranks, rann = rann)
  
  # KNN from the second distance matrix
  message("Finding k-",
          if (kfn) "farthest" else "nearest",
          " neighbors for original data")
  nn_mat2 <- MakeNnMatrix(orig_data, kfn = kfn, k = k, metric = metric,
                          ranks = ranks, rann = (rann&&metric=="euclidean"))
  
  sleepwalk::sleepwalk(embeddings = embedding,
                       compare = "distances",
                       distances = list(nn_mat1, nn_mat2),
                       saveToFile = output_file,
                       pointSize = point_size,
                       titles = plot_names)
}

#' @title K-NN Sleepwalk view of Local Structure Preservation
#' @description Assesses local structure preservation (SP) in a 2-D embedding by
#' projecting K-nearest-neighbor ranks from the original data onto it. If
#' multiple embeddings are specified, it compares them side-by-side with respect
#' to point-wise Local SP.
#' @param embeddings Single embedding in matrix format with data points
#' as rows and two columns, or a list of multiple embeddings
#' @param plot_names What you want the comparison plots to be named. This vector
#' of strings needs to match the number of plots (eg. 2 if we're dealing with
#' 2 plots here).
#' @param nn Nearest-neighbor matrix of original data generated using 
#' `KnnSleepwalk::MakeNnMatrix` with `ranks=TRUE`
#' @param output_file The file to save your sleepwalk html page to
#' @param point_size How big you want the point on the plots to be
#' @seealso [KnnSleepwalk::MakeNnMatrix] to pre-compute the neighbor matrix
#' @examples
#' ## Load example data & UMAP embedding
#' data("example_umap") # `umap`
#' data("example_surface_markers") # `surface`
#' 
#' ## Generate alternative PCA embedding
#' pca <- prcomp(surface)$x[, 1:2]
#'
#' ## Compute neighbor matrix
#' nn <- MakeNnMatrix(surface, k = 100, ranks = TRUE, rann = TRUE)
#' 
#' ## Compare structure preservation by PCA & UMAP
#' SleepwalkSP(
#'     embeddings = list(pca, umap),
#'     plot_names = c("PCA", "UMAP"),
#'     nn = nn,
#'     point_size = 3
#' )
#' @export
SleepwalkSP <- function(embeddings,
                        nn,
                        plot_names = paste0("Emb", seq_along(embeddings)),
                        output_file = NULL,
                        point_size = 1.5,
                        metric = "euclidean") {
  
  if (is.null(attributes(nn)$ranks) || !attributes(nn)$ranks) {
    stop("`nn` must be generated using `KnnSleepwalk::MakeNnMatrix` with ",
         "`ranks` set to TRUE")
  }
  if (is.matrix(embeddings)) {
    embeddings <- list(embeddings)
  }
  if (length(unique(sapply(embeddings, nrow)))>1) {
    stop("Embeddings must have the number of rows")
  }
  if (any(sapply(embeddings, ncol)!=2)) {
    stop("Embeddings must have 2 columns")
  }
  
  if (is.null(nn)) {
    message("Finding k-",
            if (kfn) "farthest" else "nearest",
            " neighbors for original data")
    nn <- MakeNnMatrix(orig_data, kfn = kfn, k = k, metric = metric)
  }
  
  sleepwalk::sleepwalk(embeddings = embeddings,
                       same = "objects",
                       compare = "embeddings",
                       maxdists = max(nn)-1,
                       distances = nn,
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
#' @examples
#' data("example_umap")
#' data("example_surface_markers")
#' BiaxialSleepwalk(
#'   root_biax = surface[,c("CD4", "CD8")],
#'   biax_list = list(
#'     surface[,c("CD4", "CD8")],
#'     surface[,c("IgA", "IgM")],
#'     surface[,c("IgD", "HLA-DR")]
#'   )
#' )
#' @export
BiaxialSleepwalk <- function(root_biax,
                             biax_list,
                             k = 100,
                             output_file = NULL,
                             point_size = 1.5,
                             plot_names = c()) {
  
  nn_mat <- MakeNnMatrix(mat = root_biax, k = k, kfn = FALSE)
  
  sleepwalk::sleepwalk(embeddings = biax_list,
                       distances = nn_mat,
                       pointSize = point_size,
                       saveToFile = output_file,
                       titles = plot_names)
}

#' @title Plot point neighborhood in an embedding
#' @description The neighborhood of a single point is shown in an embedding.
#' This is useful for visualizing the preservation of local structure and
#' comparing that between embeddings. To that end, the neighborhood matrix
#' should be computed for the original data.
#' @param embedding A data matrix of points in a 2-D embedding, with points as
#' rows
#' @param point_index Row index of a single point within `embedding`
#' @param nn Neighborhood matrix computed using `KnnSleepwalk::MakeNnMatrix` 
#' with `ranks=TRUE`
#' @param point_size size of points in the embedding plot
#' @param preserve_aspect_ratio whether to fix a 1-to-1 aspect ratio between the
#' latent-space components
#' @returns A `ggplot` object
#' @seealso [KnnSleepwalk::GridSamplePoints] to sample representative points for
#' use in the `point_index`
#' @examples
#' ## Load example embedding & convert to matrix
#' data("example_umap") # `umap`
#' umap <- as.matrix(umap)
#' 
#' ## Sample representative points
#' pts <- GridSamplePoints(umap)
#' 
#' ## Show these points in a plot
#' plot(umap, pch = 20, cex = .5, col = "gray", asp = 1,
#'      xlab = "Component 1", ylab = "Component 2")
#' points(umap[pts, 1], umap[pts, 2],
#'        pch = 10, cex = 2, col = "darkred")
#' 
#' ## Load the original data
#' data("example_surface_markers") # `surface`
#' 
#' ## Compute local neighborhoods in high dimension
#' nn <- MakeNnMatrix(surface, k = 100, ranks = TRUE, rann = TRUE)
#' 
#' ## Create composite PNG of plots of neighborhoods for the sample points
#' neighb_plots <- lapply(pts, function(pt) 
#'   PlotNeighborhood(embedding = umap, point_index = pt, nn = nn) +
#'     theme(legend.position = "none"))
#' figure <- cowplot::plot_grid(plotlist = neighb_plots)
#' ggplot2::ggsave("./neighbors.png", figure,
#'                 width = 30, height = 30, dpi = 300)
#' 
#' @export
PlotNeighborhood <- function(
    embedding,
    point_index,
    nn,
    point_size = 1.
) {
  
  stopifnot(isTRUE(attributes(nn)$ranks))
  require(tidyverse)
  require(cowplot)
  
  d <- embedding %>%
    as_tibble() %>% 
    `colnames<-`(c('Component 1', 'Component 2')) %>% 
    mutate('Neighbor rank' =
             na_if(nn[point_index,], max(nn))) # out-of-bounds as NA~gray
  
  ggplot(d, aes(
    x = `Component 1`, y = `Component 2`, color = `Neighbor rank`
  )) +
    geom_point( # background points
      data = dplyr::filter(d, is.na(`Neighbor rank`)),
      size = point_size
    ) + 
    geom_point( # neighbors
      data = dplyr::filter(d, !is.na(`Neighbor rank`)),
      size = point_size
    ) + 
    geom_point( # original point
      data = dplyr::filter(d, row_number() == point_index),
      color = 'orange',
      shape = 10,
      size = point_size*1.5
    ) + 
    scale_color_continuous(na.value = '#cccccc') +
    coord_fixed()
}

#' @title Sample representative points on a grid
#' @description Returns indices of representative points within a 2-D
#' embedding, sampled uniformly on a `xdim`-by-`ydim` grid. Points closest to 
#' the center of each non-empty tile are selected.
#' @param mat A data matrix of points in a 2-D embedding, with points as
#' rows
#' @param xdim Grid width
#' @param ydim Grid height
#' @returns Numeric vector
#' @seealso [KnnSleepwalk::PlotNeighborhood] to plot local neighborhoods of 
#' @examples
#' ## Load example embedding & convert to matrix
#' data("example_umap") # `umap`
#' umap <- as.matrix(umap)
#' 
#' ## Sample representative points
#' pts <- GridSamplePoints(umap)
#' 
#' ## Show these points in a plot
#' plot(umap, pch = 20, cex = .5, col = "gray", asp = 1,
#'      xlab = "Component 1", ylab = "Component 2")
#' points(umap[pts, 1], umap[pts, 2],
#'        pch = 10, cex = 2, col = "darkred")
#' 
#' ## Load the original data
#' data("example_surface_markers") # `surface`
#' 
#' ## Compute local neighborhoods in high dimension
#' nn <- MakeNnMatrix(surface, k = 100, ranks = TRUE, rann = TRUE)
#' 
#' ## Create composite PNG of plots of neighborhoods for the sample points
#' neighb_plots <- lapply(pts, function(pt) 
#'   PlotNeighborhood(embedding = umap, point_index = pt, nn = nn) +
#'     theme(legend.position = "none"))
#' figure <- cowplot::plot_grid(plotlist = neighb_plots)
#' ggplot2::ggsave("./neighbors.png", figure,
#'                 width = 30, height = 30, dpi = 300)
#' @export
GridSamplePoints <- function(
    mat,
    xdim = 10,
    ydim = 10
) {
  
  ## Create grid of points
  x <- seq(min(mat[, 1]), max(mat[, 1]), length.out = xdim+1)
  y <- seq(min(mat[, 2]), max(mat[, 2]), length.out = ydim+1)
  diffx <- diff(x)[1] # stepsize
  diffy <- diff(y)[1]
  x <- head(x, xdim) # trim edge
  y <- head(y, ydim)
  grid <- as.matrix(expand.grid(x, y))
  
  ## Sample point per tile
  nearest_points <- unlist(
    apply(
      grid, 1, function(g) {
        
        ## Find points within tile
        x1 <- g[1]; x2 <- g[1]+diffx; y1 <- g[2]; y2 <- g[2]+diffy
        mask <- mat[,1]>x1 & mat[,2]>y1 & mat[,1]<=x2 & mat[,2]<=y2
        if (sum(mask)==0) NULL
        
        ## Find point closest to tile center
        dists <- (mat[mask, 1]-x1-diffx/2)^2 + (mat[mask, 2]-y1-diffy/2)^2
        which(mask)[which.min(dists)]
      })
  )
  return(nearest_points)
}
