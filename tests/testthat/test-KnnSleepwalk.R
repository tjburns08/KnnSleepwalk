library(testthat)
library(KnnSleepwalk)

# Load the datasets
data("example_surface_markers", package = "KnnSleepwalk")
data("example_umap", package = "KnnSleepwalk")

test_that("A matrix is produced in MakeNnMatrix", {
  testthat::expect_true(is.matrix(KnnSleepwalk:::MakeNnMatrix(umap)))
})

test_that("The knn matrix is thresholded at 2 values", {
  testthat::expect_true(length(table(KnnSleepwalk:::MakeNnMatrix(umap))) == 2)
})
