library(testthat)
library(KnnSleepwalk)

# Load the datasets
data("example_surface_markers", package = "KnnSleepwalk")
data("example_umap", package = "KnnSleepwalk")

# TODO
test_that("KnnSleepwalk runs with valid inputs", {
  expect_silent(2 + 2)
})
