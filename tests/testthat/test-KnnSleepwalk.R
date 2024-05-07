test_that("KnnSleepwalk runs with valid inputs", {
  mat1 <- matrix(rnorm(1000), ncol = 10)
  mat2 <- matrix(rnorm(1000), ncol = 10)
  embedding <- matrix(rnorm(200), ncol = 2)
  expect_silent(KnnSleepwalk(mat1, mat2, embedding, k = 5))
})
