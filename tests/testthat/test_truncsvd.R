context("test trunc_svd")

test_that("trunc_svd works", {
  set.seed(35)
  n <- 19
  p <- 23
  r <- 5
  prop_row <- 0.5
  prop_col <- 0.5
  itermax <- 500
  max_rank <- NULL

  A <- matrix(stats::rnorm(n * r), nrow = n)
  B <- matrix(stats::rnorm(r * p), nrow = r)
  E <- matrix(stats::rnorm(n * p), nrow = n)
  Y <- A %*% B + E

  tout <- trunc_svd(Y = Y, prop_row = prop_row, prop_col = prop_col,
                    max_rank = max_rank, itermax = itermax)

  plot(tout)
}
)

test_that("soft_svd works", {
  set.seed(126)
  n <- 19
  p <- 23
  r <- 5
  prop_row <- 0.5
  prop_col <- 0.5
  itermax <- 200
  max_rank <- NULL

  A <- matrix(stats::rnorm(n * r), nrow = n)
  B <- matrix(stats::rnorm(r * p), nrow = r)
  E <- matrix(stats::rnorm(n * p), nrow = n)
  Y <- A %*% B + E

  tout <- soft_svd(Y = Y, prop_row = prop_row, prop_col = prop_col,
                   itermax = itermax)

  plot(tout)
}
)
