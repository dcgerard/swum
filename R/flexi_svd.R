## Flexible SVD

#' Hyper-flexible SVD.
#'
#' @inheritParams trunc_svd
#'
#' @author David Gerard
#'
flexi_svd <- function(Y, prop_row = 0.5, prop_col = 0.5,
                      itermax = 500, max_rank = NULL) {
  ## Test Input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(Y))
  assertthat::assert_that(prop_row > 0, prop_row < 1)
  assertthat::assert_that(prop_col > 0, prop_col < 1)
  assertthat::assert_that(itermax > 0)

  ## Get Dimensions ----------------------------------------------------------
  n <- nrow(Y)
  p <- ncol(Y)
  n1 <- round(n * prop_row)
  n2 <- n - n1
  p1 <- round(p * prop_col)
  p2 <- p - p1

  assertthat::assert_that(n1 > 0, n2 > 0, p1 > 0, p2 > 0)

  if (is.null(max_rank)) {
    max_rank <- min(c(n1, n2, p1, p2))
  } else {
    assertthat::assert_that(max_rank <= n1, max_rank <= n2, max_rank <= p1, max_rank <= p2)
  }

  d_samp <- matrix(NA, nrow = itermax, ncol = max_rank)
  for (index in 1:itermax) {
    which_row <- rep(FALSE, length = n)
    which_col <- rep(FALSE, length = p)
    which_row[sample(x = 1:n, size = n1)] <- TRUE
    which_col[sample(x = 1:p, size = p1)] <- TRUE

    Y11 <- Y[which_row, which_col]
    Y12 <- Y[which_row, !which_col]
    Y21 <- Y[!which_row, which_col]
    Y22 <- Y[!which_row, !which_col]

    svdY11 <- svd(Y11)

    Y21V <- Y21 %*% svdY11$v
    UtY12 <- crossprod(svdY11$u, Y12)

    A <- Y21V[, 1:max_rank]
    B <- UtY12[1:max_rank, ]

    svhat <- sort(abs(1 / diag(solve(crossprod(A)) %*% A %*% Y22 %*% crossprod(B, solve(tcrossprod(B))))), decreasing = TRUE)

    d_samp[index, ] <- svhat
  }

  apply(d_samp, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
  colMeans(d_samp)

}
