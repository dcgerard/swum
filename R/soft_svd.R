
#' Soft-thresholding SVD estimator.
#'
#' @inheritParams trunc_svd
#'
#' @author David Gerard
#'
#' @export
#'
#' @return A list with the following elements
#' \itemize{
#'   \item{"rank_dist"}{The distribution of ranks.}
#'   \item{"d"}{The empirical singular values.}
#'   \item{"dhat"}{The estimated singular values}
#'   \item{"dlower"}{The lower bound for the 95\% interval}
#'   \item{"dupper"}{The upper bound for the 95\% interval}
#'   \item{"lambda_samp"}{The samples of the shrinkage parameter.}
#' }
#'
#'
soft_svd <- function(Y, prop_row = 0.5, prop_col = 0.5,
                     itermax = 500) {
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

  lambda_samp <- rep(NA, length = itermax)
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

    oout <- stats::optim(par = 1, fn = soft_obj, method = "Brent", lower = 0,
                         upper = max(svdY11$d), Y22 = Y22, Y21V = Y21V, UtY12 = UtY12,
                         d = svdY11$d)

    lambda_samp[index] <- oout$par
  }

  ## Summaries
  final_d <- svd(x = Y, nu = 0, nv = 0)$d
  sample_d <- outer(final_d, lambda_samp, FUN = `-`)
  sample_d[sample_d < 0] <- 0
  rank_dist <- table(colSums(sample_d > 10 ^ -6))
  possible_ranks <- names(rank_dist)
  rank_dist <- as.numeric(rank_dist / sum(rank_dist))
  names(rank_dist) <- possible_ranks
  dhat <- apply(sample_d, 1, stats::median)
  dlower <- apply(sample_d, 1, stats::quantile, probs = 0.025)
  dupper <- apply(sample_d, 1, stats::quantile, probs = 0.975)

  outlist <- list(rank_dist = rank_dist, d = final_d, dhat = dhat, dlower = dlower,
                  dupper = dupper, lambda_samp = lambda_samp)
  class(outlist) <- "swumobj"
  return(outlist)
}

#' The objective function to minimize in \code{\link{soft_svd}}.
#'
#' @param lambda The tuning parameter.
#' @param Y22 The lower right matrix.
#' @param Y21V The estimate of U2
#' @param UtY12 The estimate of V2^t
#' @param d The singular values of Y11.
#'
#' @seealso \code{\link{soft_svd}}
#'
#' @author David Gerard
#'
#' @return The sums of squares difference from Y22 of the estimate of Y22.
soft_obj <- function(lambda, Y22, Y21V, UtY12, d) {

  ## Check Input -------------------------------------------------------------
  assertthat::are_equal(nrow(Y22), nrow(Y21V))
  assertthat::are_equal(ncol(Y22), ncol(UtY12))
  assertthat::are_equal(ncol(Y21V), nrow(UtY12))
  assertthat::are_equal(length(d), nrow(UtY12))
  if (lambda < 0) {
    return(Inf) ## Only allow for shrinkage
  }

  new_d <- d - lambda
  current_rank <- sum(new_d > 0)
  D_mat <- diag(x = 1 / new_d[new_d > 0], nrow = current_rank, ncol = current_rank)
  err <- mean((Y22 - Y21V[, 1:current_rank] %*% D_mat %*% UtY12[1:current_rank, ]) ^ 2)
  return(err)
}




