
#' Truncated SVD with uncertainty measures.
#'
#' @param Y A matrix of data. There can be no missing values (yet).
#' @param prop_row The proportion of rows in Y11.
#' @param prop_col The proportion of columns in Y11.
#' @param max_rank The maximum rank to consider.
#' @param itermax The maximum number of repetitions to run.
#'
#' @author David Gerard
#'
#' @export
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item{"rank_dist"}{The distribution of ranks.}
#'   \item{"d"}{The empirical singular values.}
#' }
#'
trunc_svd <- function(Y, prop_row = 0.5, prop_col = 0.5,
                      max_rank = NULL, itermax = 500) {
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
    max_rank <- min(n1, p1)
  } else {
    assertthat::assert_that(n1 >= max_rank, p1 >= max_rank)
  }

  ## Run Algorithm -----------------------------------------------------------
  rank_dist <- rep(0, length = max_rank + 1)
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

    error_vec <- rep(NA, length = max_rank + 1)
    error_vec[1] <- mean(Y22 ^ 2)
    for (trunc_index in 1:max_rank) {
      new_d <- svdY11$d[1:trunc_index]
      error_vec[trunc_index + 1] <-
        mean((Y22 - Y21V[, 1:trunc_index, drop = FALSE] %*%
                diag(x = 1 / new_d, nrow = trunc_index, ncol = trunc_index) %*%
                UtY12[1:trunc_index, ]) ^ 2)
    }

    new_rank_loc <- which.min(error_vec)
    rank_dist[new_rank_loc] <- rank_dist[new_rank_loc] + 1
  }

  final_d <- svd(x = Y, nu = 0, nv = 0)$d

  rank_dist <- rank_dist / sum(rank_dist)
  names(rank_dist) <- 0:max_rank
  outlist <- list(rank_dist = rank_dist, d = final_d)
  class(outlist) <- "swumobj"
  return(outlist)
}