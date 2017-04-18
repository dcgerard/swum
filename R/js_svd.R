
#' Josse-Sardy svd.
#'
#' @inheritParams trunc_svd
#'
#' @author David Gerard
#'
#' @export
#'
#' @references Josse, Julie, and Sylvain Sardy. "Adaptive shrinkage of singular values." Statistics and Computing 26.3 (2016): 715-724.
#'
js_svd <- function(Y, prop_row = 0.5, prop_col = 0.5,
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
  gamma_samp <- rep(NA, length = itermax)

  gamma_seq <- 2 ^ seq(0, 8, length = 20)
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

    lambda_seq <- c(svdY11$d, 0)

    lambda_best <- lambda_seq[1]
    gamma_best <- gamma_seq[1]
    err_best <- js_obj(par = c(lambda_best, gamma_best),
                       Y22 = Y22, Y21V = Y21V, UtY12 = UtY12,
                       d = svdY11$d)
    for (lambda_index in 1:length(lambda_seq)) {
      for (gamma_index in 1:length(gamma_seq)) {
        err <- js_obj(par = c(lambda_seq[lambda_index], gamma_seq[gamma_index]),
                      Y22 = Y22, Y21V = Y21V, UtY12 = UtY12,
                      d = svdY11$d)
        if (err < err_best) {
          err_best <- err
          gamma_best <- gamma_seq[gamma_index]
          lambda_best <- lambda_seq[lambda_index]
        }
      }
    }

    # oout <- stats::optim(par = c(lambda_best, gamma_best), fn = js_obj,
    #                      Y22 = Y22, Y21V = Y21V, UtY12 = UtY12,
    #                      d = svdY11$d)

    lambda_samp[index] <- lambda_best
    gamma_samp[index] <- gamma_best
  }

  final_d <- svd(x = Y, nu = 0, nv = 0)$d

  sample_d <- pmax(t(1 - outer(lambda_samp, final_d, FUN = `/`) ^ gamma_samp) * final_d, 0)
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
#' @param par The tuning parameters. The first is the shrinkage parameter.
#'     The other is the magnitude parameter.
#' @inheritParams soft_obj
#'
#' @seealso \code{\link{soft_svd}}
#'
#' @author David Gerard
#'
#' @return The sums of squares difference from Y22 of the estimate of Y22.
js_obj <- function(par, Y22, Y21V, UtY12, d) {

  ## Check Input -------------------------------------------------------------
  assertthat::are_equal(nrow(Y22), nrow(Y21V))
  assertthat::are_equal(ncol(Y22), ncol(UtY12))
  assertthat::are_equal(ncol(Y21V), nrow(UtY12))
  assertthat::are_equal(length(d), nrow(UtY12))

  lambda <- par[1]
  gamma <- par[2]
  if (lambda < 0 | gamma < 0) {
    return(Inf) ## Only allow for shrinkage
  }


  new_d <- d * (1 - (lambda / d) ^ gamma)
  current_rank <- sum(new_d > 0)

  if (current_rank == 0) {
    return(mean(Y22 ^ 2))
  }

  D_mat <- diag(x = 1 / new_d[new_d > 0], nrow = current_rank, ncol = current_rank)
  err <- mean((Y22 - Y21V[, 1:current_rank] %*% D_mat %*% UtY12[1:current_rank, ]) ^ 2)
  return(err)
}

