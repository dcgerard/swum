## plotting functions

#' Plots output from \code{\link{trunc_svd}}.
#'
#' @param x The output of \code{\link{trunc_svd}}.
#' @param ... Not used.
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso \code{\link{trunc_svd}}
plot.truncsvd <- function(x, ...) {
  assertthat::assert_that(class(x) == "truncsvd")

  dat <- data.frame(Probability = x$rank_dist, Rank = as.factor(0:(length(x$rank_dist) - 1)))

  pl1 <- ggplot2::ggplot(data = dat, mapping = ggplot2::aes_string(x = "Rank")) +
    ggplot2::geom_linerange(mapping = ggplot2::aes_string(ymax = "Probability"), ymin = 0) +
    ggplot2::theme_bw() +
    ggplot2::ylim(0, max(x$rank_dist)) +
    ggplot2::ggtitle("Rank Distribution")

  print(pl1)

  cat ("Press [enter] to continue")
  line <- readline()

  dat <- data.frame(Value = x$d, SV = as.factor(1:length(x$d)))

  pl2 <- ggplot2::ggplot(data = dat, mapping = ggplot2::aes_string(x = "SV")) +
    ggplot2::geom_linerange(mapping = ggplot2::aes_string(ymax = "Value"), ymin = 0) +
    ggplot2::theme_bw() +
    ggplot2::ylim(0, max(x$d)) +
    ggplot2::ggtitle("Scree Plot")

  print(pl2)
}