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
plot.swumobj <- function(x, ...) {
  assertthat::assert_that(class(x) == "swumobj")

  ## All objects will have the empirical sv's --------------------------------
  dat <- data.frame(Value = x$d, SV = as.factor(1:length(x$d)), stringsAsFactors = FALSE)
  pl2 <- ggplot2::ggplot(data = dat, mapping = ggplot2::aes_string(x = "SV")) +
    ggplot2::geom_linerange(mapping = ggplot2::aes_string(ymax = "Value"), ymin = 0) +
    ggplot2::theme_bw() +
    ggplot2::ylim(0, max(x$d)) +
    ggplot2::ggtitle("Scree Plot")
  print(pl2)

  ## Rank distribution -------------------------------------------------------
  if (!is.null(x$rank_dist)) {
    cat ("Press [enter] to continue")
    line <- readline()
    dat <- data.frame(Probability = x$rank_dist,
                      Rank = as.factor(as.numeric(names(x$rank_dist))),
                      stringsAsFactors = FALSE)
    pl1 <- ggplot2::ggplot(data = dat, mapping = ggplot2::aes_string(x = "Rank")) +
      ggplot2::geom_linerange(mapping = ggplot2::aes_string(ymax = "Probability"), ymin = 0) +
      ggplot2::theme_bw() +
      ggplot2::ylim(0, max(x$rank_dist)) +
      ggplot2::ggtitle("Rank Distribution")
    print(pl1)
  }

  if (!is.null(x$dhat)) {
    cat ("Press [enter] to continue")
    line <- readline()
    dat <- data.frame(Empirical = x$d, Estimated = x$dhat, SV = as.factor(1:length(x$d)),
                      stringsAsFactors = FALSE)
    longdat <- tidyr::gather(data = dat, key = "Type", value = "Value", 1:2)
    pl3 <- ggplot2::ggplot(data = longdat, mapping = ggplot2::aes_string(x = "SV", lty = "Type")) +
      ggplot2::geom_linerange(mapping = ggplot2::aes_string(ymax = "Value"),
                              ymin = 0, position = ggplot2::position_dodge(0.4)) +
      ggplot2::ylim(0, max(x$d)) +
      ggplot2::ggtitle("Empirical and Estimated SV's") +
      ggplot2::theme_bw()
    print(pl3)
  }

  if (!is.null(x$dlower) & !is.null(x$dupper)) {
    cat ("Press [enter] to continue")
    line <- readline()
    dat2 <- data.frame(Lower = x$dlower, Upper = x$dupper, Estimated = x$dhat,
                       Empirical = x$d,
                       SV = as.factor(1:length(x$dhat)),
                       stringsAsFactors = FALSE)
    longdat <- tidyr::gather(data = dat2, key = "Type", value = "Value", 3:4)
    pl3 <- ggplot2::ggplot(data = longdat,
                           mapping = ggplot2::aes_string(x = "SV", y = "Value", pch = "Type")) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(mapping = ggplot2::aes_string(ymin = "Lower", ymax = "Upper"),
                             width = 0.2, lty = 2) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Value") +
      ggplot2::ylim(0, max(x$d)) +
      ggplot2::ggtitle("Point and Interval Estimates of SV's")
    print(pl3)
  }







}