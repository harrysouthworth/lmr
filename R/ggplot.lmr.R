#' QQ-plot for residuals from a model using ggplot2
#'
#' @param o A fitted model with a \code{resid} component.
#' @return a \code{ggplot2} plot representing the Gaussian QQ-plot of the residuals
#'           with a fitted line.
#; @keywords hplot
#' @export ggqqplot
ggqqplot <- function(o) {
  y <- quantile(resid(o)[!is.na(resid(o))], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]

  d <- data.frame(r=resid(o))

  p <- ggplot(d, aes(sample = r)) +
    stat_qq(alpha = 0.5, color="blue") +
    geom_abline(slope = slope, intercept = int, color="orange")

  return(p)
}

#' ggplot for rlm results
#' @aliases ggplot.lmr
#' @param data an object of class \code{rlm}.
#' @param hist.scale A scaleing factor for bin widths in the histogram. The default bin width
#'        will be \code{range(resid(data))/hist.scale)}. Defaults to \code{hist.scale=10},
#'        which is smaler than the \code{ggplot} default of 30.
#' @param plot. Logical indicating whether to plot. If FALSE, the function
#'   returns a list of ggplot objects that can be passed to \code{gridExtra::grid.arrange}.
#' @param caption A caption to appear at the bottom left of the plot
#' @param ... Additional arguments passed to \code{ggplot}. Currently unused.
#' @importFrom gridExtra grid.arrange
#' @method plot lmr
#' @export
plot.lmr <- function(x=NULL, hist.scale=10, plot.=TRUE, caption = NULL, ...){
  data <- x
  d <- data.frame(r = resid(data), f=fitted(data), o=resid(data) + fitted(data))

  qq <- ggqqplot(data)

  bin <- diff(range(d$r, na.rm=TRUE) / hist.scale)
  hist <- ggplot() +
    geom_histogram(data=d, aes(r, after_stat(density)), fill="blue", binwidth=bin) +
    geom_density(data=d, aes(r, after_stat(density)), color="light blue" ) +
    geom_rug(data=d, aes(r), color="orange" ) +
    scale_x_continuous("Residuals")

  fr <- ggplot(d, aes(f, r)) +
    geom_point(color="blue") +
    stat_smooth(method="loess", color="orange", formula = y ~ x,
                method.args=list(span=2/3, family="symmetric", degree=1)) +
    scale_x_continuous("Fitted values") +
    scale_y_continuous("Residuals")

  fo <-  ggplot(d, aes(f, o)) +
    geom_point(color="blue") +
    stat_smooth(method="loess", color="orange", formula = y ~ x,
                method.args=list(span=2/3, family="symmetric", degree=1)) +
    scale_x_continuous("Fitted values") +
    scale_y_continuous("Observed values")

  if (!is.null(caption)){
    fo <- fo + labs(caption = caption)
  }

  if (plot.){
    gridExtra::grid.arrange(qq, hist, fr, fo)
  } else {
    list(qq, hist, fr, fo)
  }
}

#' ggplot method
#' @method ggplot lmr
#' @export
ggplot.lmr <- plot.lmr
