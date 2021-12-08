## Requirements:
##   Use MASS::rlm, robustbase::lmrob;
##   Allow initial estimates from Peno-Yohai;
##   Enable RFPE;
##   Have summary function not produce pages of crap;
##   Consistent, nice defaults (bisquare, 85% efficiency);
##   Consistent, non-whiney behaviour of predict;
##   Consistent computation of p-values;
##   Analysis of variance type function.
##   Fast bootstrap;

##   Consistent R^2 returned;

## Nice to haves:
##   Use robust::lmRob
##   Easy switch options to RobStatTM2019 proposals.

#' Pena-Yohai regression estimates
#' @details This is a simple wrapper to \code{pyinit::pyinit} that imposes
#'   some defaults and provides a formula/data interface.
#' @param formula A formula describing a linear model.
#' @param data A data frame.
#' @param resid_keep_method,psc_keep,resid_keep_prop,... All passed through to
#'   \code{pyinit::pyinit}. Basically, this just provides defaults that
#'   \code{pyinit::pyinit} insists the user decides upon.
#' @export
py <- function(formula, data, resid_keep_method = "proportion",
               psc_keep = .5, resid_keep_prop = .5, ...){
  mf <- model.frame(formula, data) %>% droplevels()
  mm <- model.matrix(formula, droplevels(data))
  resp <- model.response(mf)

  ## Whether an intercept is present is determined by the formula
  inits <- pyinit::pyinit(x = mm, y = resp, intercept = FALSE,
                          resid_keep_method = resid_keep_method,
                          psc_keep = psc_keep,
                          resid_keep_prop = resid_keep_prop, ...)

  list(coefficients = inits$coefficients[, 1], scale = inits$objective[1])
}

#' Robust final predcition error
#' @details This is a modified copy-paste of \code{RobStatTM::lmrobdetMM.RFPE}.
#' @param object An object of class 'lmr'.
#' @param scale The scale estimate. This must be provided: RFPE is not comparable
#'   across models with different scale estimates. The scale estimate should be
#'   the one from the model with the full set of predictors.
#' @export
rfpe <- function(object, scale){
  if (!object$converged){
    warning("The algorithm did not converge, inference is not recommended.")
  }

  if (missing(scale)){
    stop("RFPE only works if the same scale estimate is used for each model. Therefore,
         you must provide the scale.")
  }

  p <- length(object$coef)

  res <- residuals(object) / scale

  if (inherits(object, "lmrob")){
    fam <- object$control$psi
    cnst <- object$control$tuning.psi
  } else if (inherits(object, "rlm")){
    fam <- object$psi
    cnst <- object$cc
  } else {
    stop("lmr ought to be fitting using lmrob or rlm")
  }

  a2 <- mean(RobStatTM::rho(u = res, family = fam, cc = cnst, standardize=TRUE))
  b2 <- p * mean(RobStatTM::rhoprime(u=res, family = fam, cc = cnst, standardize=TRUE)^2)
  d2 <- mean(RobStatTM::rhoprime2(u=res, family = fam, cc = cnst, standardize=TRUE))

  if (d2 <= 0){
    NA
  } else {
    a2 + b2/d2/length(res)
  }
}

#' Fast robust bootstrap
#' @details This is a simple wrapper for \code{FRB::frb}.
#' @param lmrob.object,nboot,return.coef As for \code{FRB::frb}.
#' @export
frb <- function(lmrob.object, nboot = 1000, return.coef = FALSE){
  if (inherits(lmrob.object, "rlm")){
    stop("engine must be lmrob, not rlm")
  }

  FRB::frb(lmrob.object, nboot, return.coef)
}

#' Robust regression using MM-estimation
#'
#' Robust regression using MM-estimation with 85\% efficiency for Gaussian data.
#'
#' @param formula A formula describing a linear model.
#' @param data An appropriate data frame.
#' @param weights Not used. This is only here because \code{ggplot2::geom_smooth}
#'   appears to require any custom smoother to take the argument.
#' @param psi The psi function to use. The default is to use bisquare weight functions.
#'   Note that if \code{engine = "rlm"} then psi is set to \code{MASS::psi.bisquare}
#'   (i.e. a function) and if \code{engine = "lmrob"}, psi is set to \code{"bisqare"}
#'   (i.e. a character string). See the help files for \code{rlm} and \code{lmrob}
#'   for further information on available values for psi.
#' @param method The robust fitting method. Defaults to \code{method="MM"}.
#' @param c Tuning parameter to the MM-algorithm. Defaults to \code{c=3.443689} giving 85\% efficiency for Gaussian data.
#' @param engine Character string specifying either 'rlm' in which case \code{MASS::rlm} is used,
#'   or 'lmrob' in which case \code{robustbase::lmrob} is used. In the latter case, a robust
#'   version of R^2 is provided, but the default output produces p-values based on t-distributions
#'   that have no theoretical justification. Bootstrapping would be much better.
#' @param maxit The maximum number of iterations to perform. Defaults to \code{maxit = 40}
#' @param ... Other arguments to be passed to \code{rlm} (like contrasts)
#' @return An object of class 'rlm', fit by the function in the MASS package.
#'         \code{lmr} is just a simple wrapper to \code{rlm}. The returned
#'         object has an additional component, \code{cov}.
#' @return A list containing the same elements as an object of class 'lmr' but with
#'        additional elements containing the covariance matrix of the parameter
#'        estimates ('cov'), the data provided in the call to \code{lmr} ('data'),
#'        the tuning constant for the bisquare loss functions ('c'), and the
#'        residual degrees of freedom.
#' @details The tuning constant for the bisquare function defaults to
#'        \code{c=3.443689} providing 85\% efficiency for Gaussian data.
#'        Maronna et al suggest bisquare weight functions and 85\% efficiency
#'        with MM-estimation in Sections 5.9 and 11.2 of their book. In the setting of
#'        eliminating baseline effects from clinical trial data, the models
#'        considered are fairly simple and these defaults appear to work well.
#'        The value of 3.443689 is 'borrowed' from \code{lmRob} and its support
#'        functions in the 'robust' package. Rounded values for various Gaussian
#'        efficiencies appear in Seciton 2.2 of Maronna et al.
#'
#'        Following the Second Edition of the book, the function now works with
#'        Peno-Yohai initial estimates with specified defaults. If the Peno-Yohai
#'        procedure fails, least trimmed squares is attempted.
#'
#'        Note that \code{rlm} produces an object
#'        that inherits from \code{lm}, but \code{lmr} does not in order to avoid
#'        \code{lm} methods being used on it. Also, \code{lmr} attaches the
#'        residual degrees of freedom to the object; \code{rlm} deliberately sets
#'        it to NA to avoid erroneous application of \code{lm} methods.
#' @references Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0
#'             Maronna, R. A, Martin, R. D and Yohai, V. J. (2006) Robust Statistics: Theory and Methods, Wiley
#' @keywords models
#' @importFrom MASS rlm
#' @export lmr
lmr <- function(formula, data, weights, psi = "bisquare", method = "MM",
                c = 3.443689, engine = "lmrob", maxit = 50, ...){
  thecall <- match.call()

  inits <- try(py(formula = formula, data = data), silent = TRUE)
  if (class(inits) == "try-error"){
    message("pyinit failed, trying lts")
    inits <- MASS::ltsreg(formula, data)
    if (engine == "rlm"){
      inits <- coef(inits)
    } else if (engine == "lmrob"){
      inits <- list(coefficients = coef(inits), scale = inits$scale[1])
    }
  }

  if (engine == "rlm") {
    if (psi == "bisquare") psiFun <- MASS::psi.bisquare

    res <- MASS::rlm(formula, data, psi = psiFun, method = method, c = c,
                     maxit = maxit, init = inits, ...)

    s <- summary(res)
    res$cov <- s$cov.unscaled * s$sigma^2
    res$data <- data # used by boxplot.rlm
    res$c <- c
    res$call <- thecall

    res$df.residual <- length(res$residuals) - length(res$coefficients)
    res$formula <- formula
    res$psiFun <- res$psi
    res$psi <- psi
    res$cc <- c
    res$scale <- res$s
    res$s <- NULL
    res$init <- inits
  } else if (engine == "lmrob"){
    psi <- if (is.null(psi)) psi <- "bisquare"

    res <- robustbase::lmrob(formula, data, init = inits,
                             control=robustbase::lmrob.control(tuning.psi = c,
                                                               max.it = maxit))
    res$formula <- formula

    res$data <- data

  } else {
    stop("engine should be either 'rlm' or 'lmrob'")
  }

  if (engine == "rlm") {
    class(res) <- c("lmr", "rlm") # drop "lm" because it can lead to errors
  } else {
    class(res) <- c("lmr", "lmrob")
  }
  res
}

#' @method summary lmr
#' @export
summary.lmr <- function(object, ...){
  if ("rlm" %in% class(object)){
    object$s <- object$scale
    object$psi <- object$psiFun
    res <- MASS:::summary.rlm(object, ...)
    co <- res$coefficients
    p <- 2 * pt(abs(co[, 3]), df = object$df.residual, lower.tail = FALSE)
    co <- cbind(co, p)
    colnames(co)[4] <- "Pr(>|t|)"
    res$coefficients <- co
  } else if ("lmrob" %in% class(object)){
    res <- robustbase:::summary.lmrob(object, ...)
  } else {
    stop("the object class should include either 'rlm' or 'lmrob'")
  }

  class(res) <- c("summary.lmr", class(res))
  res
}

#' @method print summary.lmr
#' @export
print.summary.lmr <- function(x, ...){
  if (inherits(x, "summary.lmrob")){
    ## Stop summarizing weights and algorithmic stuff
    x$rweights <- NULL
    robustbase:::print.summary.lmrob(x, showAlgo = FALSE, ...)
  } else if (inherits(x, "summary.rlm")){
    MASS:::print.summary.rlm(x, ...)
  }
  invisible()
}

#' @method predict lmr
#' @export
predict.lmr <- function(object, newdata = NULL, interval = "conf", level = 0.95, ...){
  if (inherits(object, "lmrob")){
    robustbase:::predict.lmrob(object, newdata = newdata, interval = interval, level = level, ...)
  } else if (inherits(object, "rlm")){
    if (is.null(newdata)){
      newdata <- object$data
    }
    object$s <- object$scale
    object$qr <- qr(sqrt(object$weights) * object$x)
    class(object) <- c(class(object), "lm")
    suppressWarnings(stats::predict.lm(object, newdata = newdata, interval = interval, level = level, ...))
  }
}

#' @method update lmr
#' @export
update.lmr  <- function(object, data, ...){
  lmr(object$formula, data = data)
}
