#' Robust analysis of variance type tests.
#' @param object,... Robust linear models fit by \code{robustbase::lmrob}.
#' @details The first object must be the model with the least terms.The results of
#'   the robust F-test and robust chi-squared test are both returned. The second
#'   of these should, generally, be preferred. The function is a simple
#'   wrapper to \link{\code{lmrobLinTest}}.
#' @method anova lmr
#' @export
anova.lmr <- function(object, ...){
  o <- rev(list(object, ...))
  res <- list()

  if (inherits(object, "lmrobdetMM")){
    o <- lapply(o, function(X){
      X$control$psi <- X$control$family
      X
    })
  } else if (inherits(object, "rlm")){
    o <- lapply(o, function(X){
      X$control <- list(psi = X$psi, tuning.psi = X$cc)
      X
    })
  }

  for (i in 2:length(o)){
    res[[i - 1]] <- lmrobLinTest(o[[i - 1]], o[[i]])
  }
  res[[i]] <- list(test = NA, chisq.pvalue = NA, F.pvalue = NA,
                   df = c(NA, length(resid(object))))

  res <- as.data.frame(t(sapply(res, unlist))[, c(4, 5, 1, 3, 2)])

  res[, 1] <- as.integer(res[, 1])
  res[, 2] <- as.integer(res[, 2])
  res <- res[nrow(res):1, ]

  fos <- lapply(o, formula)

  fos <- sapply(fos, function(X){
    paste(X[2], X[3], sep = " ~ ")
  })

  rownames(res) <- rev(fos)
  colnames(res)[4:5] <- c("P(>F)", "P(>Chi)")

  res
}


#' Robust likelihood ratio test for linear hypotheses
#'
#' This function computes a robust likelihood ratio test for linear hypotheses.
#'
#' @export  lmrobLinTest
#' @rdname lmrobLinTest
#'
#' @param object1 an \code{lmrob} object with the fit corresponding to the complete model
#' @param object2 an \code{lmrob} object with the fit corresponding to the model
#' restricted under the null linear hypothesis.
#'
#' @return A list with the following components: c("test","chisq.pvalue","f.pvalue","df")
#' \item{test}{The value of the F-statistic}
#' \item{f.pvalue}{p-value based on the F distribution}
#' \item{chisq.pvalue}{p-value based on the chi-squared distribution}
#' \item{df}{degrees of freedom}
#'
#' @author Victor Yohai, \email{vyohai@gmail.com}
#' @references \url{http://www.wiley.com/go/maronna/robust}
lmrobLinTest <- function (object1, object2){
  rho <- RobStatTM::rho
  rhoprime <- RobStatTM::rhoprime
  rhoprime2 <- RobStatTM::rhoprime2

  p <- length(object1$coeff)
  q <- length(object1$coeff) - length(object2$coeff)
  n <- length(na.omit(object1$resid))

  if (q <= 0){
    stop("Contrary to other anova functions, the largest model must be the first.")
  }

  if (length(na.omit(object2$resid)) != n){
    stop("Models were not fit to the same data (missing values?)")
  }

  family1 <- object1$control$psi
  family2 <- object2$control$psi

  cc1 <- object1$control$tuning.psi
  cc2 <- object2$control$tuning.psi

  same.name <- (family1 == family2)
  same.cc <- isTRUE(all.equal(cc1, cc2))

  if (!(same.name && same.cc)){
    stop("object1 and object2 do not have the same rho family")
  }

  s <- object1$scale

  a <- sum(rho(object2$resid/s, family = family1, cc = cc1,
               standardize = TRUE))
  b <- sum(rho(object1$resid/s, family = family1, cc = cc1,
               standardize = TRUE))
  c <- sum(rhoprime2(object1$resid/s, family = family1, cc = cc1,
                     standardize = TRUE))
  d <- sum(rhoprime(object1$resid/s, family = family1, cc = cc1,
                    standardize = TRUE)^2)

  test <- (2 * (a - b) * c)/d
  chisq.pvalue <- 1 - pchisq(test, q)
  F.pvalue <- 1 - pf(test/q, q, n - p)
  df <- c(q, n - p)
  output <- list(test = test, chisq.pvalue = chisq.pvalue,
                 F.pvalue = F.pvalue, df = df)
  output
}
