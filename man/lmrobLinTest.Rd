% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aov.lmr.R
\name{lmrobLinTest}
\alias{lmrobLinTest}
\alias{anova.lmr}
\title{Robust likelihood ratio test for linear hypotheses}
\usage{
lmrobLinTest(object1, object2)
}
\arguments{
\item{object1}{an \code{lmrob} object with the fit corresponding to the complete model}

\item{object2}{an \code{lmrob} object with the fit corresponding to the model
restricted under the null linear hypothesis.}
}
\value{
A list with the following components: c("test","chisq.pvalue","f.pvalue","df")
\item{test}{The value of the F-statistic}
\item{f.pvalue}{p-value based on the F distribution}
\item{chisq.pvalue}{p-value based on the chi-squared distribution}
\item{df}{degrees of freedom}
}
\description{
Robust likelihood ratio test for linear hypotheses
}
\details{
This function computes a robust likelihood ratio test for linear hypotheses.
  The first object must be the model with the least terms.The results of
  the robust F-test and robust chi-squared test are both returned. The second
  of these should, generally, be preferred. The function \code{anova.lmr} is a simple
  wrapper to \code{lmrobLinTest} in package RobStatTM. \code{anova.lmr}
  will take multiple model objects, not just 2, and it rearranges the output.
}
\references{
\url{http://www.wiley.com/go/maronna/robust}
}
\author{
Victor Yohai, \email{vyohai@gmail.com}
}
