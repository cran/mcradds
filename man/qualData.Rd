% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{qualData}
\alias{qualData}
\title{Simulated Qualitative Data}
\format{
A \link{qualData} data set contains 200 observations and 3 variables.
\describe{
\item{Sample}{Sample id}
\item{ComparativeN}{Measurements from comparative analyzer with \code{1=positive} and \code{0=negative}}
\item{CandidateN}{Measurements from candidate analyzer with \code{1=positive} and \code{0=negative}}
}
}
\usage{
qualData
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

This simulated data \code{\link{qualData}} can be used to calculate the
qualitative performance such as sensitivity and specificity.
}
\seealso{
\link{platelet} that contains quantitative data comparing
platelet results from two analyzers.
}
\keyword{datasets}
