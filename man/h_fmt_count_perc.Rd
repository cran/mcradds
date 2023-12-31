% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{h_fmt_count_perc}
\alias{h_fmt_count_perc}
\title{Format count and percent}
\usage{
h_fmt_count_perc(cnt, perc = NULL, format, ...)
}
\arguments{
\item{cnt}{(\code{numeric})\cr numeric vector for count.}

\item{perc}{(\code{numeric})\cr numeric vector for percent, if Null only format count.}

\item{format}{(\code{string})\cr formatting string from \code{formatters::list_valid_format_labels()}
for \code{formatters::format_value()} function.}

\item{...}{other arguments to be passed to \link[formatters:format_value]{formatters::format_value}.}
}
\value{
A character vector of formatted counts and percents.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

Help function to format the count and percent into one string.
}
\examples{
h_fmt_count_perc(cnt = c(5, 9, 12, 110, 0), format = "xx")
h_fmt_count_perc(
  cnt = c(5, 9, 12, 110, 0),
  perc = c(0.0368, 0.0662, 0.0882, 0.8088, 0),
  format = "xx (xx.x\%)"
)
}
