% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GPCalculate.R
\name{informdeviance}
\alias{informdeviance}
\title{Inform Deviance}
\usage{
informdeviance(
  values,
  k = 2,
  createplot = FALSE,
  saveplot = FALSE,
  location = "",
  outputall = FALSE,
  plotresolution = 100
)
}
\arguments{
\item{values}{Marked values.}
}
\value{
Deviance threshold.
}
\description{
Try to determine an appropriate deviance threshold for JOSEPH based on the distribution of marked values.
}
