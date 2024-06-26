% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allfunctions.R
\name{computeIOUfromjoseph}
\alias{computeIOUfromjoseph}
\title{Compute IOU from JOSEPH.}
\usage{
computeIOUfromjoseph(vertices_1, josephresults)
}
\arguments{
\item{vertices_1}{Data frame or matrix. Vertices of first set of hulls. Represented as a matrix with three columns: x, y and vertex label.}

\item{josephresults}{Data frame or matrix. Output from JOSEPH clustering.}
}
\value{
IOU score.
}
\description{
Computes the IOU between a defined set of hulls and the hulls identified from JOSEPH clustering.
}
