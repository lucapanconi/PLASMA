% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MPPSimulate.R
\name{mppsimulate}
\alias{mppsimulate}
\title{MPP Simulate.}
\usage{
mppsimulate(
  number_of_frames = 100,
  number_of_points = 500,
  domain_size_x = 1000,
  domain_size_y = NULL,
  search_radius = 50,
  L_target = 100,
  D_max = 50,
  tomato_threshold = NULL,
  grid_size = 100,
  number_of_categories = 1,
  number_of_resimulated_points = NULL,
  use_csr = FALSE,
  probability_of_picking_categories = c(0.7),
  precision_alpha = 0,
  precision_beta = 0,
  category_mean = c(1),
  category_sd = c(0),
  background_mean = 0,
  background_sd = 0,
  createplots = FALSE,
  saveplots = FALSE,
  location = "",
  plotname = "Plot",
  heatmapname = "Heatmap",
  outputall = FALSE
)
}
\arguments{
\item{number_of_frames}{Number of frames to run initial simulation over.}

\item{number_of_points}{Number of points in initial simulation.}

\item{domain_size_x}{Domain size in x direction.}

\item{domain_size_y}{Domain size in y direction. Taken as same as x if omitted.}

\item{search_radius}{Search radius defining neighbourhood.}

\item{L_target}{Target L value.}

\item{D_max}{Diffusivity multiplier for offset.}

\item{tomato_threshold}{Threshold used in ToMATo clustering to determine domains. Will be automatically estimated if omitted.}

\item{grid_size}{Number of bins in each direction of grid which new point cloud is overlaid onto.}

\item{number_of_categories}{Number of different domain types (not including background).}

\item{number_of_resimulated_points}{Number of points in final MPP, will be same as number_of_points if omitted.}

\item{use_csr}{Boolean value dictating whether to use completely spatially random distribution. If TRUE, probabilities are ignored.}

\item{probability_of_picking_categories}{Vector of probabilities of assigning points to domain types.}

\item{precision_alpha}{Alpha parameter used in Gamma distribution to simulate localisation error.}

\item{precision_beta}{Beta parameter used in Gamma distribution to simulate localisation error.}

\item{category_mean}{Vector of mean values for different domain types}

\item{category_sd}{Vector of standard deviations for different domain types}

\item{background_mean}{Mean background value.}

\item{background_sd}{Background standard deviation.}
}
\value{
Offset.
}
\description{
Simulate a marked point pattern given a range of input parameters for greater control of output.
}
