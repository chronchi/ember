% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_base_plot.R
\name{get_base_plot}
\alias{get_base_plot}
\title{Plot PCA embedding of the 3 big cohorts to serve as a base plot}
\usage{
get_base_plot(
  x = "PC3",
  y = "PC4",
  color = "pam50",
  size_dots = 2,
  alpha_val = 0.1,
  size_legend = 4,
  base_size = 10
)
}
\arguments{
\item{x}{A string. PC for x-axis}

\item{y}{A string. PC for y-axis}

\item{color}{A string. Which column to use when coloring
the dots.}

\item{size_dots}{An integer. Size of the points in the plot}

\item{alpha_val}{A number. Alpha value of the points}

\item{size_legend}{An integer. Size of the color legend}

\item{base_size}{An integer. Parameter for ggplot2::theme_bw}
}
\value{
A plot of the TCGA, SCANB and METABRIC embedding
}
\description{
Plot PCA embedding of the 3 big cohorts to serve as a base plot
}
