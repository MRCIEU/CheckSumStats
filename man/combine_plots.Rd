% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_plots.R
\name{combine_plots}
\alias{combine_plots}
\title{Make cow plot}
\usage{
combine_plots(
  Plot_list = NULL,
  out_file = NULL,
  return_plot = FALSE,
  width = 800,
  height = 1000,
  Title = "",
  Xlab = "",
  Ylab = "",
  Title_size = 0,
  Title_axis_size = 0,
  by2cols = TRUE,
  Ncol = 2,
  Tiff = FALSE
)
}
\arguments{
\item{Plot_list}{plots to combine. Can either be vector of character strings giving the names of plot objects or a list of plot objects.}

\item{out_file}{filepath to save the plot}

\item{return_plot}{logical argument. If TRUE, plot is returned and is not save to out_file}

\item{width}{width of plot}

\item{height}{height of plot}

\item{Title}{plot title}

\item{Xlab}{label for X axis}

\item{Ylab}{label for Y axis}

\item{Title_size}{size of title}

\item{Title_axis_size}{size of x axis title}

\item{by2cols}{logical argument. If true, forces plot to have 2 columns}

\item{Ncol}{number of columns}

\item{Tiff}{save plot in tiff format. Default is set to FALSE. If set to FALSE, the  plot is saved in png format. Not applicable if return_plot is set to TRUE.}
}
\value{
plot
}
\description{
Combine all plots into a single plot using the cowplot package
}
