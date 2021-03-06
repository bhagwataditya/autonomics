% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6_transform.R
\name{explore_transformations}
\alias{explore_transformations}
\title{Explore transformations}
\usage{
explore_transformations(
  object,
  subgroup = subgroup,
  transformations = c("quantnorm", "zscore", "invnorm"),
  method = "pca",
  xdim = 1,
  ydim = 2,
  ...
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{subgroup}{subgroup (sym)}

\item{transformations}{vector}

\item{method}{'pca', 'pls', 'sma', or 'lda'}

\item{xdim}{number (default 1)}

\item{ydim}{number (default 2)}

\item{...}{passed to plot_data}
}
\value{
grid object
}
\description{
Explore transformations
}
\examples{
file <- download_data('billing16.proteingroups.txt')
invert <- c('EM_E', 'EM_BM', 'BM_E')
object <- read_proteingroups(file, invert_subgroups = invert, plot=FALSE)
explore_transformations(object)
}
