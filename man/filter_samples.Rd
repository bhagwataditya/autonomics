% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_filter.R
\name{filter_samples}
\alias{filter_samples}
\title{Filter samples on condition}
\usage{
filter_samples(object, condition, verbose = FALSE, record = TRUE)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{condition}{filter condition}

\item{verbose}{TRUE or FALSE (default)}

\item{record}{TRUE (default) or FALSE}
}
\value{
filtered SummarizedExperiment
}
\description{
Filter samples on condition
}
\examples{
file <- download_data('atkin18.metabolon.xlsx')
object <- read_metabolon(file, plot=FALSE)
filter_samples(object, Group != 't0', verbose = TRUE)
}
