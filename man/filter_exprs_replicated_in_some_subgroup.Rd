% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_filter.R
\name{filter_exprs_replicated_in_some_subgroup}
\alias{filter_exprs_replicated_in_some_subgroup}
\title{Filter features with replicated expression in some subgroup}
\usage{
filter_exprs_replicated_in_some_subgroup(
  object,
  subgroupvar = "subgroup",
  comparator = if (contains_ratios(object)) "!=" else ">",
  lod = 0,
  verbose = TRUE
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{subgroupvar}{subgroup svar}

\item{comparator}{'>' or '!='}

\item{lod}{number: limit of detection}

\item{verbose}{TRUE or FALSE}
}
\value{
Filtered SummarizedExperiment
}
\description{
Filter features with replicated expression in some subgroup
}
\examples{
require(magrittr)
file <- download_data('atkin18.metabolon.xlsx')
object <- read_metabolon(file, plot=FALSE)
object \%<>\% filter_exprs_replicated_in_some_subgroup(subgroupvar = 'Group')
filter_exprs_replicated_in_some_subgroup(object, character(0))
filter_exprs_replicated_in_some_subgroup(object, NULL)
}
