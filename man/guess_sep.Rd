% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_fit.R
\name{guess_sep}
\alias{guess_sep}
\alias{guess_sep.character}
\alias{guess_sep.factor}
\alias{guess_sep.SummarizedExperiment}
\title{Guess separator}
\usage{
guess_sep(x, ...)

\method{guess_sep}{character}(x, separators = c(".", "_"), verbose = FALSE, ...)

\method{guess_sep}{factor}(x, ...)

\method{guess_sep}{SummarizedExperiment}(x, var = "sample_id", separators = c(".", "_"), verbose = FALSE, ...)
}
\arguments{
\item{x}{character vector or SummarizedExperiment}

\item{...}{used for proper S3 method dispatch}

\item{separators}{character vector: possible separators to look for}

\item{verbose}{TRUE or FALSE}

\item{var}{svar or fvar}
}
\value{
separator (string) or NULL (if no separator could be identified)
}
\description{
Guess separator
}
\examples{
# charactervector
   x <- c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]', 'PERM_NON.R3[H/L]')
   guess_sep(x)

   x <- c('WT untreated 1', 'WT untreated 2', 'WT treated 1')
   guess_sep(x)

   x <- c('group1', 'group2', 'group3.R1')
   guess_sep(x)

# SummarizedExperiment
   # file <- download_data('halama18.metabolon.xlsx')
   # object <- read_metabolon(file, plot=FALSE)
   # guess_sep(object)

   # file <- download_data('billing16.proteingroups.txt')
   # object <- read_proteingroups(file, plot=FALSE)
   # guess_sep(object)
}
