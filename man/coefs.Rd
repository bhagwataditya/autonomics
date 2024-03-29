% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_fit.R
\name{coefs}
\alias{coefs}
\alias{coefs.factor}
\alias{coefs.data.table}
\title{Get coefs}
\usage{
coefs(featuredt, ...)

\method{coefs}{factor}(featuredt, ...)

\method{coefs}{data.table}(featuredt, fit = fits(featuredt), svars = NULL, ...)
}
\arguments{
\item{featuredt}{data.table or factor}

\item{...}{required for s3 dispatch}

\item{fit}{string: 'limma', 'lm', 'lme', 'lmer'}

\item{svars}{NULL or charactervector (svar for which to return coefs)}
}
\value{
character vector
}
\description{
Get coefs
}
\examples{
# Factor
    x <- factor(c('A', 'B', 'C'))
    coefs(x)
    coefs(code(x, contr.treatment.explicit))
    coefs(code(x, code_control))
    
# SummarizedExperiment
    file <- download_data('atkin.metabolon.xlsx')
    object <- read_metabolon(file, fit = 'limma')
    coefs(fdt(object))
}
