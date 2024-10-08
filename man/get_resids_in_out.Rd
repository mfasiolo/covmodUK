% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_resids_in_out.R
\name{get_resids_in_out}
\alias{get_resids_in_out}
\title{Get uniform residuals from GAMLSS models on test and training sets}
\usage{
get_resids_in_out(mod, final, logp)
}
\arguments{
\item{mod}{The output of the fitted GAMLSS models}

\item{final}{Same as above but for the models fitted to the whole data}

\item{logp}{if TRUE the uniform residuals are on log scale}
}
\value{
A list of in-sample and out-of-sample residuals
}
\description{
Get uniform residuals from GAMLSS models on test and training sets
}
