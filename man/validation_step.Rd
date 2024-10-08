% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation_step.R
\name{validation_step}
\alias{validation_step}
\title{Validation for to choose optimal number of effects}
\usage{
validation_step(data_boost, dat, start, nstop, ncv, ncores)
}
\arguments{
\item{data_boost}{output of boosting_step() function}

\item{dat}{data on which to do the validation}

\item{start}{Start of validation}

\item{nstop}{Sequence of number of effects on which to evaluate the performance}

\item{ncv}{number of cross-validation folds}

\item{ncores}{number of parallel threads}
}
\value{
A list containing information on the validation scores as a function of the number of effects
}
\description{
Validation for to choose optimal number of effects
}
