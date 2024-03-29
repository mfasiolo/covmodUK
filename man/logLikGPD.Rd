% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logLikGPD.R
\name{logLikGPD}
\alias{logLikGPD}
\title{Log-likelihood of Generalized Pareto Distribution}
\usage{
logLikGPD(y = NULL)
}
\arguments{
\item{np}{XXX.}
}
\description{
XXX.
}
\examples{
library(gamFactory)
n <- 1000
pars <- c(1, 0.3)
obj <- logLikGPD( )$initialize(n, pars)

# Derivatives should match exactly
fdDeriv(obj = derFunWrapper(obj$derObj), 
        param = pars, 
        ord = 1:3)

# Should look fine
der <- derivCheck(np = 100, 
                  parSim = function(n){ cbind(runif(n, 0.1, 1e6), runif(n, -0.45, 0.45)) }, 
                  obj = obj,
                  ord = 1:3, 
                  trans = function(.x){
                    si <- sign(.x)
                    return( si * sqrt(abs(.x)) )
                  }, 
                  n = 1000)

# We see some discrepancies here
par(mfrow = c(2, 2))
for(ii in 1:3) { plot(der[[ii]][ , 1] - der[[ii]][ , 2]) }

# But relative error seems acceptable 
par(mfrow = c(2, 2))
for(ii in 1:3) { plot( (der[[ii]][ , 1] - der[[ii]][ , 2]) / 
                         abs( der[[ii]][ , 2] )) }

}
