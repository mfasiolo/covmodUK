% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sel_pred_elem.R
\name{sel_elem}
\alias{sel_elem}
\title{Prediction Element Index}
\usage{
sel_elem(d)
}
\arguments{
\item{d}{The dimension of the outcome}
}
\value{
The index selected
}
\description{
Select the prediction element
}
\examples{
idx <- sel_elem(4)
idx(3,1)
}
