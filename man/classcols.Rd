% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classcols.R
\name{classcols}
\alias{classcols}
\title{Classification Colour}
\usage{
classcols(classif, alphaLevel = 1)
}
\arguments{
\item{classif}{A numeric vector for the deficiency probabilities. All elements must be between 0.0 and 1.0.}

\item{alphaLevel}{A number between 0.0 and 1.0 to control the level of the opacity of the colours.}
}
\value{
A character vector hex value RGB names.
}
\description{
The colour given to a fibre given its probability of being deficient. The
colours are on a scale of red to blue where red corresponds to deficiency i.e.
a probability of 1.0.
}
\details{
The colour gradient is created using \code{\link[=ramp_BlueRed]{ramp_BlueRed()}}.
}
