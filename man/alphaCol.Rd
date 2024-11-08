% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alphaCol.R
\name{alphaCol}
\alias{alphaCol}
\title{Transparent Colour Creator}
\usage{
alphaCol(col, alpha = 0.5)
}
\arguments{
\item{col}{A valid colour name, as listed by \link{colors} function. OR a positive
integer indexing the colour listed in the \link{palette} function. OR a string of
the form "#rrggbb", see \link{rgb}.}

\item{alpha}{The level of transparency, a \code{numeric} in the range \link{0,1}.}

\item{maxValue}{The maximum colour value, default to 255.}
}
\value{
A character vector of the hexidecimal value defining the colour.
}
\description{
Defining a transparent version of a given colour.
}
\details{
Either a valid colour name or RGB values must be passed to the function. If
both \code{name} and \code{rgbVals} are passed to the function, \code{name} is used to create
the colour.
}
\examples{
colOne = alphaCol(col="blue", alpha=0.5)
colTwo = alphaCol(col=6, alpha=0.1)
colThree = alphaCol(rgb(0.6, 0.1, 0.6), alpha=0.8)


}
