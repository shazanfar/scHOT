% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/higherOrderFunctions.R
\name{weightedZIKendall}
\alias{weightedZIKendall}
\title{weightedZIKendall}
\usage{
weightedZIKendall(x, y, w = 1)
}
\arguments{
\item{x}{x and y are non-negative data vectors}

\item{y}{x and y are non-negative data vectors}

\item{w}{weight vector, values should be between 0 and 1}
}
\value{
\code{numeric} weighted Tau* association value between x and y
}
\description{
the weightedZIKendall function calculates weighted Tau*,
where Tau* is described in Pimentel et al (2015)
doi:10.1016/j.spl.2014.09.002. This association measure
is defined for zero-inflated, non-negative random variables.
}
\examples{
x = pmax(0,rnorm(100))
y = pmax(0,rnorm(100))
w = runif(100)
weightedZIKendall(x,y,w)

}
