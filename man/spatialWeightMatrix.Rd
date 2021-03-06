% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{spatialWeightMatrix}
\alias{spatialWeightMatrix}
\title{spatialWeightMatrix}
\usage{
spatialWeightMatrix(x, span = NULL)
}
\arguments{
\item{x}{a matrix with rows corresponding to cells and columns
corresponding to dimensions to calculate Euclidean distance}

\item{span}{proportion of samples to include on either side,
default is 13/(number of rows in `x`), corresponding roughly
to points within a diamond shape distance away}
}
\value{
A weighted matrix
}
\description{
Create weight matrix for spatial data
}
\examples{

spat_x <- rnorm(50)
spat_y <- rnorm(50)
spat_coord <- cbind(spat_x, spat_y)
W <- spatialWeightMatrix(spat_coord)

}
