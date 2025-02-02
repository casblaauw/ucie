% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data_transformations.R
\name{ScalingGuess}
\alias{ScalingGuess}
\title{Estimate an initial scaling parameter}
\usage{
ScalingGuess(dataset, polygon)
}
\arguments{
\item{dataset}{Data to fit inside the polygon.
Only uses the minimum and maximum of each of 3 dimensions, so inputting convex hull
vertices of the data (like from [DataConvex()] returns the same result.}

\item{polygon}{Vertices of the polygon to fit inside (e.g. CIELAB convex hull vertices)}
}
\value{
A named integer S to use as a first scaling value.
}
\description{
Takes the range (max - min) in all 3 dimensions for both data and polygon,
does simple 'alignment' by matching smallest to smallest and biggest to biggest
dimensions, then returns the smallest ratio between ranges: i.e. the scaling
factor necessary for the first dimension to 'bump' into the polygon limits.
}
\keyword{internal}
