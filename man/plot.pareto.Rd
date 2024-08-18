\name{pareto}
\alias{pareto}
\title{plot the change of values across iterations}
\description{\code{plot} for monitoring.
}
\usage{
pareto(object, scaled=TRUE,pch=20, xlim, \dots)
}
\arguments{
 \item{object}{model object of class \code{"evolafit"}}
 \item{scaled}{a logical value to specify the scale of the y-axis (gain in merit).}
 \item{pch}{symbol for plotting points as desribed in par}
 \item{xlim}{upper and lower bound in the x-axis}
 \item{\dots}{Further arguments to be passed to the plot function.}
}
\value{
  vector of plot
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{plot}}, \code{\link{evolafit}}}
\keyword{models}

