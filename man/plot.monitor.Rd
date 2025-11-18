\name{evolmonitor}
\alias{evolmonitor}
\title{plot the change of values across iterations}
\description{\code{plot} for monitoring.
}
\usage{
evolmonitor(object, kind, \dots)
}
\arguments{
 \item{object}{model object of class \code{"evolafit"}}
 \item{kind}{a numeric value indicating what to plot according to the following values:
 
 1: Average and best q'a (contribution)
 
 2. Average q'Dq and deltaC 
 
 3. Number of QTLs activated
 
 }
 \item{\dots}{Further arguments to be passed to the plot function.}
}
\value{
  trace plot
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{plot}}, \code{\link{evolafit}}}
\keyword{models}

