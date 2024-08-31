\name{stan}
\alias{stan}

\title{
Standardize a vector of values in range 0 to 1
}
\description{
Simple function to map a vector of values to the range of 0 and 1 values to have a better behavior of the algorithm.
}
\usage{
stan(x)
}

\arguments{
  \item{x}{
A vector of numeric values.
}

}
\details{

Simple function to map a vector of values to the range of 0 and 1 values to have a better behavior of the algorithm.

}
\value{
\describe{
\item{$res}{new values in range 0 to 1}
}
}

\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

}

\examples{


x <- rnorm(20, 10, 3);x
stan(x)


}
\seealso{\code{\link{evolafit}} -- the core function of the package }
