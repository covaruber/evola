\name{fitnessf}
\alias{fitnessf}

\title{
Fitness function for contribution theory
}
\description{
Simple function for fitness where an index of traits is weighted by the group relationship. 
}
\usage{
fitnessf(Y,b,d,Q)
}

\arguments{
  \item{Y}{
A matrix of traits. See details.
}
  \item{b}{
A vector of weights. See details.
}
  \item{d}{
A vector of group relationships. See details.
}
  \item{Q}{
A QTL matrix. See details.
}

}
\details{

A simple apply function of a regular index weighted by a vector of relationships.

Y\%*\%b - d

Internally, we use this function in the following way:

The Y matrix is the matrix of trait-GEBVs and b is the user-specified trait weights. 

d = xtAx * lambda; where xtAx is equal to Matrix::diag(Q\%*\%Matrix::tcrossprod(A,Q))

Notice that Q represents the marker of QTLs (columns) for all solutions (rows). The user can modify this function as needed and provide it to the evolafit function along with other arguments.

}
\value{
\describe{
\item{$res}{a vector of values}
}
}

\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

}

\examples{


Y=matrix(1:9,3,3)
b=c(1,2,3)
d=c(1,1,1)
fitnessf(Y,b,d)

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
