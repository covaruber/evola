\name{inbFun}
\alias{inbFun}

\title{
Fitness function from contribution theory using only the group relationship
}
\description{
Simple function for fitness where we only use the group relationship. 
}
\usage{
inbFun(Y,b,d,Q,D,a, lambda)
}

\arguments{
  \item{Y}{
A matrix of trait values. See details.
}
  \item{b}{
A vector of trait weights. See details.
}
  \item{d}{
A vector of group relationships. See details.
}
  \item{Q}{
A QTL matrix. See details.
}
  \item{D}{
An LD matrix. See details.
}
  \item{a}{
A named list with vectors of average allelic effects per trait. See details.
}
  \item{lambda}{
A numeric value to weight the Q'DQ portion of the objective function (to be provided by the user with the lambda argument). See details.
}

}
\details{

A simple apply function of a regular index weighted by a vector of relationships.

Matrix::diag(Q\%*\%Matrix::tcrossprod(D,Q)) of dimensions n x n

Notice that Q represents the marker of QTLs (columns) for all solutions (rows) and D the LD between QTLs. The user can modify this function as needed and provide it to the evolafit function along with other arguments.

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

Q <- matrix(1,3,3) # QTL matrix available internally
D <- diag(3) # LD matrix
inbFun(Q=Q, D=D) # group relationship

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
