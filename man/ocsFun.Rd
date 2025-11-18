\name{ocsFun}
\alias{ocsFun}

\title{
Fitness function from contribution theory
}
\description{
Simple function for fitness where an index of traits is weighted by the group relationship. 
}
\usage{
ocsFun(Y,b,Q,D,lambda,scaled=TRUE,...)
}

\arguments{
  \item{Y}{
A matrix of trait values. See details.
}
  \item{b}{
A vector of trait weights. See details.
}
  \item{Q}{
A QTL matrix for solutions. See details.
}
  \item{D}{
An LD matrix for solutions. See details.
}
  \item{lambda}{
A numeric value to weight the Q'DQ portion of the objective function (to be provided by the user with the lambda argument). See details.
}
  \item{scaled}{
A logical value to indicate if traits should be scaled prior to multiply by the weights.
}
  \item{...}{
additional arguments to pass.
}

}
\details{

A simple apply function of a regular index weighted by a vector of relationships.

Y\%*\%b - d

Internally, we use this function in the following way:

The Y matrix is the matrix of trait-GEBVs (Y = Q\%*\%alpha) for each solution, and b is the user-specified trait weights. 

d = qtDq * lambda; where qtDq is equal to Matrix::diag(Q\%*\%Matrix::tcrossprod(D,Q)) of dimensions n x n

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


Y=matrix(1:12,4,3) # 4 solutions with 3 traits

Q=matrix(sample(0:1,32,replace = TRUE),nrow=4,ncol=8) # coancestry for each solution

D=diag(8)

b=rep(1,3)

lambda=0.5

ocsFun(Y=Y,Q=Q,D=D,b=b, lambda=lambda) # Yb - d where d is QAQ' and A is the LD between QTNs


}
\seealso{\code{\link{evolafit}} -- the core function of the package }
