\name{regFun}
\alias{regFun}

\title{
Fitness function from linear regressions based on mean squared error.
}
\description{
Simple function for fitness where the mean squared error is computed when the user provides y and X and b are the average allelic effects of the population in the genetic algorithm. 
}
\usage{
regFun(Y,b,d,Q,a,X,y)
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
  \item{a}{
A named list with vectors of average allelic effects per trait. See details.
}
  \item{X}{
A matrix of covariates or explanatory variables. See details.
}
  \item{y}{
A vector of the response variable. See details.
}

}
\details{

A simple apply function of a regular mean squared error.

( y - X\%*\%b ) ^ 2

Internally, we use this function in the following way:

The y vector and X matrix are provided by the user and are fixed values that do not change.
The evolutionary algorithm optimizes the b values which are the QTLs and associated average
allelic effects that are evolving.

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

y <- rnorm(40) # 4 responses
X=matrix(rnorm(120),40,3) # covariates
Q=matrix(0,40,30) # QTL matrix with 30 QTLs
for(i in 1:nrow(Q)){Q[i,sample(1:ncol(Q),3)]=1}
a <- list(trait1=rnorm(30)) # 30 average allelic effects in trait 1

mse = regFun(y=y, X=X, Q=Q, a=a, # used
             # ignored, Y is normally available in the evolafit routine
             Y=X, b=NA, d=NA) 

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
