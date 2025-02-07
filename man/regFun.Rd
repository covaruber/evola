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
A matrix of trait values (internal matrix). See details.
}
  \item{b}{
A vector of trait weights (provided by users). See details.
}
  \item{d}{
A vector of group relationships (internal matrix). See details.
}
  \item{Q}{
A QTL matrix (internal matrix). See details.
}
  \item{a}{
A named list with vectors of average allelic effects per trait (internal matrix). See details.
}
  \item{X}{
A matrix of covariates or explanatory variables (to be provided by the user in the ... arguments). See details.
}
  \item{y}{
A vector of the response variable (to be provided by the user in the ... arguments). See details.
}

}
\details{

A simple apply function of a regular mean squared error.

\code{( y - X\%*\%b ) ^ 2}

Internally, we use this function in the following way:

The y vector and X matrix are provided by the user and are fixed values that do not change across iterations.
The evolutionary algorithm optimizes the b values which are the QTLs and associated average
allelic effects that are evolving. The 'b' coefficients in the formula come from the GA and are computed as:

\code{b[j] = a[[1]][p[[j]]]}

where a[[1]] is the list of QTL average allelic effects per trait provided in the original dataset, whereas p is a list (with length
equal to the number of solutions) that indicates which QTLs are activated in each solution and the j variable
is just a counter so each solution is tested.

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
