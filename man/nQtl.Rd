\name{nQtl}
\alias{nQtl}

\title{
Matrix of number of activated QTLs
}
\description{
Makes a matrix indicating how many QTLs were activated for each solution. 
}
\usage{
nQtl(object)
}

\arguments{
  \item{object}{
Object returned by the evolafit function.
}

}
\details{

A simple apply function to count the number of active QTLs per solution (row) per trait (columns).

}
\value{
\describe{
\item{$res}{a matrix}
}
}

\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

}

\examples{

set.seed(1)

# Data
Gems <- data.frame(
  Color = c("Red", "Blue", "Purple", "Orange",
            "Green", "Pink", "White", "Black", 
            "Yellow"),
  Weight = round(runif(9,0.5,5),2),
  Value = round(abs(rnorm(9,0,5))+0.5,2),
  Times=c(rep(1,8),0)
)
head(Gems)

 \donttest{

# Task: Gem selection. 
# Aim: Get highest combined value.
# Restriction: Max weight of the gem combined = 10. 

# simple specification
res00<-evolafit(formula=cbind(Weight,Value)~Color, dt= Gems,
               # constraints on traits: if greater than this ignore
               constraintsUB = c(10,Inf), nGenerations = 10
)
nQtl(res00)

 }

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
