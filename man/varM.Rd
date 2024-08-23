\name{varM}
\alias{varM}

\title{
Extract the variance existing in the genome solutions
}
\description{
Extracts the variance found across the M element of the resulting object of the evolafit() function which contains the different solution and somehow represents the genome of the population. 
}
\usage{
varM(object)
}

\arguments{
  \item{object}{
A resulting object from the function evolafit.
}

}
\details{

A simple apply function looking at the variance in each column of the M element of the resulting object of the evolafit function.

}
\value{
\describe{
\item{$res}{a value of variance}
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
res0<-evolafit(cbind(Weight,Value)~Color, dt= Gems,
               # constraints: if greater than this ignore
               constraintsUB = c(10,Inf), 
               # constraints: if smaller than this ignore
               constraintsLB= c(-Inf,-Inf), 
               # weight the traits for the selection
               traitWeight = c(0,1), 
               # population parameters
               nCrosses = 100, nProgeny = 20, recombGens = 1, 
               # coancestry parameters
               A=NULL, lambda=c(0,0), nQTLperInd = 1, 
               # selection parameters
               propSelBetween = .9, propSelWithin =0.9, 
               nGenerations = 50
) 

varM(res0)
}

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
