\name{bestSol}
\alias{bestSol}

\title{
Extract the index of the best solution
}
\description{
Extracts the index of the best solution for all traits under the constraints specified. 
}
\usage{
bestSol(object, selectTop=TRUE)
}

\arguments{
  \item{object}{
A resulting object from the function evolafit.
}
  \item{selectTop}{
Selects highest values for the fitness value if TRUE. Selects lowest values if FALSE.
}

}
\details{

A simple apply function looking at the fitness value of all the solution in the last generation to find the maximum value.

}
\value{
\describe{
\item{$res}{the vector of best solutions in M for each trait in the problem}
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

bestSol(res0)
}

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
