\name{drift}
\alias{drift}

\title{
Drift of positive allele calculation
}
\description{
This function takes a population and the simulation parameters from a .Pop class and 
calculates the current frequencies of the positive alleles for each trait.
}
\usage{
drift(pop, simParam)
}

\arguments{

\item{pop}{an object of class \code{"Pop"}}
\item{simParam}{the simulation parameters stored in the evolaMod object.}

}
\details{

A simple apply function to get all the QTLs, segregation sites, positive alleles
and frecuency of the positive allele for each trait in a population.

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

res00<-evolafit(formula=cbind(Weight,Value)~Color, dt= Gems,
               # constraints on traits: if greater than this ignore
               constraintsUB = c(10,Inf), nGenerations = 2
)

drift(pop=res00$pop, simParam=res00$simParam )

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
