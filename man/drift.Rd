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

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)

#Make 10 crosses
pop2 = randCross(pop, 10, simParam=SP)

drift(pop=pop, simParam=SP )

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
