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
drift(pop, simParam, solution=NULL, traits=1)
}

\arguments{

\item{pop}{an object of class \code{"Pop"}}
\item{simParam}{the simulation parameters stored in the evolaMod object.}
\item{solution}{an alternative RRBLUP solution so the frequencies of positive alleles are calculated for a given SNP chip. If NULL the frequencies are calculated for the real QTLs. }
\item{traits}{traits considered.}

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
founderPop = quickHaplo(nInd=10, nChr=1, segSites=30)

#Set simulation parameters
SP = SimParam$new(founderPop)

SP$addTraitA(10)
SP$setVarE(h2=0.5)
SP$addSnpChip(15)

#Create population
pop = newPop(founderPop, simParam=SP)

# drift for QTLs
drift(pop, simParam = SP ) 

\donttest{

ans = RRBLUP(pop, simParam=SP)

# drift for average allelic effects
drift(pop, simParam = SP, solution = ans)

}

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
