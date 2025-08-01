\name{DT_cpdata}
\alias{DT_cpdata}
\alias{A}
\docType{data}
\title{
Genotypic and Phenotypic data for a CP population
}
\description{
A CP population or F1 cross is the designation for a cross between 2 highly heterozygote individuals; i.e. humans, fruit crops, bredding populations in recurrent selection.

This dataset contains phenotpic data for 363 siblings for an F1 cross. These are averages over 2 environments evaluated for 4 traits; color, yield, fruit average weight, and firmness. The columns in the CPgeno file are the markers whereas the rows are the individuals. The CPpheno data frame contains the measurements for the 363 siblings, and as mentioned before are averages over 2 environments.

}
\usage{data("DT_cpdata")}
\format{
  The format is:
 chr "DT_cpdata"
}
\source{
This data was simulated for fruit breeding applications.
}
\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.

}
\examples{

data(DT_cpdata)
DT <- DT_cpdata
A <- A[DT$id,DT$id]

\donttest{
# get best 20 individuals weighting variance by ~0.5=(30*pi)/180
res<-evolafit(formula=cbind(Yield, occ)~id, dt= DT, 
              # constraints: if sum is greater than this ignore 
              constraintsUB = c(Inf,20), 
              # constraints: if sum is smaller than this ignore
              constraintsLB= c(-Inf,-Inf), 
              # weight the traits for the selection
              b = c(1,0), 
              # population parameters
              nCrosses = 100, nProgeny = 10, 
              recombGens=1, nChr=1, mutRateAllele=0,
              # coancestry parameters
              D=A, lambda= (30*pi)/180 , nQtlStart = 20, 
              # selection parameters
              propSelBetween = 0.5, propSelWithin =0.5, 
              nGenerations = 40) 
              
Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
best = bestSol(res$pop)[,"fitness"];best
qa = (Q \%*\% DT$Yield)[best,]; qa 
qDq = Q[best,] \%*\% A \%*\% Q[best,]; qDq 
sum(Q[best,]) # total # of inds selected

pmonitor(res)

plot(DT$Yield, col=as.factor(Q[best,]), 
     pch=(Q[best,]*19)+1)

pareto(res)
 
}

}
\keyword{datasets}
