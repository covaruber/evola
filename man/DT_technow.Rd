\name{DT_technow}
\alias{DT_technow}
\alias{M_technow}
\docType{data}
\title{
Genotypic and Phenotypic data from single cross hybrids (Technow et al.,2014)
}
\description{
This dataset contains phenotpic data for 2 traits measured in 1254 single cross hybrids coming from the cross of Flint x Dent heterotic groups. In addition contains the genotipic data (35,478 markers) for each of the 123 Dent lines and 86 Flint lines. The purpose of this data is to demosntrate the prediction of unrealized crosses (9324 unrealized crosses, 1254 evaluated, total 10578 single crosses). We have added the additive relationship matrix (A) but can be easily obtained using the A.mat function on the marker data. Please if using this data for your own research cite Technow et al. (2014) publication (see References).

}
\usage{data("DT_technow")}
\format{
  The format is:
 chr "DT_technow"
}
\source{
This data was extracted from Technow et al. (2014).
}
\references{
If using this data for your own research please cite:

Technow et al. 2014. Genome properties and prospects of genomic predictions of hybrid performance in a Breeding program of maize. Genetics 197:1343-1355.

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.

}
\examples{

data(DT_technow)
DT <- DT_technow
DT$occ <- 1; DT$occ[1]=0
M <- M_technow

\donttest{
  A <- A.mat(M)
  # run the genetic algorithm
  res<-evolafit(formula = c(GY, occ)~hy,
                dt= DT, 
                # constraints: if sum is greater than this ignore
                constraintsUB = c(Inf,100), 
                # constraints: if sum is smaller than this ignore
                constraintsLB= c(-Inf,-Inf),
                # weight the traits for the selection
                traitWeight = c(1,0), 
                # population parameters
                nCrosses = 100, nProgeny = 10, 
                # coancestry parameters
                A=A, lambda=c(0.4,0), nQTLperInd = 70, 
                # selection parameters
                propSelBetween = 0.5, propSelWithin =0.5, 
                nGenerations = 20) 
  
  best = bestSol(res)[1]
  xa = (res$M \%*\% DT$GY)[best,]; xa 
  xAx = res$M[best,] \%*\% A \%*\% res$M[best,]; xAx 
  sum(res$M[best,]) # total # of inds selected
  
  pmonitor(res)
  plot(DT$GY, col=as.factor(res$M[best,]), 
       pch=(res$M[best,]*19)+1)
       
  pareto(res)


}

}
\keyword{datasets}
