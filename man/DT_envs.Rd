\name{DT_envs}
\alias{DT_envs}
\alias{W}
\docType{data}
\title{
Simulated phenotypic records for a set of 100 environments and 200 entries.
}
\description{

This is a simulated dataset that aims to show how to pick the most representative
locations and how to optimize sparse testing allocation.

}
\usage{data("DT_envs")}
\format{
  The format is:
 chr "DT_envs"
}
\source{
This data was simulated for experimental design applications.
}
\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.

}
\examples{

data(DT_envs)
# scale passport data for environments
S=scale(W)
# build a relationship matrix between environments
G=tcrossprod(S)/ncol(S)
G[1:4,1:4]
# build the dataset for the genetic algorithm
DT <- data.frame(cov=rnorm(nrow(G)), occ=1, id=rownames(G))
nLocs=10 # desired number of representative locations

\donttest{

res<-evolafit(formula=cbind(cov, occ)~id, dt= DT, 
              # constraints: if sum is greater than this ignore 
              constraintsUB = c(Inf, Inf), 
              # constraints: if sum is smaller than this ignore
              constraintsLB= c(-Inf, nLocs), 
              # weight the traits for the selection
              b = c(1,0), 
              # population parameters
              nCrosses = 100, nProgeny = 20, 
              recombGens=1, nChr=1, mutRate=0,
              # coancestry parameters
              D=G, lambda= 0 , nQtlStart = nLocs, 
              # selection parameters
              fitnessf = inbFun, selectTop = FALSE,
              # fixNumQtlPerInd = TRUE,
              # propSelBetween = 0.5, propSelWithin =0.5, 
              nGenerations = 50, verbose = TRUE) 

Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
best <- bestSol(res$pop, selectTop = FALSE)[,"fitness"]
solution <- Q[best,]
sum(solution) # total # of inds selected
names(solution[which(solution>0)])

pmonitor(res)
pareto(res)

svdW <- svd(G, nu = 2, nv = 2)
PCW <- G \%*\% svdW$v
rownames(PCW) <- rownames(G)
plot(PCW[,1], PCW[,2], col = solution+1, 
    pch=as.character(solution+1),  xlab = "pc1", ylab = "pc2")
 
}

}
\keyword{datasets}
