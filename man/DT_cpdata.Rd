\name{DT_cpdata}
\alias{DT_cpdata}
\alias{GT_cpdata}
\alias{MP_cpdata}
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
GT <- GT_cpdata
MP <- MP_cpdata

\dontrun{
## library sommer is needed for the relationship matrix
library(sommer)
## look at the data and fit the model
DT$Yield <- scale(sommer::imputev(DT$Yield))
DT$occ <- 1; DT$occ[1]=0
head(DT)
## create the variance-covariance matrix
A <- sommer::A.mat(GT) # additive relationship matrix
A <- A + diag(1e-4, ncol(A), ncol(A))
# get best 20 individuals weighting variance by 0.5
res<-evolafit(cbind(Yield, occ)~id, dt= DT, 
              # constraints: if sum is greater than this ignore 
              constraintsUB = c(Inf,20), 
              # constraints: if sum is smaller than this ignore
              constraintsLB= c(-Inf,-Inf), 
              # weight the traits for the selection
              traitWeight = c(1,0), 
              # population parameters
              nCrosses = 50, nProgeny = 40, 
              # coancestry parameters
              A=A, lambda=c(0.5,0), nQTLperInd = 20, 
              # selection parameters
              propSelBetween = 0.9, propSelWithin =0.9, 
              nGenerations = 30) 

best = which(res$pheno[,1]==max(res$pheno[,1]))[1];best
xa = (res$M \%*\% DT$Yield)[best,]; xa 
xAx = res$M[best,] \%*\% A \%*\% res$M[best,]; xAx 
sum(res$M[best,]) # total # of inds selected

plot(NULL, xlim =c(1,nrow(res$score)), ylim = c(0, max(res$score)), 
     ylab = "Value", xlab = "Generation")
legend("bottomright", 
       c("Best x'a", "Average x'a", "Average x'Ax", "Mean nQTL"), 
       pch = 20, col = c(1, 2, 3, 4),bty = "n")
lines(res$score[,2])
lines(res$score[,1], col = "red")
lines(res$score[,3], col = "green")
lines(res$score[,4], col = "blue")

plot(DT$Yield, col=as.factor(res$M[best,]), 
     pch=(res$M[best,]*19)+1)
 
}

}
\keyword{datasets}
