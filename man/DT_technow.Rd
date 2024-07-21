\name{DT_technow}
\alias{DT_technow}
\alias{Ad_technow}
\alias{Af_technow}
\alias{Md_technow}
\alias{Mf_technow}
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
Md <- (Md_technow*2) - 1
Mf <- (Mf_technow*2) - 1

\dontrun{


## first get all possible hybrids
res1 <- sommer::build.HMM(Md, Mf, 
                          return.combos.only = TRUE)
head(res1$data.used)
## build the marker matrix for the first 50 hybrids
pred <- which(res1$data.used$hybrid %in% DT$hy)
res2 <- sommer::build.HMM(Md, Mf,
                          custom.hyb = res1$data.used[pred,]
)
res2$HMM.add[1:5,1:5]
res2$HMM.dom[1:5,1:5]
## now you can use the A.mat()
M <- res2$HMM.add
A <- sommer::A.mat(M)
DT$GY <- imputev(DT$GY)
DT$occ <- 1; DT$occ[1]=0
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
              nCrosses = 50, nProgeny = 30, 
              # coancestry parameters
              A=A, lambda=c(0.4,0), nQTLperInd = 100, 
              # selection parameters
              propSelBetween = 0.9, propSelWithin =0.9, 
              nGenerations = 20) 

best = which(res$pheno[,1]==max(res$pheno[,1]))[1];best
xa = (res$M %*% DT$GY)[best,]; xa 
xAx = res$M[best,] %*% A %*% res$M[best,]; xAx 
sum(res$M[best,]) # total # of inds selected

plot(NULL, xlim =c(1,nrow(res$score)), 
     ylim = c(0, max(res$score) ), 
     ylab = "Value", xlab = "Generation")
legend("bottomright", 
       c("Best x'a", "Average x'a", "Average x'Ax", "Mean nQTL"), 
       pch = 20, col = c(1, 2, 3, 4),bty = "n")
lines(res$score[,2])
lines(res$score[,1], col = "red")
lines(res$score[,3], col = "green")
lines(res$score[,4], col = "blue")
plot(DT$GY, col=as.factor(res$M[best,]), 
     pch=(res$M[best,]*19)+1)


}

}
\keyword{datasets}
