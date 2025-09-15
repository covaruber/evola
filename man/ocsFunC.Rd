\name{ocsFunC}
\alias{ocsFunC}

\title{
Fitness function from contribution theory
}
\description{
Simple function for fitness where an index of traits is weighted by the group relationship. 
}
\usage{
ocsFunC(Q,
        SNP, solution,
        weightsTraitFreq=NULL,
        useAlphaWeight=FALSE,
        ...)
}

\arguments{

  \item{Q}{
A QTL matrix for solutions. See details.
}
  \item{SNP}{
A SNP marker matrix for the QTLs if the QTLs are individuals.
}
  \item{solution}{
an object of RRsol-class.
}
  \item{weightsTraitFreq}{
A vector of trait weights for the solution object. See details.
}
  \item{useAlphaWeight}{
A logical value to indicate if the weights for the beneficial allele frequencies 
should be the average allelic effects (TRUE) or the proposed 1-benef.freq.base.pop
as proposed originaly by Covarrubias-Pazaran and Pita (2026).
}
  \item{...}{
additional arguments to pass.
}

}
\details{

A function of the type:

sum( freq.solution * (1 - freq.base.pop) )

the frequency of beneficial alleles is computed for each solution in Q and then
it is weighted by 1 - freq.base.pop which is the frequency of beneficial alleles
in the base population.

}
\value{
\describe{
\item{$res}{a vector of values}
}
}

\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

}

\examples{

# mickey mouse example

Q=matrix(sample(0:1,32,replace = TRUE),nrow=4,ncol=8) # possible solutions

SNP=matrix(sample(0:2,80,replace = TRUE),nrow=8,ncol=10) # SNPs of 8 individuals with 10 markers

solution = matrix(rnorm(10)) # marker effects

ocsFunC(SNP=SNP,Q=Q, solution=solution) # sum of frequency of positive alleles

 \donttest{

  # using the function in an optimal contribution

  data(DT_technow)
  DT <- DT_technow
  DT$occ <- 1; DT$occ[1]=0
  M <- M_technow+1

  A <- A.mat(M)
  A <- A[DT$hy,DT$hy]

  mixm = rrBLUP::mixed.solve(y=DT$GY, Z=M)
  solution = matrix(mixm$u[colnames(M)])

  sum(freqPosAllele(M,solution[,1])) # starting frequencies

  # run the genetic algorithm
  # we assig a weight to x'Dx of (20*pi)/180=0.34
  res<-evolafit(formula = c(GY, occ)~hy,
              dt= DT, 
              # constraints: if sum is greater than this ignore
              constraintsUB = c(Inf,100), 
              # constraints: if sum is smaller than this ignore
              constraintsLB= c(-Inf,-Inf),
              # weight the traits for the selection
              b = c(1,0), 
              # population parameters
              nCrosses = 100, nProgeny = 10, 
              recombGens=1, nChr=1, mutRateAllele=0,
              # coancestry parameters
              D=A, lambda= (20*pi)/180 , nQtlStart = 90, 
              # selection parameters
              propSelBetween = 0.5, propSelWithin =0.5, 
              fitnessf = ocsFunC, 
              SNP=M, solution=solution,
              nGenerations = 20) 

  Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
  best = bestSol(res$pop)[,"fitness"]
  qa = (Q \%*\% DT$GY)[best,]; qa 
  qAq = Q[best,] \%*\% A \%*\% Q[best,]; qAq 
  sum(Q[best,]) # total # of inds selected
  
  sum(freqPosAllele(diag(Q[best,])\%*\%M,solution[,1]))

  pmonitor(res)
  plot(DT$GY, col=as.factor(Q[best,]), 
     pch=(Q[best,]*19)+1)

  pareto(res)

}

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
