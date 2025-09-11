\name{ocsFunC}
\alias{ocsFunC}

\title{
Fitness function from contribution theory
}
\description{
Simple function for fitness where an index of traits is weighted by the group relationship. 
}
\usage{
ocsFunC(Y,b,Q,D,a,lambda,scaled=TRUE,
        SNP, solution, alphaLog=1, 
        weightsTraitFreq=NULL,...)
}

\arguments{
  \item{Y}{
A matrix of trait values. See details.
}
  \item{b}{
A vector of trait weights. See details.
}
  \item{Q}{
A QTL matrix for solutions. See details.
}
  \item{omega}{
a numerical value between 0 and 1 to specify the trade-off between frequency-based contribution and
breeding value contribution when using the \code{\link{ocsFunC}} fitness function.
}
  \item{scaled}{
A logical value to indicate if traits should be scaled prior to multiply by the weights.
}
  \item{SNP}{
A SNP marker matrix for the QTLs if the QTLs are individuals.
}
  \item{solution}{
an object of RRsol-class.
}
 \item{alphaLog}{
an integer value to be applied to function \code{\link{logspace}} when using the \code{\link{ocsFunC}}
fitness function. Such fitness function is part of the evola package.
}
  \item{weightsTraitFreq}{
A vector of trait weights for the solution object. See details.
}
  \item{...}{
additional arguments to pass.
}

}
\details{

A function of the type:

Y\%*\%b + Y2\%*\%b2

Internally, we use this function in the following way:

The Y matrix is the matrix of trait-GEBVs (Y = Q\%*\%alpha) for each solution, and b is the user-specified trait weights. 

The Y2 matrix is the matrix of positive-allele-trait-frequencies for each solution, and b2 is the user-specified trait-solution weights. 

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


data(DT_technow)
DT <- DT_technow
DT$occ <- 1; DT$occ[1]=0
M <- M_technow+1

A <- A.mat(M)
A <- A[DT$hy,DT$hy]

mixm = rrBLUP::mixed.solve(y=DT$GY, Z=M)
solution = matrix(mixm$u[colnames(M)])

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
              fitnessf = ocsFunC, omega=1,
              SNP=M, solution=solution,
              nGenerations = 20) 

  Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
  best = bestSol(res$pop)[,"fitness"]
  qa = (Q \%*\% DT$GY)[best,]; qa 
  qAq = Q[best,] \%*\% A \%*\% Q[best,]; qAq 
  sum(Q[best,]) # total # of inds selected

pmonitor(res)
plot(DT$GY, col=as.factor(Q[best,]), 
     pch=(Q[best,]*19)+1)

pareto(res)

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
