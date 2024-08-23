\name{evolafit}
\alias{evolafit}

\title{
Fits a genetic algorithm for a set of traits and constraints.
}
\description{
Using the AlphaSimR machinery it recreates the evolutionary forces applied to a problem where possible solutions replace individuals and combinations of variables in the problem replace the genes. Then evolutionary forces are applied to find a close-to-optimal solution. 
}
\usage{

evolafit(formula, dt, 
     constraintsUB, constraintsLB, traitWeight,
     nCrosses=50, nProgeny=40,nGenerations=30, recombGens=1,
     nQTLperInd=NULL, A=NULL, lambda=NULL,
     propSelBetween=1,propSelWithin=0.5,
     fitnessf=NULL, verbose=TRUE, dateWarning=TRUE)

}

\arguments{

  \item{formula}{
Formula of the form y~x where 'y' refers to the traits or features involved in selecting or putting constraints to the solutions, and 'x'refers to the classifiers of the problem to define the chromosome regions.
}
  \item{dt}{
A dataset containing the features and classifiers.
}
  \item{constraintsUB}{
A numeric vector specifying the upper bound constraints in the traits/features (y). The length is equal to the number of traits/features. If missing is assume infinite for all traits.
}
  \item{constraintsLB}{
A numeric vector specifying the lower bound constraints in the traits/features (y). The length is equal to the number of traits/features. If missing is assume minus infinite for all traits.
}
  \item{traitWeight}{
A numeric vector specifying the weights that each trait has in the final selection index. The length is equal to the number of traits/features. If missing is assumed equal weight for all traits.
}
  \item{nCrosses}{
A numeric value indicating how many crosses should occur in the population of solutions at every generation.
}
  \item{nProgeny}{
A numeric value indicating how many progeny each cross should generate in the population of solutions at every generation.
}
  \item{nGenerations}{
The number of generations that the evolutionary process should run for.
}
  \item{recombGens}{
The number of recombination generations that should occur before selection is applied. This is in case the user wants to allow for more recombination before selection operates.
}
  \item{nQTLperInd}{
The number of levels corresponding to the classifier (x) that should show up as present at the begginning of the simulation. If not specified it will be equal to the number of rows in the dataset/5. See details section.
}
  \item{A}{
A relationship matrix between the levels of the classifier variable (x).
}
  \item{lambda}{
A numeric value indicating the weight assigned to the relationship between levels of the classifiers with respect to the selection index. If not specified is assumed to be 0.
}
  \item{propSelBetween}{
A numeric value between 0 and 1 indicating the proportion of families/crosses that should be selected.
}
  \item{propSelWithin}{
A numeric value between 0 and 1 indicating the proportion of individuals within families/crosses that should be selected.
}
  \item{fitnessf}{
An alternative named list defining the fitness function for a given trait. If NULL the default function will be xa - (lam * x'Ax) where x is the contribution vector to the solution, a are the original values assigned to the trait-merit, and A is the covariance between solutions.
}
  \item{verbose}{
A logical value indicating if we should print logs.
}
  \item{dateWarning}{
A logical value indicating if you should be warned when there is a new version on CRAN.
}

}
\details{

Using the AlphaSimR machinery (runMacs) it recreates the evolutionary forces applied to a problem where possible solutions replace individuals and combinations of variables in the problem replace the genes. Then evolutionary forces are applied to find a close-to-optimal solution. The number of solutions are controlled with the nCrosses and nProgeny parameters, whereas the number of initial combinations present for the classifier/genes is controlled by the nQTLperInd parameter. This of course will increase if has an effect in the fitness. The drift force can be controlled by the recombGens parameter.

}
\value{

\describe{

\item{$M}{the matrix of haplotypes/solutions after selection.}
\item{$score}{a matrix with scores for different metrics () across generations of evolution.}
\item{$pheno}{the matrix of phenotypes of individuals/solutions present in the last generation.}
\item{indivPerformance}{the matrix of x'a, x'Ax, deltaC, nQTLs per solution per generation.}
\item{pop}{AlphaSimR object used for the evolutionary algorithm.}

}

}

\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.

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
head(Gems)
#     Color Weight Value
# 1    Red   4.88  9.95
# 2   Blue   1.43  2.73
# 3 Purple   1.52  2.60
# 4 Orange   3.11  0.61
# 5  Green   2.49  0.77
# 6   Pink   3.53  1.99
# 7  White   0.62  9.64
# 8  Black   2.59  1.14
# 9 Yellow   1.77 10.21

 \donttest{

# Task: Gem selection. 
# Aim: Get highest combined value.
# Restriction: Max weight of the gem combined = 10. 
res0<-evolafit(cbind(Weight,Value)~Color, dt= Gems,
               # constraints: if greater than this ignore
               constraintsUB = c(10,Inf), 
               # constraints: if smaller than this ignore
               constraintsLB= c(-Inf,-Inf), 
               # weight the traits for the selection
               traitWeight = c(0,1), 
               # population parameters
               nCrosses = 100, nProgeny = 20, recombGens = 1, 
               # coancestry parameters
               A=NULL, lambda=c(0,0), nQTLperInd = 1, 
               # selection parameters
               propSelBetween = .9, propSelWithin =0.9, 
               nGenerations = 50
) 

best = bestSol(res0)[2]
xa = res0$M[best,] \%*\% as.matrix(Gems[,c("Weight","Value")]); xa

res0$M[best,]
res0$score[nrow(res0$score),]

# $`Genes`
# Red   Blue Purple Orange  Green   Pink  White  Black Yellow 
# 1      1      0      0      1      0      0      1      0 
# 
# $Result
# Weight  Value 
# 8.74  32.10 
pmonitor(res0)
pareto(res0)

 }

}
\seealso{\code{\link{evolafit}} -- the information of the package }
