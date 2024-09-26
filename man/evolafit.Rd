\name{evolafit}
\alias{evolafit}

\title{
Fits a genetic algorithm for a set of traits and constraints.
}
\description{
Using the AlphaSimR machinery it recreates the evolutionary forces applied to a problem where possible solutions replace individuals and combinations of variables to optimize in the problem replace the genes or QTLs. Then evolutionary forces (mutation, selection and drift) are applied to find a close-to-optimal solution. 
}
\usage{

evolafit(formula, dt, 
     constraintsUB, constraintsLB, traitWeight,
     nCrosses=50, nProgeny=40,nGenerations=30, 
     recombGens=1, nChr=1, mutRate=0,
     nQTLperInd=NULL, A=NULL, lambda=NULL,
     propSelBetween=1,propSelWithin=0.5,
     fitnessf=NULL, verbose=TRUE, dateWarning=TRUE,
     selectTop=TRUE, tolVarG=1e-6, keepBest=FALSE,
     ...)

}

\arguments{

  \item{formula}{
Formula of the form y~x where 'y' refers to the traits or features defining the average allelic substitution effects of the QTLs, and 'x'refers to the variable defining the genes or QTLs to be combined in the possible solutions.
}
  \item{dt}{
A dataset containing the features/traits and classifiers/genes/QTLs.
}
  \item{constraintsUB}{
A numeric vector specifying the upper bound constraints in the traits/features (y). The length is equal to the number of traits/features in the formula. If missing is assume an infinite value for all traits. Solutions with higher value than the upper bound are assigned a -infinite value if the argument selectTop=TRUE and to +infinite when selectTop=FALSE, which is equivalent to reject/discard a solution based on the fitness function.
}
  \item{constraintsLB}{
A numeric vector specifying the lower bound constraints in the traits/features (y). The length is equal to the number of traits/features in the formula. If missing is assume a -infinite value for all traits. Solutions with lower value than the lower bound are assigned a +infinite value if the argument selectTop=TRUE and to -infinite when selectTop=FALSE, which is equivalent to reject/discard a solution based on the fitness function.
}
  \item{traitWeight}{
A numeric vector specifying the weights that each trait has in the fitness function (e.g., a selection index). The length is equal to the number of traits/features. If missing is assumed equal weight (1) for all traits.
}
  \item{nCrosses}{
A numeric value indicating how many crosses should occur in the population of solutions at every generation.
}
  \item{nProgeny}{
A numeric value indicating how many progeny (solutions) each cross should generate in the population of solutions at every generation.
}
  \item{nGenerations}{
The number of generations that the evolutionary process should run for.
}
  \item{recombGens}{
The number of recombination generations that should occur before selection is applied. This is in case the user wants to allow for more recombination before selection operates. The default is 1.
}
  \item{nChr}{
The number of chromosomes where features/genes should be allocated to. The default value is 1 but this number can be increased to mimic more recombination events at every generation. STILL NOT FUNCTIONAL.
}
  \item{mutRate}{
A value between 0 and 1 to indicate the proportion of random QTLs that should mutate in each individual. For example, a value of 0.1 means that a random 10\% of the QTLs will mutate in each individual randomly taking values of 0 or 1. Is important to notice that this implies that the stopping criteria based in variance will never be reached because we keep introducing variance through random mutation. 
}
  \item{nQTLperInd}{
The number of QTLs/genes (classifier x) that should be fixed for the positive allele at the begginning of the simulation. If not specified it will be equal to the 20\% of the QTLs (number of rows in the dt over 5). See details section.
}
  \item{A}{
A relationship matrix between the levels of the classifier variable (x or QTLs; not between the solutions). It is a kind of a linkage disequilibrium matrix. This function can be used or ignored in the fitness function.
}
  \item{lambda}{
A numeric value indicating the weight assigned to the relationship between levels of the classifiers in comparison with the trait value. If not specified is assumed to be 0. This can be used or ignored in your customized fitness function.
}
  \item{propSelBetween}{
A numeric value between 0 and 1 indicating the proportion of families/crosses that should be selected. The default is 1, meaning all crosses are selected.
}
  \item{propSelWithin}{
A numeric value between 0 and 1 indicating the proportion of individuals within families/crosses that should be selected. The default value is 0.5, meaning that 50\% of the top individuals are selected.
}
  \item{fitnessf}{
An alternative fitness function for a linear combination of the traits. If NULL the default function will be:

function(Y,b,d,Q)\{(Y\%*\%b) - d\}

where Y\%*\%b is equivalent to xa in contribution theory, and d is equal to x'Ax, being x the contribution vector to the solution, a are the QTL effects, and A is the covariance between QTLs, Q is the QTL matrix for the solution. If you provide your own fitness function please keep in mind that the variables Y, b, d and Q are already reserved and should always be added to your function (even if not used) in addition to your new variables.

}
  \item{verbose}{
A logical value indicating if we should print logs.
}
  \item{dateWarning}{
A logical value indicating if you should be warned when there is a new version on CRAN.
}
  \item{selectTop}{
Selects highest values for the fitness value if TRUE. Selects lowest values if FALSE.
}
  \item{tolVarG}{
A stopping criteria (tolerance for genetic variance) when the variance across traits is smaller than this value, which is equivalent to assume that all solutions having the same QTL profile (depleted variance). The default value is 1e-6 and is computed as the sum of the diagonal values of the genetic variance covariance matrix between traits.
}
  \item{keepBest}{
A TRUE/FALSE value to indicate if we should store the QTL matrix and pedigree of the solutions selected across generations. This can be useful if we want to recreate the path to the best solution (e.g., best crossing path to a product).
}
\item{\dots}{Further arguments to be passed to the fitness function.}

}
\details{

Using the AlphaSimR machinery (runMacs) it recreates the evolutionary forces applied to a problem where possible solutions replace individuals and combinations of variables in the problem replace the genes. Then evolutionary forces are applied to find a close-to-optimal solution. The number of solutions are controlled with the nCrosses and nProgeny parameters, whereas the number of initial QTLs activated in a solution is controlled by the nQTLperInd parameter. The number of activated QTLs of course will increase if has a positive effect in the fitness of the solutions. The drift force can be controlled by the recombGens parameter. The mutation rate can be controlled with the mutRate parameter. The recombination rate can be controlled with the nChr argument.

}
\value{

\describe{

\item{$M}{the matrix of haplotypes/solutions at the end of the run.}
\item{$Mb}{the matrix of top (parents) haplotypes/solutions at the end of the run.}
\item{$score}{a matrix with scores for different metrics across n generations of evolution.}
\item{$pheno}{the matrix of phenotypes of individuals/solutions present in the last generation.}
\item{$phenoBest}{the matrix of phenotypes of top (parents) individuals/solutions present in the last generation.}
\item{indivPerformance}{the matrix of x'a, x'Ax, deltaC, nQTLs per solution per generation.}
\item{pop}{AlphaSimR object used for the evolutionary algorithm at the last iteration.}
\item{best}{AlphaSimR object corresponding to the best parental haplotypes/solutions selected for new crosses across all cycles.}
\item{pedBest}{if the argument keepBest=TRUE this contains the pedigree of the selected solutions across iterations.}
\item{traits}{a character vector corresponding to the name of the variables used in the fitness function.}

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
res0<-evolafit(formula=cbind(Weight,Value)~Color, dt= Gems,
               # constraints: if greater than this ignore
               constraintsUB = c(10,Inf), 
               # constraints: if smaller than this ignore
               constraintsLB= c(-Inf,-Inf), 
               # weight the traits for the selection
               traitWeight = c(0,1), 
               # population parameters
               nCrosses = 100, nProgeny = 20, 
               # genome parameters
               recombGens = 1, nChr=1, mutRate=0, nQTLperInd = 1, 
               # coancestry parameters
               A=NULL, lambda=0, 
               # selection parameters
               propSelBetween = .9, propSelWithin =0.9, 
               nGenerations = 50
) 

best = bestSol(res0)["pop","Value"]
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
