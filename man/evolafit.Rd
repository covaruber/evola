\name{evolafit}
\alias{evolafit}

\title{
Fits a genetic algorithm for a set of traits and constraints.
}
\description{
Using the AlphaSimR machinery it recreates the evolutionary forces applied to a 
problem where possible solutions replace individuals and combinations of variables 
to optimize in the problem replace the genes or QTLs. Then evolutionary forces 
(mutation, selection and drift) are applied to find a 
close-to-optimal solution. Although multiple traits are enabled it is assumed 
that same QTLs are behind all the traits, differing only in their average allelic effects.
}
\usage{

evolafit(formula, dt, 
     constraintsUB, constraintsLB, constraintW=NULL, 
     b, nCrosses=50, nProgeny=20,nGenerations=20, 
     recombGens=1, nChr=1, mutRateAllele=0,  mutRateAlpha=0,
     nQtlStart=NULL, D=NULL, lambda=0,
     propSelBetween=NULL,propSelWithin=NULL,
     fitnessf=NULL, verbose=TRUE, dateWarning=TRUE,
     selectTop=TRUE, tolVarG=1e-6, 
     Ne=50, initPop=NULL, simParam=NULL,
     fixNumQtlPerInd=FALSE, traceDelta=TRUE, topN=10,
     includeSet=NULL, excludeSet=NULL,
     ...)

}

\arguments{

  \item{formula}{
Formula of the form \emph{y~x where} \code{y} refers to the average allelic 
substitution effects of the QTLs (alpha) for each trait, and \code{x} refers to 
the variable defining the genes or QTLs to be combined in the possible solutions.
}
  \item{dt}{
A dataset containing the average allelic effects (a) and classifiers/genes/QTLs.
}
  \item{constraintsUB}{
A numeric vector specifying the upper bound constraints for the breeding values 
applied at each trait. The length is equal to the number of traits/features in 
the formula. If missing is assume an infinite value for all traits. Solutions 
(individuals in the population) with higher value than the upper bound are 
assigned a -infinite value if the argument \code{selectTop=TRUE} and to 
+infinite when \code{selectTop=FALSE}, which is equivalent to reject/discard a 
solution based on the fitness function.
}
  \item{constraintsLB}{
A numeric vector specifying the lower bound constraints for the breeding values 
applied at each trait. The length is equal to the number of traits/features 
in the formula. If missing is assume a -infinite value for all traits. Solutions 
with lower value than the lower bound are assigned a +infinite value if the 
argument \code{selectTop=TRUE} and to -infinite when \code{selectTop=FALSE}, 
which is equivalent to reject/discard a solution based on the fitness function.
}
  \item{constraintW}{
A numeric vector of length equal to the number of generations to specify the weights 
to be applied to the lower and upper bound constraint values to relax the constraints
if needed. If values is NULL a vector of 1s with length equal to the number of
generations is created.
}
  \item{b}{
A numeric vector specifying the weights that each trait has in the fitness function 
(i.e., a selection index). The length should be equal to the number of traits/features. 
If missing is assumed equal weight (1) for all traits.
}
  \item{nCrosses}{
A numeric value indicating how many crosses should occur in the population of 
solutions at every generation.
}
  \item{nProgeny}{
A numeric value indicating how many progeny (solutions) each cross should generate 
in the population of solutions at every generation.
}
  \item{nGenerations}{
The number of generations that the evolutionary process should run for.
}
  \item{recombGens}{
The number of recombination generations that should occur before selection is
applied. This is in case the user wants to allow for more recombination before 
selection operates. The default is 1.
}
  \item{nChr}{
The number of chromosomes where features/genes should be allocated to. The default 
value is 1 but this number can be increased to mimic more recombination events at 
every generation and avoid linkage disequilibrium. 
}
  \item{mutRateAllele}{
A value between 0 and 1 to indicate the proportion of random QTL dosage that should 
mutate in each individual. For example, a value of 0.1 means that a random 10\% 
of the QTLs will mutate in each individual randomly taking values of 0 or 1. 
Is important to notice that this implies that the stopping criteria based in 
variance will never be reached because we keep introducing variance through 
random mutation. 
}
  \item{mutRateAlpha}{
A value between 0 and 1 to indicate the proportion of random QTL average allelic effects that should 
mutate in each individual. For example, a value of 0.1 means that a random 10\% 
of the QTLs will mutate in each individual randomly taking values of 0 or 1. 
Is important to notice that this implies that the stopping criteria based in 
variance will never be reached because we keep introducing variance through 
random mutation. 
}
  \item{nQtlStart}{
The number of QTLs/genes (classifier \code{x} in the formula) that should be 
fixed for the positive allele at the begginning of the simulation. If not specified
it will be equal to the 20\% of the QTLs (calculated as the number of rows in the 
dt argument over 5). This is just an initial value and will change as the population
evolve under the constraints specified by the user. See details section.
}
  \item{D}{
A relationship matrix between the QTLs (a kind of linkage disequilibrium) specified 
in the right side of the formula (levels of the \code{x} variable). This matrix 
can be used or ignored in the fitness function. By default the weight to the q'Dq 
component is 0 though the lambda argument, where x is an individual in the population 
of a solution.
}
  \item{lambda}{
A numeric value indicating the weight assigned to the relationship between QTLs 
in the fitness function. If not specified is assumed to be 0. This can be used 
or ignored in your customized fitness function.
}
  \item{propSelBetween}{
A numeric value between 0 and 1 indicating the proportion of families/crosses of 
solutions/individuals that should be selected. The default is 1, meaning all 
crosses are selected or passed to the next generation.
}
  \item{propSelWithin}{
A numeric value between 0 and 1 indicating the proportion of individuals/solutions 
within families/crosses that should be selected. The default value is 0.5, meaning 
that 50\% of the top individuals are selected.
}
  \item{fitnessf}{
An alternative fitness function to be applied at the level of individuals or 
solutions. It could be a linear combination of the trait breeding values. 
The available variables internally are:

\strong{Y}: matrix of trait breeding values for the individuals/solutions. Of 
dimensions s x t, s soultions and t traits.

\strong{b}: vector of trait weights, specified in the 'b' argument. Of dimensions 
t x 1, t traits by 1

\strong{Q}: matrix with QTLs for the individuals/solutions. Of dimensions s x p,
s solutions and p QTL columns. Although multiple traits are enabled it is assumed 
that same QTLs are behind all the traits, differing only in their average allelic 
effects.

\strong{D}: matrix of relationship between the QTLs, specified in the 'D' argument. 
Of dimensions p x p, for p QTL columns

\strong{lambda}: a numeric value indicating the weight assigned to the relationship 
between QTLs in the fitness function. If not specified is assumed to be 0. This can 
be used or ignored in your customized fitness function.

\strong{a}: list of vectors with average allelic effects for a given trait. Of 
dimensions s x 1, s solutions by 1 column

If \code{fitnessf=NULL}, the default function will be the \code{\link{ocsFun}} function:

\code{ function(Y,b,d,Q,D,a,lambda)\{(Y\%*\%b) - d\} }

where \code{(Y\%*\%b)} is equivalent to [(Q'a)b] in genetic contribution theory, 
and \code{d} is equal to the diagonal values from Q'DQ from contribution theory, 

\strong{\emph{If you provide your own fitness function please keep in mind that 
the variables Y, b, Q, D, a, and lambda are already reserved and these variables 
should always be added to your function (even if you do not use them) in addition 
to your new variables so the machinery runs.}}

An additional fitness function for accounting only for the group relationship is
\code{\link{inbFun}} when the user wants to find solutions that maximize the 
representativeness of a sample and the \code{D} argument is not NULL. 
You will need to select the solutions with lower values ( \code{selectTop=TRUE} ) 
which indicate solutions with more representativeness and you may need to indicate
lower bound constraints ( \code{constraintsLB} ).

An additional fitness function available for regression problems is \code{\link{regFun}}
but is not the default since it would require additional arguments not available 
in a regular genetic algorithm problem (e.g., \code{y} and \code{X} to compute 
\code{y-Xb} ).

}
  \item{verbose}{
A logical value indicating if we should print logs.
}
  \item{dateWarning}{
A logical value indicating if you should be warned when there is a new version
on \code{CRAN}.
}
  \item{selectTop}{
Selects highest values for the fitness value if \code{TRUE}. Selects lowest 
values if \code{FALSE}.
}
  \item{tolVarG}{
A stopping criteria (tolerance for genetic variance) when the variance across 
traits is smaller than this value, which is equivalent to assume that all solutions 
having the same QTL profile (depleted variance). The default value is \code{1e-6} 
and is computed as the sum of the diagonal values of the genetic variance covariance 
matrix between traits.
}
  \item{Ne}{
initial number of founders in the population (will be important for long term 
sustainability of genetic variance).
}
  \item{initPop}{
an object of \code{Pop-class}.
}
  \item{simParam}{
  an object of \code{SimParam}.
}
  \item{fixNumQtlPerInd}{
A \code{TRUE/FALSE} value to indicate if we should fix the argument nQtlStart 
across all generations. This should be used with care since this is not how 
usually genetic algorithms work and in my experience only using GA for regression 
problems is a special case where this argument should be set to TRUE. The behavior 
assumes that if set to TRUE and a particular solution has more QTLs active than 
nQtlStart some QTLs will be set to 0 and if a solution has less QTLs active than 
nQtlStart some QTLs will be activated. All activations or deactivations are done 
at random. This could be an alternative to use a counting trait to restrain the number
of QTLs active in a solution but is slower.
}
  \item{traceDelta}{
  a logical value indicating if we should compute the rate of coancestry Q'DQ at each iteration. 
  This metric is used by the pareto plot but is not needed for the evolutionary process
  and it can take a considerable amount of time when the number of QTLs is big.
}
  \item{topN}{
  an integer value indicating the maximum number of solutions to keep in each generation.
}
  \item{includeSet}{
  a numeric vector with 0s and 1s of length equal to the number of QTLs (number of rows in your 
  dataset) to indicate which individuals (1s) should be forced to be activated.
}
  \item{excludeSet}{
   a numeric vector with 0s and 1s of length equal to the number of QTLs (number of rows in your 
  dataset) to indicate which individuals (1s) should be forced to be deactivated.
}
\item{\dots}{Further arguments to be passed to the fitness function if required.}

}
\details{

Using the \code{AlphaSimR} machinery (runMacs) it recreates the evolutionary 
forces applied to a problem where possible solutions replace individuals and 
combinations of variables in the problem replace the genes. Then evolutionary 
forces are applied to find a close-to-optimal solution. The number of solutions
are controlled with the nCrosses and nProgeny parameters, whereas the number of 
initial QTLs activated in a solution is controlled by the nQtlStart parameter. 
The number of activated QTLs of course will increase if has a positive effect in 
the fitness of the solutions. The drift force can be controlled by the recombGens
parameter. The mutation rate can be controlled with the  mutRateAllele parameter. The 
recombination rate can be controlled with the nChr argument.

The \code{indivPerformance} output slot contains the columns id, fitness, 
generation, nQTLs, and deltaC. These mean the following:

\strong{In fitness} : represents the fitness function value of a solution.

\strong{In deltaC} : it represents the change in coancestry (e.g., inbreeding), 
it can be thought as the rate of coancestry. It is calculated as q'Dq where 'q' 
represents the contribution vector, 'D' is the linkage disequilibrium matrix 
between QTNs (whatever the QTNs represent for your specific problem). In practice 
we do QAQ' and extract the diagonal values.

\strong{In generation} : it represents the generation at which this solution appeared.

\strong{In nQTNs} : it represent the final number of QTNs that are activated in 
homozygote state for the positive effect.

During the run the columns printed in the console mean the following:

\strong{generation}: generation of reproduction

\strong{constrainedUB}: number of solutions constrained by the upper bound specified

\strong{constrainedLB}: number of solutions constrained by the lower bound specified

\strong{varG}: genetic variance present in the population due to the QTNs

\strong{propB}: proportion of families selected during that iteration

\strong{propW}: proportion of individuals within a family selected in that iteration

\strong{time}: the time when the iteration has finished.


}
\value{

\describe{

\item{indivPerformance}{the matrix of fitness, deltaC, generation, nQTNs per 
                        solution per generation. See details section above. }
\item{pedBest}{contains the pedigree of the selected solutions across iterations.}
\item{$score}{a matrix with scores for different metrics across n generations of 
              evolution.}
\item{$pheno}{the matrix of phenotypes of individuals/solutions present in the 
              last generation.}
\item{pop}{AlphaSimR object used for the evolutionary algorithm at the last iteration.}
\item{constCheckUB}{A matrix with as many rows as solutions and columns as traits 
                    to be constrained. 0s indicate that such trait went beyond 
                    the bound in that particular solution.}
\item{constCheckLB}{A matrix with as many rows as solutions and columns as traits
                    to be constrained. 0s indicate that such trait went beyond 
                    the bound in that particular solution.}
\item{traits}{a character vector corresponding to the name of the variables used 
              in the fitness function.}

}

}

\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for 
complex problems. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R 
package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. 
https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence
Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.

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

# simple specification
res00<-evolafit(formula=cbind(Weight,Value)~Color, dt= Gems,
               # constraints on traits: if greater than this ignore
               constraintsUB = c(10,Inf), nGenerations = 10
)
best = bestSol(res00$pop)[,"fitness"]
Q <- pullQtlGeno(res00$pop, simParam = res00$simParam, trait=1); Q <- Q/2
qa = Q[best,] \%*\% as.matrix(Gems[,c("Weight","Value")]); qa


# more complete specification
res0<-evolafit(formula=cbind(Weight,Value)~Color, dt= Gems,
               # constraints on traits: if greater than this ignore
               constraintsUB = c(10,Inf), 
               # constraints on traits: if smaller than this ignore
               constraintsLB= c(-Inf,-Inf), 
               # weight the traits for the selection (fitness function)
               b = c(0,1), 
               # population parameters
               nCrosses = 100, nProgeny = 20, 
               # genome parameters
               recombGens = 1, nChr=1,  mutRateAllele=0, nQtlStart = 2, 
               # coancestry parameters
               D=NULL, lambda=0, 
               # selection parameters
               propSelBetween = .9, propSelWithin =0.9, 
               nGenerations = 50
) 

Q <- pullQtlGeno(res0$pop, simParam = res0$simParam, trait=2); Q <- Q/2
best = bestSol(res0$pop)[,"fitness"]
qa = Q[best,] \%*\% as.matrix(Gems[,c("Weight","Value")]); qa
Q[best,]

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
