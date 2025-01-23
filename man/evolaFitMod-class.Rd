\name{evolaFitMod-class}
\docType{class}
\alias{evolaFitMod-class}
\alias{update,evolaFitMod-method}
\title{Genetic algorithm for a set of traits and constraints}
\description{A genetic algorithm fit by \code{\link{evolaFitMod}}.
  This class extends class \code{"\link[AlphaSimR:Pop-class]{Pop}"} class and includes some
  additional slots.}
\section{Objects from the Class}{Objects are created by calls to the
  \code{\link{evolafit}} function.}
\section{Slots}{
  \describe{
  
  \item{indivPerformance}{the matrix of q'a (score), deltaC, q'Dq, generation, nQTNs per solution per generation. See details section above. All other slots are inherited from class \code{"\link[AlphaSimR:Pop-class]{Pop}"}. }
  \item{pedBest}{if the argument keepBest=TRUE this contains the pedigree of the selected solutions across iterations. All other slots are inherited from class \code{"\link[AlphaSimR:Pop-class]{Pop}"}. }
  \item{$score}{a matrix with scores for different metrics across n generations of evolution. All other slots are inherited from class \code{"\link[AlphaSimR:Pop-class]{Pop}"}. }
  \item{$pheno}{the matrix of phenotypes of individuals/solutions present in the last generation. All other slots are inherited from class \code{"\link[AlphaSimR:Pop-class]{Pop}"}. }
  \item{$phenoBest}{the matrix of phenotypes of top (parents) individuals/solutions present in the last generation. All other slots are inherited from class \code{"\link[AlphaSimR:Pop-class]{Pop}"}. }
  \item{constCheckUB}{A matrix with as many rows as solutions and columns as traits to be constrained. 0s indicate that such trait went beyond the bound in that particular solution. All other slots are inherited from class \code{"\link[AlphaSimR:Pop-class]{Pop}"}. }
  \item{constCheckLB}{A matrix with as many rows as solutions and columns as traits to be constrained. 0s indicate that such trait went beyond the bound in that particular solution. All other slots are inherited from class \code{"\link[AlphaSimR:Pop-class]{Pop}"}. }
  \item{traits}{a character vector corresponding to the name of the variables used in the fitness function. All other slots are inherited from class \code{"\link[AlphaSimR:Pop-class]{Pop}"}. }

  }
}
\section{Extends}{
Class \code{"\link[AlphaSimR:Pop-class]{Pop}"}, directly.
}
\section{Methods}{
  \describe{
    \item{update}{\code{signature(object = "evolaFitMod")}: also a
      non-method for the same reason as \code{update}}
  }
}
\seealso{\code{\link{evolaFitMod}}
}
\examples{
showClass("evolaFitMod")
}
\keyword{classes}
