\name{freqPosAllele}
\alias{freqPosAllele}

\title{
Extract frequency of positive alleles
}
\description{
Computes the frequency of positive alleles given certain marker effects. 
}
\usage{
freqPosAllele(M, alpha)
}

\arguments{
  \item{M}{
Marker matrix.
}
  \item{alpha}{
Vector of marker effects.
}

}
\details{

A simple apply function to compute the frequency of the positive alleles.

}
\value{
\describe{
\item{$res}{a matrix}
}
}

\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

}

\examples{

#random population of 10 lines with 5 markers
M <- matrix(rep(0,10*5),10,5)
for (i in 1:10) {
  M[i,] <- ifelse(runif(5)<0.5,0,2)
}
alpha <- sample(c(-1,1),5, replace=TRUE)

freqPosAllele(M, alpha)

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
