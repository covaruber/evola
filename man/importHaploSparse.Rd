\name{importHaploSparse}
\alias{importHaploSparse}

\title{
Import haplotypes
}
\description{
Formats haplotype in a matrix format to an AlphaSimR population that can be used 
to initialize a simulation. This function serves as wrapper for newMapPop that
utilizes a more user friendly input format.
}
\usage{
importHaploSparse(haplo, genMap, ploidy = 2L, ped = NULL)
}

\arguments{
  \item{haplo}{
a sparse matrix of haplotypes
}
  \item{genMap}{
genetic map as a data.frame. The first three columns must be: marker name, 
chromosome, and map position (Morgans). Marker name and chromosome are coerced
using as.character. See importGenMap.
}
  \item{ploidy}{
ploidy level of the organism.
}
  \item{ped}{
an optional pedigree for the supplied genotypes. See details.
}

}
\details{

The optional pedigree can be a data.frame, matrix or a vector. If the object is
a data.frame or matrix, the first three columns must include information in the
following order: id, mother, and father. All values are coerced using as.character.
If the object is a vector, it is assumed to only include the id. In this case, 
the mother and father will be set to "0" for all individuals.

}
\value{
\describe{
\item{$res}{a MapPop-class if ped is NULL, otherwise a NamedMapPop-class}
}
}

\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

}

\examples{

haplo <- Matrix::Matrix(0, nrow=4, ncol=5)
for (i in 1:4) {
  haplo[i,] <- ifelse(runif(5)<0.2,0,1)
}
colnames(haplo) = letters[1:5]

genMap = data.frame(markerName=letters[1:5],
                    chromosome=c(1,1,1,2,2),
                    position=c(0,0.5,1,0.15,0.4))

ped = data.frame(id=c("a","b"),
                 mother=c(0,0),
                 father=c(0,0))

founderPop = importHaploSparse(haplo=haplo, 
                         genMap=genMap,
                         ploidy=2L,
                         ped=ped)


}
\seealso{\code{\link{evolafit}} -- the core function of the package }
