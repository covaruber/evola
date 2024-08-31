\name{Jc}
\alias{Jc}

\title{
Matrix of ones
}
\description{
Makes a matrix of ones with a single row and nc columns. 
}
\usage{
Jc(nc)
}

\arguments{
  \item{nc}{
Number of columns to create.
}

}
\details{

A simple apply function to make a matrix of one row and nc columns.

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

Jc(5)

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
