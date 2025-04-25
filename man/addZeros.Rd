\name{addZeros}
\alias{addZeros}

\title{
Function to add zeros before and after a numeric vector to have the same number
of characters.
}
\description{
Function to add zeros before and after a numeric vector to have the same number
of characters. 
}
\usage{
addZeros(x, nr=2)
}

\arguments{
  \item{x}{
Numeric vector.
}

  \item{nr}{
number of digits to keep to the right.
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

addZeros(5)

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
