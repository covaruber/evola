\name{Jr}
\alias{Jr}

\title{
Matrix of ones
}
\description{
Makes a matrix of ones with a single column and nr rows. 
}
\usage{
Jr(nr)
}

\arguments{
  \item{nr}{
Number of rows to create.
}

}
\details{

A simple apply function to make a matrix of one column and nr rows.

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

Jr(5)

}
\seealso{\code{\link{evolafit}} -- the core function of the package }
