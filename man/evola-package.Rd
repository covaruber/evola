\name{evola-package}
\alias{evola}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{

\strong{EVOL}utionary \strong{A}lgorithm
\cr

}
\description{

The evola package is nice wrapper of the AlphaSimR package that enables the use of the evolutionary algorithm to solve complex questions in a simple manner.

The \code{\link{evolafit}} function is the core function of the package which allows the user to specify the problem and constraints to find a close-to-optimal solution using the evolutionary forces.

}

\section{Keeping evola updated}{

The evola package is updated on CRAN every 4-months due to CRAN policies but you can find the latest source at https://github.com/covaruber/evola. This can be easily installed typing the following in the R console:

library(devtools)

install_github("covaruber/evola")

This is recommended if you reported a bug, was fixed and was immediately pushed to GitHub but not in CRAN until the next update.
}

\section{Tutorials}{

For tutorials on how to perform different analysis with evola please look at the vignettes by typing in the terminal:

\strong{vignette("evola.intro")}

}

\section{Getting started}{
The package has been equiped with several datasets to learn how to use the evola package: 

* \code{\link{DT_technow}} datasets to perform optimal cross selection. 

* \code{\link{DT_wheat}} TBD.

* \code{\link{DT_cpdata}} dataset to perform optimal individual. 

}

\section{Models Enabled}{

The machinery behind the scenes is AlphaSimR.

}

\section{Bug report and contact}{

If you have any questions or suggestions please post it in https://stackoverflow.com or https://stats.stackexchange.com

I'll be glad to help or answer any question. I have spent a valuable amount of time developing this package. Please cite this package in your publication. Type 'citation("evola")' to know how to cite it.


}
\references{

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.

}
\keyword{ R package }
