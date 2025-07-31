\name{summary.Pop}
\alias{summary.Pop}
\title{summary form an evolafit model}
\description{\code{summary} method for class \code{"Pop"}.
}
\usage{
\method{summary}{Pop}(pop, \dots)
}
\arguments{
 \item{pop}{an object of class \code{"Pop"}}
 \item{\dots}{Further arguments to be passed}
}
\value{
  an updated model
 }
\author{Giovanny Covarrubias }
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
# Task: Gem selection. 
# Aim: Get highest combined value.
# Restriction: Max weight of the gem combined = 10. 

# simple specification
res<-evolafit(formula=cbind(Weight,Value)~Color, dt= Gems,
               # constraints on traits: if greater than this ignore
               constraintsUB = c(10,Inf), nGenerations = 2
)
summary(res$pop)
}
\keyword{models}

