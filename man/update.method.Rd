\name{update.evolaFitMod}
\alias{update.evolaFitMod}
\title{update form an evolafit model}
\description{\code{update} method for class \code{"evolaFitMod"}.
}
\usage{
\method{update}{evolaFitMod}(object, formula., evaluate = TRUE, \dots)
}
\arguments{
 \item{object}{an object of class \code{"evolaFitMod"}}
 \item{formula.}{an optional formula}
 \item{evaluate}{a logical value to indicate if an evaluation of the call should be done or not}
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
resUp=update(res)
}
\keyword{models}

