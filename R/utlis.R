
##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.5.0"))
    stop("This package requires R 3.5.0 or later")
  if(interactive()) {
    packageStartupMessage(blue(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(blue(paste("[] Evolutionary Algorithm in R (evola) 1.0.1 (2024-07)              []",sep="")),appendLF=TRUE)
    packageStartupMessage(paste0(blue("[] Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen(white(" ")), bgWhite(magenta("*")), bgRed(white(" "))),"                        []")),appendLF=TRUE)
    packageStartupMessage(blue("[] Dedicated to the University of Chapingo and UW-Madison           []"),appendLF=TRUE)
    # packageStartupMessage(blue("[] Type 'vignette('evola.intro')' for a short tutorial             []"),appendLF=TRUE)
    packageStartupMessage(blue(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(blue("evola is updated on CRAN every 4-months due to CRAN policies"),appendLF=TRUE)
    packageStartupMessage(blue("Source code is available at https://github.com/covaruber/evola"),appendLF=TRUE)
  }
  invisible()
}

##################################################################################################
##################################################################################################

varM <- function(object){
  sum(apply(object$M,2,var, na.rm=TRUE))
}

bestSol <- function(object){
  apply(object$pheno,2,function(y){
    y <- y[which( y < Inf)]
    y <- y[which( y > -Inf)]
    best = which(y==max(y,na.rm=TRUE))[1]
    return(best)
  })
} 

A.mat <- function (X, min.MAF = NULL) 
{
  
  X <- as.matrix(X)
  n <- nrow(X)
  frac.missing <- apply(X, 2, function(x) {
    length(which(is.na(x)))/n
  })
  missing <- max(frac.missing) > 0
  freq <- apply(X + 1, 2, function(x) {
    mean(x, na.rm = missing)
  })/2
  MAF <- apply(rbind(freq, 1 - freq), 2, min)
  if (is.null(min.MAF)) {
    min.MAF <- 1/(2 * n)
  }
  max.missing <- 1 - 1/(2 * n)
  markers <- which((MAF >= min.MAF) & (frac.missing <= max.missing))
  m <- length(markers)
  var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
  one <- matrix(1, n, 1)
  mono <- which(freq * (1 - freq) == 0)
  X[, mono] <- 2 * tcrossprod(one, matrix(freq[mono], length(mono), 
                                          1)) - 1
  freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
  W <- X[, markers] + 1 - 2 * freq.mat
  A <- tcrossprod(W)/var.A/m
  return(A)
}
