
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
    packageStartupMessage(blue(paste("[] Evolutionary Algorithm in R (evola) 1.0.7 (2025-12)              []",sep="")),appendLF=TRUE)
    packageStartupMessage(paste0(blue("[] Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen(white(" ")), bgWhite(magenta("M")), bgRed(white(" ")),"  ", bgRed(bold(yellow(" (") )),bgRed(bold(white("W"))), bgRed(bold(yellow(") "))) ) ,"                 []")),appendLF=TRUE)
    packageStartupMessage(blue("[] Dedicated to the University of Chapingo and UW-Madison           []"),appendLF=TRUE)
    packageStartupMessage(blue("[] Type 'vignette('evola.intro')' for a short tutorial              []"),appendLF=TRUE)
    packageStartupMessage(blue(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(blue("evola is updated on CRAN every 4-months due to CRAN policies"),appendLF=TRUE)
    packageStartupMessage(blue("Source code is available at https://github.com/covaruber/evola"),appendLF=TRUE)
  }
  invisible()
}

.onLoad <- function(library, pkg){
  Sys.setenv("OMP_THREAD_LIMIT"=2)
}
##################################################################################################
##################################################################################################



##################################################################################################
##################################################################################################

addZeros<- function (x, nr=2) {
  ## left
  nz <- nchar(round(max(x)))
  add <- abs(nchar(round(x)) - nz)
  toAdd <- gsub("1", "", ifelse(add > 0, 10^(add), ""))
  newX <- paste0(toAdd, as.character(x))
  newX <- gsub("\\..*","",newX)
  ## right
  if(nr > 0){
    toAdd <- format( round(x-floor(x), nr) , nsmall=nr) 
    badCalls <- grep("1[.]",toAdd)
    if(length(badCalls) > 0){toAdd[badCalls] <- paste0("0.",paste(rep("9",nr), collapse="") )}
    newX2 <- paste(newX, gsub("0[.]","",toAdd) , sep="." )
  }else{
    newX2 <- newX
  }
  return(newX2)
}

ocsFun <- function (Y, b, Q, D, lambda, scaled=TRUE, ...) {
  # (q'a)b - l(q'Dq)
  if(scaled){
    Yb <- apply(Y,2,function(x){
      if(var(x)>0){return(scale(x))}else{return(x*0)}
    }) %*% b
  }else{
    Yb <- Y %*% b
  }
  QtDQ <- Matrix::diag(Q %*% Matrix::tcrossprod(D, Q))
  if(var(Yb) > 0){Yb <- stan(Yb)}
  if(var(QtDQ) > 0){QtDQ <- stan(QtDQ)}
  return( Yb -  lambda*QtDQ )
}

ocsFunC <- function (Q, 
                     SNP, solution, #alphaLog=1, 
                     wtf=NULL, wbaf='base', ...) {
  # wtf is weights for trait frequencies ;) 

  if(!missing(solution)){
    if(is.null(solution)){stop("This function cannot have a solution argument as NULL. Please correct.", call. = FALSE)}
    if(inherits(solution,"RRsol")){
      solution <- do.call(cbind, lapply(solution@bv,function(x){x@addEff}))
    }
  }
  # allele frequency breeding value
  total_alleles <- 2 * colSums(!is.na(SNP))
  # Alternate allele count: 2 * homozygous alt (2) + 1 * heterozygous (1)
  alt_alleles <- colSums(SNP, na.rm = TRUE)
  # Frequencies of the base population
  freq_alt <- alt_alleles / total_alleles
  nans <- which(is.nan(freq_alt))
  if(length(nans) > 0){freq_alt[nans]=0}
  freq_ref <- 1 - freq_alt
  # Combine into matrix
  traitFreqs <- list()
  for(iTrait in 1:ncol(solution)){ # iTrait=1
    freqsPos <- ifelse(solution[,iTrait]>0,freq_alt, freq_ref)
    if(wbaf == 'alpha'){ # we use allelic effects as weights for the frequencies
      traitFreqs[[iTrait]] = apply(Q,1,function(x){
        sum(freqPosAllele(SNP[which(x>0),,drop=FALSE], alpha = solution[, iTrait]) * abs(solution[, iTrait]), na.rm=TRUE )
      })
    }else if(wbaf == 'base'){ # we use the basic version of 1-freq.base
      traitFreqs[[iTrait]] = apply(Q,1,function(x){
        sum(freqPosAllele(SNP[which(x>0),,drop=FALSE], alpha = solution[, iTrait]) * (1-freqsPos), na.rm=TRUE )
      })
    }else if(wbaf == 'none'){ # no weights
      traitFreqs[[iTrait]] = apply(Q,1,function(x){
        sum(freqPosAllele(SNP[which(x>0),,drop=FALSE], alpha = solution[, iTrait]) , na.rm=TRUE )
      })
    }else{
      stop("Method not available", call. = FALSE)
    }
    # traitFreqs[[iTrait]] = Q%*%SNP%*%sign(solution[,iTrait])%*%(1-freqsPos)
  }
  Y2 = do.call(cbind, traitFreqs)
  if(is.null(wtf)){
    wtf = rep(1,ncol(Y2))
  }
  Yb2 = Y2 %*% wtf # all trait frequencies are equally important
  
  return(Yb2)
}

regFun <- function ( Q, a, X, y, ...) {
  n <- ncol(X)
  p <- apply(Q, 1, function(z) {
    which(z > 0)
  })
  if (is.matrix(p)) {
    p <- lapply(seq_len(ncol(p)), function(i) p[, i])
  }
  nq <- unlist(lapply(p, length))
  v <- 1:nrow(X)
  mse = vector("numeric", length(v))
  for (j in v) {
    if (nq[j] == n) {
      mse[j] = sum(((y[v]) - (as.matrix(X[v, , drop = FALSE]) %*% 
                                a[ p[[j]], 1 ] ))^2)
    }
    else {
      mse[j] = Inf
    }
  }
  return(mse)
}

inbFun <- function (Q, D, ...) {
  return( Matrix::diag(Q %*% Matrix::tcrossprod(D, Q))  )
}

varQ <- function(object){
  nTraits <- length(object$pop@traits)
  n <- numeric()
  for(iTrait in 1:nTraits){
    Q <- pullQtlGeno(object$pop, simParam = object$simParam, trait=iTrait); Q <- Q/2
    n[iTrait] <- sum(apply(Q,2,var, na.rm=TRUE))
  }
  names(n) <- object$pop@traits
  return(n)
}

nQtl <- function(object){
  nTraits <- length(object$pop@traits)
  n <- matrix(NA, nrow = nInd(object$pop), ncol = nTraits)
  for(iTrait in 1:nTraits){
    Q <- pullQtlGeno(object$pop, simParam = object$simParam, trait=iTrait); Q <- Q/2
    n[,iTrait] <- apply(Q, 1, function(x){length(which(x > 0))})
  }
  colnames(n) <- object$pop@traits
  rownames(n) <- object$pop@id
  return(n)
}

freqPosAllele <- function(M, alpha){
  # M should not be centered
  total_alleles <- 2 * colSums(!is.na(M))
  # Alternate allele count: 2 * homozygous alt (2) + 1 * heterozygous (1)
  alt_alleles <- colSums(M, na.rm = TRUE)
  # Frequencies
  freq_alt <- alt_alleles / total_alleles
  nans <- which(is.nan(freq_alt))
  if(length(nans) > 0){freq_alt[nans]=0}
  freq_ref <- 1 - freq_alt
  # Combine into matrix
  freqsPos <- ifelse(alpha>0,freq_alt, freq_ref)
  return(freqsPos)
}

Jc <- function(nc){
  matrix(1,nrow=1,ncol=nc)
}

Jr <- function(nr){
  matrix(1,nrow=nr,ncol=1)
}

bestSol <- function (object, selectTop = TRUE, n = 1){
  if (!inherits(object, c("Pop", "evolaMod"))) {
    stop("Object of type Pop or evolaMod expected", call. = FALSE)
  }
  if (nInd(object) > 0) {
    dd = as.data.frame(cbind(object@gv, object@fitness))
    rownames(dd) <- object@id
    picked <- list()
    if (selectTop) {
      for (i in 1:ncol(dd)) {
        dd = dd[order(-dd[, i]), , drop = FALSE]
        selected = rownames(dd)[1:n]
        picked[[i]] <- which(object@id %in% selected)
      }
    }
    else {
      for (i in 1:ncol(dd)) {
        dd = dd[order(dd[, i]), , drop = FALSE]
        selected = rownames(dd)[1:n]
        picked[[i]] <- which(object@id %in% selected)
      }
    }
    res1 <- do.call(cbind, picked)
    colnames(res1) <- c(object@traits, "fitness")
    return(res1)
  }
  else {
    stop("No individuals in the object provided")
  }
}

importHaploSparse <- function (haplo, genMap, ploidy = 2L, ped = NULL) 
{
  if (!is.null(ped)) {
    if (is.vector(ped)) {
      id = as.character(ped)
      stopifnot(length(id) == (nrow(haplo)/ploidy), !any(duplicated(id)))
      mother = father = rep("0", length(id))
    }
    else {
      id = as.character(ped[, 1])
      stopifnot(length(id) == (nrow(haplo)/ploidy), !any(duplicated(id)))
      mother = as.character(ped[, 2])
      father = as.character(ped[, 3])
    }
  }
  genMap = importGenMap(genMap)
  if (is.data.frame(haplo)) {
    haplo = as.matrix(haplo)
  }
  markerName = colnames(haplo)
  # haplo = matrix(as.raw(haplo), ncol = ncol(haplo))
  # stopifnot(haplo == as.raw(0) | haplo == as.raw(1))
  haplotypes = vector("list", length = length(genMap))
  for (i in seq_len(length(genMap))) {
    mapMarkers = names(genMap[[i]])
    take = match(mapMarkers, markerName)
    if (any(is.na(take))) {
      genMap[[i]] = genMap[[i]][is.na(take)]
      stopifnot(length(genMap[[i]]) >= 1L)
      genMap[[i]] = genMap[[i]] - genMap[[i]] - genMap[[i]][1]
      take = na.omit(take)
    }
    haplotypes[[i]] = haplo[, take, drop = FALSE]
  }
  founderPop = newMapPop(genMap = genMap, haplotypes = haplotypes, 
                         ploidy = ploidy)
  if (!is.null(ped)) {
    founderPop = new("NamedMapPop", id = id, mother = mother, 
                     father = father, founderPop)
  }
  return(founderPop)
}

drift <- function(pop, simParam, solution=NULL, traits=1){
  # traits <- 1:simParam$nTraits
  currentFreqPositive <- list()
  counter=1
  for(iTrait in traits){ # iTrait=1
    if(is.null(solution)){
      alpha = simParam$traits[[iTrait]]@addEff
      Qtl = pullQtlGeno(pop, simParam = simParam, trait = iTrait)
    }else{ # user provided a solution model
      alpha = solution@gv[[iTrait]]@addEff
      Qtl = pullSnpGeno(pop, simParam = simParam)
    }
    m = matrix(0,nrow=1,ncol=3); colnames(m) <- c(0,1,2)
    freqsG = apply(Qtl,2,function(x){
      tt = table(x)
      m[,names(tt)] = tt
      return(m)
    })
    freqsG = freqsG/apply(freqsG,2,sum)
    rownames(freqsG) <- 0:2
    freqsA = apply(freqsG,2, function(x){
      matrix(c( x[1]+(0.5*x[2]), x[3]+(0.5*x[2]) ), nrow = 1, ncol=2)
    })
    rownames(freqsA) <- c(0,2)
    desiredAllele = ifelse( sign(alpha) > 0 , 2, 0 )
    prov <- sapply(1:length(desiredAllele), function(x){freqsA[as.character(desiredAllele[x]) ,x]})
    
    if(is.null(solution)){
      out <- cbind(getQtlMap(trait = iTrait, simParam=simParam), desiredAllele,prov)
    }else{
      out <- cbind( getSnpMap(snpChip = 1, simParam = simParam), desiredAllele,prov)
    }
    rownames(out) <- colnames(Qtl) ; colnames(out) <- c("id","chr","ss","pos","a+","freq")
    currentFreqPositive[[counter]] <- out
    counter <- counter+1
  }
  names(currentFreqPositive) <- simParam$traitNames[traits]
  return(currentFreqPositive)
}

importHaplo = function(haplo, genMap, ploidy=2L, ped=NULL){
  # Extract pedigree, if supplied
  if(!is.null(ped)){
    if(is.vector(ped)){
      id = as.character(ped)
      stopifnot(length(id)==(nrow(haplo)/ploidy),
                !any(duplicated(id)))
      mother = father = rep("0", length(id))
    }else{
      id = as.character(ped[,1])
      stopifnot(length(id)==(nrow(haplo)/ploidy),
                !any(duplicated(id)))
      mother = as.character(ped[,2])
      father = as.character(ped[,3])
    }
  }
  
  genMap = importGenMap(genMap)
  
  # Get marker names
  if(is.data.frame(haplo)){
    haplo = as.matrix(haplo)
  }
  markerName = colnames(haplo)
  
  # Convert haplotypes to raw
  haplo = matrix(as.raw(haplo), ncol=ncol(haplo))
  stopifnot(haplo==as.raw(0) | haplo==as.raw(1))
  
  # Create haplotype list
  haplotypes = vector("list", length=length(genMap))
  
  # Order haplotypes by chromosome
  for(i in seq_len(length(genMap))){
    mapMarkers = names(genMap[[i]])
    take = match(mapMarkers, markerName)
    if(any(is.na(take))){
      genMap[[i]] = genMap[[i]][is.na(take)]
      stopifnot(length(genMap[[i]]) >= 1L)
      genMap[[i]] = genMap[[i]] - genMap[[i]]-genMap[[i]][1]
      take = na.omit(take)
    }
    haplotypes[[i]] = haplo[,take,drop=FALSE]
  }
  
  founderPop = newMapPop(genMap=genMap,
                         haplotypes=haplotypes,
                         ploidy=ploidy)
  
  if(!is.null(ped)){
    founderPop = new("NamedMapPop",
                     id=id,
                     mother=mother,
                     father=father,
                     founderPop)
  }
  
  return(founderPop)
}
