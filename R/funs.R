evolafit <- function(formula, dt, 
                 constraintsUB, constraintsLB, traitWeight,
                 nCrosses=50, nProgeny=40,nGenerations=30, recombGens=1,
                 nQTLperInd=1, A=NULL, lambda=1,
                 propSelBetween=1,propSelWithin=0.5,
                 fitnessf=NULL){
  
  if(missing(formula)){stop("Please provide the formula to know traits and classifiers.", call. = FALSE)}
  elements <- strsplit(as.character(formula), split = "[+]")#[[1]]
  elements <- lapply(elements[-c(1)], function(x){all.vars(as.formula(paste("~",x)))})
  traits <- elements[[1]]
  # add the fitness function (current options are xa=-1 to ln>=0 )
  if(is.null(fitnessf)){fitnessf <- rep( list( function(xa, xtAx.lam, dt, pop, Q, alpha){xa - xtAx.lam} ), length(traits)) ; names(fitnessf) <- traits }else{
    if(is.list(fitnessf)){
      for(iTrait in traits){if(iTrait %in% names(fitnessf)){}else{fitnessf <- c(fitnessf, function(xa, xtAx.lam, dt, pop, Q, alpha){xa - xtAx.lam} ); names(fitnessf)[length(fitnessf)] <- iTrait}}
    }else{
      stop("The argument fitnessf should be a list of functions", call. = FALSE)
    }
  };  fitnessf <- fitnessf[traits]
  classifiers <- elements[[2]]
  if(length(constraintsUB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(length(constraintsLB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(length(traitWeight) != length(traits)){stop(paste0("Weights need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(length(lambda) != length(traits)){stop(paste0("Lambda need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(is.null(A)){A <- diag(nrow(dt))}
  
  # 1) initialize the population with customized haplotypes to ensure a single QTL per individual
  haplo = matrix(0, nrow=nrow(dt)*2, ncol = nrow(dt)) # rbind( diag(nrow(dt)), diag(nrow(dt)) )
  for (i in 1:nrow(haplo)) {
    haplo[i,sample(1:ncol(haplo), nQTLperInd )] <- 1
  }
  colnames(haplo) = paste0("H", 1:ncol(haplo))
  genMap = data.frame(markerName=colnames(haplo),
                      chromosome=1,
                      position=1:ncol(haplo))
  ped = data.frame(id=paste0("I", 1:nrow(dt)),
                   mother=0, father=0)
  founderPop = importHaplo(haplo=haplo, 
                           genMap=genMap,
                           ploidy=2L,
                           ped=ped)
  # founderPop = quickHaplo(nInd=Ne,nChr=1,segSites=nrow(dt), inbred = TRUE)
  SP = SimParam$new(founderPop)
  # 2) add the traits (columns from user) to take the values (rows) as marker effects
  for(iTrait in 1:length(traits)){
    SP$importTrait(markerNames = names(SP$genMap$`1`), addEff = dt[,traits[iTrait]]/2) # over 2 because QTL data is diplodized
  }
  alpha = do.call(cbind,lapply(SP$traits, function(x){x@addEff}))
  # 3) set the population
  pop = newPop(founderPop, simParam = SP)
  pop = randCross(pop, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)
  variances = diag(varG(pop))
  pop = setPheno(pop,h2=rep(.98,length(which(variances>0))), simParam = SP, traits = which(variances > 0) ) # ignore h2 since we will replace it in line 56
  # ***) creating the frame for the plot
  Performance <- matrix(NA, nrow=nGenerations,ncol=4) # to store results
  colnames(Performance) <- c("Average.xa","Best.xa","Average.xAx","nQTL.mu")
  rownames(Performance) <- paste("Generation",seq(nrow(Performance)))
  # 4) Starting the Generational process
  for (j in 1:nGenerations) { # for each generation we breed # j=1
    message(paste("generation",j))
    if(j > 1){
      ## apply selection between and within
      pop <- selectFam(pop=pop,nFam = round(nCrosses*propSelBetween), trait = selIndex, b=traitWeight, use = "pheno", simParam = SP)
      pop <- selectWithinFam(pop = pop, nInd = round(nProgeny*propSelWithin), trait = selIndex,  b=traitWeight, use = "pheno", simParam = SP)
      ## create new progeny
      for(k in 1:recombGens){
        # pop0<-pop
        pop <- randCross(pop=pop, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)
      }
      # if(carryParents){pop <- c(pop,pop0)}
      pop <- makeDH(pop=pop, nDH = 1, simParam = SP)
      ## compute constrained traits
      pop <- setPheno(pop=pop, h2=rep(0.98,length(which(variances>0))), simParam = SP, traits = which(variances > 0))  # ignore h2 since we will replace it in line 90
    }
    Q <- pullQtlGeno(pop, simParam = SP, trait = iTrait)  # plot((apply(Q,2,sum)/2)/nrow(Q))  # plot((apply(Q,1,sum)/2)/ncol(Q))
    for(iTrait in 1:length(traits)){ # iTrait=5
      xa <- Q %*% SP$traits[[iTrait]]@addEff
      if((max(xa)-min(xa)) > 0){
        xa = (xa-min(xa))/(max(xa)-min(xa)) # standardized xa
      }
      xtAx <- diag(Q%*%A%*%t(Q))
      if((max(xtAx)-min(xtAx)) > 0){
        xtAx = (xtAx-min(xtAx))/(max(xtAx)-min(xtAx)) # standardized xAx
      }
      xtAx.lam = xtAx * lambda[iTrait]
      pop@pheno[,iTrait] <- do.call(fitnessf[[iTrait]], list(dt=dt, xa=xa, xtAx.lam=xtAx.lam, pop=pop, Q=Q, alpha=alpha)) # xa - (lambda[iTrait] * xtAx ) # breeding value + coancestry
      # we only apply constraints to traits that will account for the total merit
      if(traitWeight[iTrait] != 0){ # if the trait will be used for selection
        ## pass each trait through all constraints
        for(iConUB in 1:length(constraintsUB)){ ## upper bound constraints ## for each trait constraint
          xaPrime <- Q %*% SP$traits[[iConUB]]@addEff
          if((max(xaPrime)-min(xaPrime)) > 0){ # there is variation
            ub = ( constraintsUB[iConUB]-min(xaPrime) )/(max(xaPrime)-min(xaPrime))
            xaPrime = (xaPrime-min(xaPrime))/(max(xaPrime)-min(xaPrime)) # standardized xa
          }else{
            ub=constraintsUB[iConUB]
          }
          constr <- xaPrime  #+ (lambda[iConUB] * xtAx )  # breeding value of the trait to be used for constraints
          pop@pheno[,iTrait] <- pop@pheno[,iTrait] * ifelse(constr[,1] > ub, ifelse(traitWeight[iTrait] > 0, -Inf, Inf) , 1) # adjust the BV of the trait after observing the constraints
          nan0 <- which(is.nan( pop@pheno[,iTrait]))
          if(length(nan0) > 0){ pop@pheno[nan0,iTrait] = ifelse(traitWeight[iTrait] > 0, -Inf, Inf) }
        }
        for(iConLB in 1:length(constraintsLB)){ # lower bound constraint
          xaPrime <- Q %*% SP$traits[[iConLB]]@addEff
          if((max(xaPrime)-min(xaPrime))  >0 ){
            lb = (constraintsLB[iConLB]- -min(xaPrime) )/(max(xaPrime)-min(xaPrime)) 
            xaPrime = (xaPrime-min(xaPrime))/(max(xaPrime)-min(xaPrime)) # standardized xa
          }else{
            lb=constraintsLB[iConLB]
          }
          constr <- xaPrime # + (lambda[iConLB] * xtAx )  # breeding value of the trait to be used for constraints
          pop@pheno[,iTrait] <- pop@pheno[,iTrait] * ifelse(constr[,1] < lb, ifelse(traitWeight[iTrait] > 0, -Inf, Inf) , 1) # adjust the BV of the trait after observing the constraints
          nan0 <- which(is.nan( pop@pheno[,iTrait]))
          if(length(nan0) > 0){ pop@pheno[nan0,iTrait] = ifelse(traitWeight[iTrait] > 0, -Inf, Inf) }
        } 
      }
    }
    #store the performance of the jth generation for plots
    xaFinal <- list()
    for(kk in 1:length(traits)){xaFinal[[kk]] <- Q %*% SP$traits[[1]]@addEff}
    score <- do.call(cbind, xaFinal) %*% traitWeight
    Performance[j,] <- c( mean(score), max(score) , mean(diag(Q%*%A%*%t(Q))),  mean(apply(Q/2,1,sum)) )
  }
  # 7) retrieve output of best combinations
  M <- pullQtlGeno(pop, simParam = SP, trait=1); M <- M/2
  colnames(M) <- apply(data.frame(dt[,classifiers]),1,function(x){paste(x,collapse = "_")})
  return(list(M=M, score=Performance, pheno=pop@pheno, pop=pop))
}


pmonitor <- function(object,...){
  x <- object$score#[,"Best.xa"]
  mmin <- min(c(0, x[,"Average.xa"]))
  mmax <- max(x[,"Best.xa"])
  plot(x[1,], type="o", ylim=c(mmin,mmax), xlab="Generation", ylab="Value")
  if(ncol(x) > 1){
    for(i in 2:ncol(x)){
      par(new=TRUE)
      plot(x[i,], col=i, ylim=c(mmin,mmax), ylab="",xlab="", type="o",...)
    }
  }
  legend("topright",legend = colnames(x), col=1:(ncol(x)), bty="n", lty=1)
}

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