evolafit <- function(formula, dt, 
                 constraintsUB, constraintsLB, traitWeight,
                 nCrosses=50, nProgeny=40,nGenerations=30, recombGens=1,
                 nQTLperInd=NULL, A=NULL, lambda=NULL,
                 propSelBetween=1,propSelWithin=0.5,
                 fitnessf=NULL, verbose=TRUE, dateWarning=TRUE){
  
  my.date <- "2024-11-01"
  your.date <- Sys.Date()
  ## if your month is greater than my month you are outdated
  if(dateWarning & verbose){
    if (your.date > my.date) {
      warning("Version out of date. Please update evola to the newest version using:\ninstall.packages('evola') in a new session\n Use the 'dateWarning' argument to disable the warning message.", call. = FALSE)
    }
  }
  
  if(missing(formula)){stop("Please provide the formula to know traits and classifiers.", call. = FALSE)}
  elements <- strsplit(as.character(formula), split = "[+]")#[[1]]
  elements <- lapply(elements[-c(1)], function(x){all.vars(as.formula(paste("~",x)))})
  traits <- elements[[1]]
  if(!all(traits%in%colnames(dt))){stop("Specified traits are not traits in the dataset. Please correct.", call. = FALSE)}
  # add the fitness function (current options are xa=-1 to ln>=0 )
  if(is.null(fitnessf)){fitnessf <- rep( list( function(xa, xtAx.lam, dt, pop, Q, alpha){xa - xtAx.lam} ), length(traits)) ; names(fitnessf) <- traits }else{
    if(is.list(fitnessf)){
      for(iTrait in traits){if(iTrait %in% names(fitnessf)){}else{fitnessf <- c(fitnessf, function(xa, xtAx.lam, dt, pop, Q, alpha){xa - xtAx.lam} ); names(fitnessf)[length(fitnessf)] <- iTrait}}
    }else{
      stop("The argument fitnessf should be a list of functions", call. = FALSE)
    }
  };  fitnessf <- fitnessf[traits]
  classifiers <- elements[[2]]
  # print(classifiers)
  if(missing(constraintsUB)){constraintsUB <- rep(Inf,length(traits))}
  if(length(constraintsUB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(missing(constraintsLB)){constraintsLB <- rep(-Inf,length(traits))}
  if(length(constraintsLB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(missing(traitWeight)){traitWeight <- rep(1,length(traits))}
  if(length(traitWeight) != length(traits)){stop(paste0("Weights need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(is.null(lambda)){lambda <- rep(0, length(traits))}
  if(length(lambda) != length(traits)){stop(paste0("Lambda need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(is.null(A)){A <- Matrix::Diagonal(nrow(dt))}
  if(is.null(nQTLperInd)){nQTLperInd <- nrow(dt)/5}
  
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
  indivPerformance <- list() # store results by generation
  averagePerformance <- matrix(NA, nrow=nGenerations,ncol=5) # to store results
  colnames(averagePerformance) <- c("Average.xa","Best.xa","Average.xAx","nQTL.mu", "deltaC.mu")
  rownames(averagePerformance) <- paste("Generation",seq(nrow(averagePerformance)))
  # 4) Starting the Generational process
  for (j in 1:nGenerations) { # for each generation we breed # j=1
    if(verbose){message(paste("generation",j))}
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
    Q <- as(Q, Class = "dgCMatrix")
    # print(str(Q))
    # print(str(A%*%t(Q)))
    # print(apply(Q,1,sum))
    # print(t(Q))
    xtAx <- Matrix::diag(Q%*%Matrix::tcrossprod(A,Q))
    # print(xtAx)
    # xtAx <- diag(Q%*%tcrossprod(A,Q))
    # print(xtAx)
    # calculate base coancestry Ct
    m <- Matrix::Matrix(1,nrow=1,ncol=ncol(A))
    mAmt <- as.vector((m%*%Matrix::tcrossprod(A,m))/(4*(ncol(A)^2)))
    # rate of coancestry xtAx/4p^2 - mtAm/4n^2
    deltaC <- ( (xtAx/(4*(apply(Q/2,1,sum)^2))) - mAmt)/(1-mAmt)
    # if there is variation in min and max values in xtAx standardize
    if((max(xtAx)-min(xtAx)) > 0){ 
      xtAx = (xtAx-min(xtAx))/(max(xtAx)-min(xtAx)) # standardized xAx
    }
    for(iTrait in 1:length(traits)){ # iTrait=1
      xa <- Q %*% SP$traits[[iTrait]]@addEff
      if((max(xa)-min(xa)) > 0){ # if there is variation
        xa = (xa-min(xa))/(max(xa)-min(xa)) # standardized xa
      }
      xtAx.lam = xtAx * lambda[iTrait]
      pop@pheno[,iTrait] <- as.vector( do.call(fitnessf[[iTrait]], list(dt=dt, xa=xa, xtAx.lam=xtAx.lam, pop=pop, Q=Q, alpha=alpha)) ) # xa - (lambda[iTrait] * xtAx ) # breeding value + coancestry
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
    #store the performance of the jth generation for plot functions
    xaFinal <- list()
    for(kk in 1:length(traits)){ # save merit of the different solutions for each trait
      xaFinal[[kk]] <- Q %*% SP$traits[[1]]@addEff
    }
    score <- do.call(cbind, xaFinal) %*% traitWeight
    indivPerformance[[j]] <- data.frame(score=as.vector(score),deltaC=as.vector(deltaC), xtAx=as.vector(xtAx), generation=j, nQTL=apply(Q/2,1,sum)) # save individual solution performance
    averagePerformance[j,] <- c( mean(score,na.rm=TRUE), max(score,na.rm=TRUE) , mean(xtAx,na.rm=TRUE),  mean(apply(Q/2,1,sum),na.rm=TRUE), mean(deltaC,na.rm=TRUE) ) # save summaries of performance
  }
  # 7) retrieve output of best combinations
  M <- pullQtlGeno(pop, simParam = SP, trait=1); M <- M/2
  colnames(M) <- apply(data.frame(dt[,classifiers]),1,function(x){paste(x,collapse = "_")})
  indivPerformance <- do.call(rbind, indivPerformance)
  return(list(M=M, score=averagePerformance, pheno=pop@pheno, pop=pop, indivPerformance=indivPerformance))
}