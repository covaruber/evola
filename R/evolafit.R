evolafit <- function(formula, dt, 
                     constraintsUB, constraintsLB, traitWeight,
                     nCrosses=50, nProgeny=40,nGenerations=30, recombGens=1,
                     nQTLperInd=NULL, A=NULL, lambda=NULL,
                     propSelBetween=1,propSelWithin=0.5,
                     fitnessf=NULL, verbose=TRUE, dateWarning=TRUE,
                     selectTop=TRUE, ...){
  
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
  if(is.null(fitnessf)){
    fitnessf <-function (Y, b, d, scale = FALSE) {
      if (scale) {
        return(scale(Y) %*% b - d)
      }
      return(Y %*% b - d)
    }
  };  #fitnessf <- fitnessf[traits]
  classifiers <- elements[[2]]
  # check that the user has provided a single value for each QTL
  checkNQtls <- table(dt[,classifiers])
  if(length(which(checkNQtls > 1)) > 0){stop("You cannot provide more than one alpha value per QTL. Make sure that your x variable has only one value.", call. = FALSE)}
  # print(classifiers)
  if(missing(constraintsUB)){constraintsUB <- rep(Inf,length(traits))}
  if(length(constraintsUB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(missing(constraintsLB)){constraintsLB <- rep(-Inf,length(traits))}
  if(length(constraintsLB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(missing(traitWeight)){traitWeight <- rep(1,length(traits))}
  if(length(traitWeight) != length(traits)){stop(paste0("Weights need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(is.null(lambda)){lambda <- 0}
  # if(length(lambda) != length(traits)){stop(paste0("Lambda need to have the same length than traits (",length(traits),")"), call. = FALSE)}
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
  ################################
  ################################
  ################################
  ## FOR GENERATION
  for (j in 1:nGenerations) { # for each generation we breed # j=1
    if(verbose){message(paste("generation",j))}
    if(j > 1){
      # group relationship
      xtAx <- Matrix::diag(Q%*%Matrix::tcrossprod(A,Q))
      # calculate base coancestry Ct
      m <- Matrix::Matrix(1,nrow=1,ncol=ncol(A))
      mAmt <- as.vector((m%*%Matrix::tcrossprod(A,m))/(4*(ncol(A)^2)))
      # rate of coancestry xtAx/4p^2 - mtAm/4n^2
      deltaC <- ( (xtAx/(4*(apply(Q/2,1,sum)^2))) - mAmt)/(1-mAmt)
      # if there is variation in min and max values in xtAx standardize
      if((max(xtAx)-min(xtAx)) > 0){ 
        xtAx = (xtAx-min(xtAx))/(max(xtAx)-min(xtAx)) # standardized xAx
      }
      ## apply selection between and within
      xtAx.lam = xtAx * lambda#[iTrait]
      names(xtAx.lam) <- pop@id
      if(j == 2){
        suppressWarnings( best <- selectInd(pop=pop, nInd =min(c(10,nInd(pop))), trait = fitnessf,  
                          b=traitWeight, d=xtAx.lam[pop@id], use = "pheno", simParam = SP, selectTop=selectTop, ...), classes = "warning")
      }else{
        best <- c(best, suppressWarnings( selectInd(pop=pop, nInd =min(c(10,nInd(pop))), trait = fitnessf,  
                                  b=traitWeight, d=xtAx.lam[pop@id], use = "pheno", simParam = SP, selectTop=selectTop, ...), classes = "warning")  )
       
      }
      suppressWarnings( popF <- selectFam(pop=pop,nFam = round(nCrosses*propSelBetween), trait = fitnessf, 
                        b=traitWeight,d=xtAx.lam[pop@id], use = "pheno", simParam = SP,selectTop=selectTop,...), classes = "warning")
      suppressWarnings( popW <- selectWithinFam(pop = pop, nInd = round(nProgeny*propSelWithin), 
                              trait = fitnessf,  b=traitWeight,d=xtAx.lam[pop@id], use = "pheno", simParam = SP,
                              selectTop=selectTop,...), classes = "warning")
      selected <- intersect(popF@id,popW@id)
      pop <- pop[which(pop@id %in% selected)]
      ## create new progeny
      for(k in 1:recombGens){
        pop <- randCross(pop=pop, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)
      }
      pop <- makeDH(pop=pop, nDH = 1, simParam = SP)
      ## compute constrained traits
      pop <- setPheno(pop=pop, h2=rep(0.98,length(which(variances>0))), simParam = SP, traits = which(variances > 0))  # ignore h2 since we will replace it in line 90
    }
    # extract solutions for the trait 1 because all traits have the same QTLs
    Q <- pullQtlGeno(pop, simParam = SP, trait = 1)  #
    Q <- as(Q, Class = "dgCMatrix")
    constCheckUB <- constCheckLB <- matrix(1, nrow=nrow(Q), ncol=length(traits))
    xaFinal <- list()
    ################################
    ################################
    ## FOR EACH TRAIT
    for(iTrait in 1:length(traits)){ # iTrait=1
      xaOr <- Q %*% SP$traits[[iTrait]]@addEff # solutions * alpha for the iTrait
      xaFinal[[iTrait]] <- xaOr
      if((max(xaOr)-min(xaOr)) > 0){ # if there is variation
        xa = (xaOr-min(xaOr))/(max(xaOr)-min(xaOr)) # standardized xa
      }
      # calculate the genetic value of solutions using the objective functions
      pop@pheno[,iTrait] <- xa[,1] # as.vector( do.call(fitnessf[[iTrait]], list(dt=dt, xa=xa, xtAx.lam=xtAx.lam, pop=pop, Q=Q, alpha=alpha,...)) ) # xa - (lambda[iTrait] * xtAx ) # breeding value + coancestry
      
      # check the contraints and trace them back
      constCheckUB[,iTrait] <- ifelse( (xaOr[,1] > constraintsUB[iTrait])  , 0 , 1) # ifelse(c1+c2 < 2, 0, 1)
      constCheckLB[,iTrait] <- ifelse( (xaOr[,1] < constraintsLB[iTrait]) , 0 , 1)
      nan0 <- which(is.nan( pop@pheno[,iTrait]))
      if(length(nan0) > 0){ constCheckUB[nan0,iTrait] = 0; constCheckLB[nan0,iTrait] = 0 }
    } # end of for each trait
    ################################
    ################################
    # sum of how many trait constraints are met
    metConstCheck <- apply(constCheckUB,1,sum) 
    didntMetConst <- which(metConstCheck < length(traits))
    # remove individuals that break the constraints UB
    if(length(didntMetConst)>0){ # 
      for(iTrait in 1:length(traits)){
        pop <- pop[setdiff(1:nInd(pop),didntMetConst)]
      }
    }
    #  remove individuals that break the constraints LB
    metConstCheck <- apply(constCheckLB,1,sum) 
    didntMetConst <- which(metConstCheck < length(traits))
    # impute with mean value the ones that do not met the constraints
    if(length(didntMetConst)>0){
      for(iTrait in 1:length(traits)){
        pop <- pop[setdiff(1:nInd(pop),didntMetConst)]
      }
    }
    #store the performance of the jth generation for plot functions
    score <- do.call(cbind, xaFinal) %*% traitWeight
    
    if(j > 1){
      indivPerformance[[j]] <- data.frame(score=as.vector(score),deltaC= as.vector(deltaC) , xtAx= as.vector(xtAx), generation=j, nQTL=apply(Q/2,1,sum)) # save individual solution performance
      averagePerformance[j,] <- c( mean(score,na.rm=TRUE), max(score,na.rm=TRUE) , mean(xtAx,na.rm=TRUE),  mean(apply(Q/2,1,sum),na.rm=TRUE), mean(deltaC,na.rm=TRUE) ) # save summaries of performance
    }
  }# end of for each generation
  ################################
  ################################
  ################################
  # 7) retrieve output of best combinations
  M <- pullQtlGeno(pop, simParam = SP, trait=1); M <- M/2
  Mb <- pullQtlGeno(best, simParam = SP, trait=1); Mb <- Mb/2
  colnames(M) <- apply(data.frame(dt[,classifiers]),1,function(x){paste(x,collapse = "_")})
  indivPerformance <- do.call(rbind, indivPerformance)
  return(list(M=M, Mb=Mb, score=averagePerformance, pheno=pop@pheno,phenoBest=best@pheno, pop=pop, best=best, indivPerformance=indivPerformance, constCheckUB=constCheckUB, constCheckLB=constCheckLB, traits=traits))
}

