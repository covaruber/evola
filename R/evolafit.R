evolafit <- function(formula, dt, 
                     constraintsUB, constraintsLB,constraintW=NULL, 
                     b, nCrosses=50, nProgeny=20,nGenerations=20, 
                     recombGens=1, nChr=1, mutRate=0,
                     nQTLperInd=NULL, D=NULL, lambda=0,
                     propSelBetween=NULL,propSelWithin=NULL,
                     fitnessf=NULL, verbose=TRUE, dateWarning=TRUE,
                     selectTop=TRUE, tolVarG=1e-6, 
                     Ne=50, initPop=NULL, simParam = NULL, 
                     fixQTLperInd=FALSE, traceDelta=TRUE, topN=10, ...){
  
  my.date <- "2025-08-01"
  your.date <- Sys.Date()
  ## if your month is greater than my month you are outdated
  if(dateWarning & verbose){
    if (your.date > my.date) {
      warning("Version out of date. Please update evola to the newest version using:\ninstall.packages('evola') in a new session\n Use the 'dateWarning' argument to disable the warning message.", call. = FALSE)
    }
  }
  # if(propSelBetween==0 | propSelWithin==0){stop("Please ensure that parameters propSelWithin and propSelBetween are different than zero.", call. = FALSE)}
  if(missing(formula)){stop("Please provide the formula to know traits and classifiers.", call. = FALSE)}
  if(is.null(propSelBetween)){
    propSelBetween <- stan(logspace(seq(1,-1, -2/nGenerations), p=3), ub=0.8, lb=0.2)
  }else{
    propSelBetween <- rep(propSelBetween, nGenerations)
  }
  if(is.null(propSelWithin)){
    propSelWithin <- stan(logspace(seq(1,-1, -2/nGenerations), p=3), ub=0.2, lb=0.8)
  }else{
    propSelWithin <- rep(propSelWithin, nGenerations)
  }
  if(propSelBetween[1]>0 & nCrosses==0){stop("If you apply selection between families you need to set nCrosses to a value > 0.", call. = FALSE)}
  if(nGenerations < 1){stop("nGenerations cannot be smaller than 1.")}
  if(is.null(constraintW)){
    constraintW <- rep(1, nGenerations)
  } # different weight to constraints at each generation
  
  mc <- match.call() # create a call
  # add the fitness function (current options are qa=-1 to ln>=0 )
  if(is.null(fitnessf)){
    fitnessf <- ocsFun
  };
  # get the name of the traits
  elements <- strsplit(as.character(formula), split = "[+]")#[[1]]
  elements <- lapply(elements[-c(1)], function(x){all.vars(as.formula(paste("~",x)))})
  traits <- elements[[1]]
  if(!all(traits%in%colnames(dt))){stop("Specified traits are not traits in the dataset. Please correct.", call. = FALSE)}
  classifiers <- elements[[2]]
  checkNQtls <- table(dt[,classifiers])
  if(length(which(checkNQtls > 1)) > 0){stop("You cannot provide more than one alpha value per QTL. Make sure that your x variable has only one value.", call. = FALSE)}
  if(missing(constraintsUB)){constraintsUB <- rep(Inf,length(traits))}
  if(length(constraintsUB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(missing(constraintsLB)){constraintsLB <- rep(-Inf,length(traits))}
  if(length(constraintsLB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(missing(b)){b <- rep(1,length(traits))}
  if(length(b) != length(traits)){stop(paste0("Weights need to have the same length than traits (",length(traits),")"), call. = FALSE)}
  if(is.null(D)){D <- Matrix::Diagonal(nrow(dt)); useD=FALSE}else{useD=TRUE}
  if(is.null(nQTLperInd)){nQTLperInd <- nrow(dt)/5}
  # check that the user has provided a single value for each QTL
  nMutations = round(mutRate * nrow(dt)) # number of mutations per individual per generation
  
  if(is.null(initPop)){
    # 1) initialize the population with customized haplotypes to ensure a single QTL per individual
    av <- 1:nrow(dt)
    haplo = Matrix::Matrix(0, nrow= Ne*2, ncol = nrow(dt)) # rbind( diag(nrow(dt)), diag(nrow(dt)) )
    for (i in seq(1,nrow(haplo),2)) {
      haplo[i,sample(av,nQTLperInd)] <- 1
      # haplo[i,] <- ifelse(runif(ncol(haplo))< (nQTLperInd/ncol(haplo)) ,1,0)
      haplo[(i+1),] <- haplo[i,]
    }
    colnames(haplo) = dt[,classifiers]
    
    nQtlPerChr = rep( floor(length(colnames(haplo))/nChr), nChr)
    nQtlPerChr[length(nQtlPerChr)] = nQtlPerChr[length(nQtlPerChr)] + ( length(colnames(haplo)) - sum(nQtlPerChr) )
    
    chromosome = as.vector(unlist(apply( data.frame(nQtlPerChr,1:nChr), 1, function(x){rep(x[2],x[1])})))
    position = as.vector(unlist(apply( data.frame(nQtlPerChr,1:nChr), 1, function(x){1:x[1]})))
    genMap = data.frame(markerName=colnames(haplo),
                        chromosome=chromosome,
                        position=position)
    ped = data.frame(id=paste0("I", 1:nrow(dt)),
                     mother=0, father=0)
    founderPop = importHaploSparse(haplo=haplo, 
                                   genMap=genMap,
                                   ploidy=2L)
    # founderPop = quickHaplo(nInd=Ne,nChr=1,segSites=nrow(dt), inbred = TRUE)
    SP = SimParam$new(founderPop)
    # 2) add the traits (columns from user) to take the values (rows) as marker effects
    for(iTrait in 1:length(traits)){
      SP$importTrait(markerNames =unlist(lapply(SP$genMap,names)), addEff = dt[,traits[iTrait]]) # over 2 because QTL data is diplodized or Q/2
    }
    # 3) set the population
    pop = newPop(founderPop, simParam = SP)
    if(nCrosses > 0){
      pop = randCross(pop, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)
      pop = makeDH(pop=pop, simParam = SP)
    }
    variances = diag(varG(pop))
    if(all(SP$varG>0)){
      pop = setPheno(pop,h2=rep(.98,length(which(variances>0))), simParam = SP, traits = which(variances > 0) )
    }else{
      pop@pheno <- apply(pop@pheno,2,function(xx){rnorm(length(xx))}) # rnorm(length(pop@pheno))
    }
  }else{
    pop=initPop
    SP=simParam
    variances = diag(varG(pop))
  }
  
  # ***) creating the frame for the plot
  indivPerformance <- list() # store results by generation
  averagePerformance <- matrix(0, nrow=nGenerations,ncol=4) # to store results
  colnames(averagePerformance) <- c("Average.fitness","Best.fitness","nQTL.mu", "deltaC.mu")
  # rownames(averagePerformance) <- paste("Generation",seq(nrow(averagePerformance)))
  # 4) Starting the Generational process
  ################################
  ################################
  ################################
  ## FOR EACH GENERATION
  j =0 # in 1:nGenerations
  nonStop=TRUE
  
  if(traceDelta){ # if user wants to trace inbreeding (default is TRUE but it can be a costly operation)
    m <- Matrix::Matrix(1,nrow=1,ncol=ncol(D))
    if(useD){
      mtDm <- as.vector((m%*%Matrix::tcrossprod(D,m))/(4*(ncol(D)^2)))
    }else{
      mtDm <- as.vector((Matrix::tcrossprod(m))/(4*(ncol(D)^2)))
    }
  }
  best <- list()
  console <-  data.frame(matrix(NA,nrow=nGenerations, ncol=8))
  nin <- nCrosses*nProgeny
  initVarG = round(sum(diag(varG(pop = pop))),3)
  while(nonStop) { # for each generation we breed # j=1
    j=j+1
    
    Q <- pullQtlGeno(pop, simParam = SP, trait = iTrait)/2 #?/2
    Q <- as(as(as( Q,  "dMatrix"), "generalMatrix"), "CsparseMatrix") # as(Q, Class = "dgCMatrix")
    rownames(Q) <- pop@id
    
    ## use mutatio rate
    if(mutRate > 0){
      if(nMutations == 1){
        pointMut = t(as.matrix(apply(data.frame(1:nInd(pop)), 1, function(x){
          sample(1:nrow(dt), nMutations, replace = FALSE)
        }) ))
      }else{
        pointMut = as.matrix(apply(data.frame(1:nInd(pop)), 1, function(x){
          sample(1:nrow(dt), nMutations, replace = FALSE)
        }) )
      }
      # 
      for(iQtl in unique(as.vector(pointMut))){
        modif=which(pointMut == iQtl, arr.ind = TRUE)[,"col"]
        allele = sample(0:1, 1)
        pop = editGenome(pop, ind=modif,chr=1, segSites=iQtl, simParam=SP, allele = allele)
      }
    }else{pointMut=as.data.frame(matrix(NA, nrow=0, ncol=1))}
    ## enf of use mutation rate
    ###########################
    ## if user wants to fix the number of QTLs activated apply the following rules
    # 1) if more than nQTLperInd we silence some
    # 2) if less than nQTLperInd we activate some
    if(fixQTLperInd){
      Qfq <- pullQtlGeno(pop, simParam = SP, trait = 1); Qfq <- Qfq/2
      for(iInd in 1:nInd(pop)){ # for each individual
        iQfq <- Qfq[iInd,]; areZeros <- which(iQfq == 0); areOnes <- setdiff(1:ncol(Qfq),areZeros)
        howMany <- sum(iQfq) # how many QTLs are activated, we're assuming is a 0/1 matrix
        toAddOrRem <- abs(howMany - nQTLperInd) # deviation from expectation
        if( howMany > nQTLperInd ){ # if exceeded silence some
          toRem <- sample(areOnes, toAddOrRem) # pick which ones will be silenced
          for(iChange in 1:toAddOrRem){
            pop = editGenome(pop, ind=iInd,chr=1, segSites=toRem[iChange], simParam=SP, allele = 0)
          }
        }else if( howMany < nQTLperInd){ # if lacked activate some
          toAdd <- sample(areZeros, toAddOrRem) # pick which ones will be activated
          for(iChange in 1:toAddOrRem){
            pop = editGenome(pop, ind=iInd,chr=1, segSites=toAdd[iChange], simParam=SP, allele = 1)
          }
        } # else do nothing
      }
    }
    ## enf of fixqtl
    
    a <- do.call(cbind, lapply(SP$traits, function(x){x@addEff}))
    colnames(a) <- traits
    pop@gv <- as.matrix(Q%*% a)
    fitnessValuePop<- do.call("fitnessf", args=list(Y=pop@gv, b=b,  Q=Q[pop@id,], 
                                                    a=a, D=D, lambda=lambda,
                                                    ... ), quote = TRUE)
    if(!is.matrix(fitnessValuePop)){
      fitnessValuePop <- Matrix::Matrix(fitnessValuePop,ncol=1) 
    }
    rownames(fitnessValuePop) <- pop@id
    nanFound <- which(is.nan( fitnessValuePop[,1] )) # check if there are nans
    if(length(nanFound) > 0){ # if so assign the worst values to those individuals/solutions
      if(selectTop){
        fitnessValuePop[nanFound,1]= -Inf
      }else{
        fitnessValuePop[nanFound,1]= Inf
      }
    }
    pop@pheno[,1] <- fitnessValuePop[,1]
    ## apply selection between and within
    structure = table(paste(pop@mother, pop@father))
    nc = length(structure)
    np = floor(median(structure))
    # Although multiple traits are enabled it is assumed that same QTLs are behind all the traits, differing only in their average allelic effects.
    if( propSelBetween[j] < 1){ 
      suppressWarnings( popF <- selectFam(pop=pop,nFam = ceiling(nc*propSelBetween[j]), trait = 1, 
                                          use = "pheno", simParam = SP, 
                                          selectTop=selectTop,... #H=H,nCities=nCities
      )@id, classes = "warning")
    }else{popF = pop@id}
    if( propSelWithin[j] < 1 ){
      suppressWarnings( popW <- selectWithinFam(pop = pop, nInd = ceiling(np*propSelWithin[j]), 
                                                trait = 1, use = "pheno", simParam = SP, 
                                                selectTop=selectTop,...#H=H,nCities=nCities
      )@id, classes = "warning")
    }else{popW=pop@id}
    ################################
    ################################
    ## FOR EACH TRAIT WE APPLY CONSTRAINTS
    constCheckUB <- constCheckLB <- matrix(1, nrow=nrow(Q), ncol=length(traits))
    
    for(iTrait in 1:length(traits)){ # iTrait=1
      # check the contraints and trace them back
      constCheckUB[,iTrait] <- ifelse( (pop@gv[,iTrait] > constraintsUB[iTrait]*(1+(1-constraintW[j])) )  , 0 , 1) # ifelse(c1+c2 < 2, 0, 1)
      constCheckLB[,iTrait] <- ifelse( (pop@gv[,iTrait] < constraintsLB[iTrait]*constraintW[j] ) , 0 , 1)
      # constCheckUB[,iTrait] <- ifelse( (pop@gv[,iTrait] > constraintsUB[iTrait])  , 0 , 1) # ifelse(c1+c2 < 2, 0, 1)
      # constCheckLB[,iTrait] <- ifelse( (pop@gv[,iTrait] < constraintsLB[iTrait]) , 0 , 1)
      nan0 <- which(is.nan( pop@gv[,iTrait]))
      if(length(nan0) > 0){ constCheckUB[nan0,iTrait] = 0; constCheckLB[nan0,iTrait] = 0 }
    } # end of for each trait
    # sum of how many trait constraints are met
    
    metConstCheck <- apply(constCheckUB,1,sum) # sum how many traits we're good to go
    didntMetConst <- which(metConstCheck < length(traits))
    
    # remove individuals that break the constraints UB
    if(length(didntMetConst)>0){ # 
      popCU <- pop@id[setdiff(1:nInd(pop),didntMetConst)]
    }else{popCU<-pop@id}
    #  remove individuals that break the constraints LB
    metConstCheckL <- apply(constCheckLB,1,sum) 
    didntMetConstL <- which(metConstCheckL < length(traits))
    
    # impute with mean value the ones that do not met the constraints
    if(length(didntMetConstL)>0){
      popCL <- pop@id[setdiff(1:nInd(pop),didntMetConstL)]
    }else{popCL<-pop@id}
    ## END OF FOR EACH TRAIT WE APPLY CONSTRAINTS
    ################################
    ################################
    parentsForSelection <- list(popF,popW,popCL, popCU)
    selected <- Reduce(intersect, parentsForSelection )
    if(length(selected) == 0){
      selected <- intersect(popCL,popCU )
    }
    if(length(selected) < 2){
      message("No legal solutions found. Selecting all individuals meeting constraints.")
      selected <- Reduce(intersect, list(popCL, popCU))
      if(length(selected) < 2){
        message("Too many constraints. No legal solutions found. Random selection p=0.5 applied.")
        selected <- pop@id[sample(1:nInd(pop), ceiling(nInd(pop)*.5) )]
      }
    }
    pop <- pop[which(pop@id %in% selected)]
    
    ## calculate trace metrics
    if(traceDelta){
      if(useD){
        qtDq <- Matrix::diag(Q%*%Matrix::tcrossprod(D,Q))
      }else{
        qtDq <- Matrix::diag(Matrix::tcrossprod(Q))
      } # [inbreedingNew (biggerVal) - inbreedingOld(smallerVal=outbred)] # qDq gets smaller and smaller with generations
      # deltaC represents that with every generation we get further away from orinal outbreeding and the lower the more inbred
      deltaC <- ( (qtDq/(4*(apply(Q/2,1,sum)^2))) - mtDm)/(1-mtDm) # numerator is equivalent to mtDm 
    }else{deltaC=NA}
    
    #################################
    # solutions selected for tracing
    best[[j]] <- selectInd(pop=pop, nInd = min(c(nInd(pop),topN)), trait = 1, 
                           use = "pheno", simParam = SP, 
                           selectTop=selectTop,... #H=H,nCities=nCities
    )
   
    mfvp =  mean(as.vector(fitnessValuePop[best[[j]]@id,]))
    indivPerformance[[j]] <- data.frame(id=best[[j]]@id, fitness=as.vector(fitnessValuePop[best[[j]]@id,]), 
                                        generation=j, nQTL=as.vector(apply(Q[best[[j]]@id,,drop=FALSE]/2,1,sum)),
                                        deltaC= as.vector(deltaC[best[[j]]@id]) ) # save individual solution performance
    
    averagePerformance[j,] <- c( mfvp , max(fitnessValuePop,na.rm=TRUE) ,  mean(apply(Q/2,1,sum),na.rm=TRUE), mean(deltaC,na.rm=TRUE) ) # save summaries of performance
    ## create new progeny
    for(k in 1:recombGens){
      if(nCrosses > 0){
        pop <- randCross(pop=pop, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)
      }
    }
    pop <- makeDH(pop=pop, nDH = 1, simParam = SP)
    
    #############################################
    ## compute phenotypes
    if(all(SP$varG>0)){
      pop = setPheno(pop,h2=rep(.98,length(which(variances>0))), simParam = SP, traits = which(variances > 0) )
    }else{
      pop@pheno <- apply(pop@pheno,2,function(xx){rnorm(length(xx))}) # rnorm(length(pop@pheno))
    }
    ##############################################################
    ##############################################################
    #store the performance of the jth generation for plot functions
    if(nrow(pop@gv) > 0){
      totalVarG = sum(diag(varG(pop = pop)))
    }else{
      totalVarG = 0
    }
    console[j,] <- c(j, length(didntMetConst), length(didntMetConstL), 
                     round(totalVarG, 3), round(mfvp, 3), round(propSelBetween[j], 2),
                     round(propSelWithin[j], 2),  Sys.time() )
    if(verbose){
      if(j==1){
        message(paste0("Population with ", nCrosses, " crosses, ", nProgeny, " progeny, and ",ncol(Q)," QTLs"))
        message(
          cat("gener  constUB  constLB  varG%   propB   propW   fit       time")
        )
      }
      sp <- paste(rep(" ", 2), collapse = "")
      message(cat(paste(
        sp, addZeros(1:nGenerations, nr=0)[j], 
        sp, addZeros(c(length(didntMetConst),nin), nr=0)[1],
        sp, addZeros(c(length(didntMetConstL),nin), nr=0)[1],
        sp, addZeros(c(round(totalVarG/initVarG, 3),1.000))[1],
        sp, addZeros(c( round(propSelBetween[j], 2) , round(propSelBetween, 2) ))[1], 
        sp, addZeros(c( round(propSelWithin[j], 2) , round(propSelWithin, 2) ))[1], 
        sp, addZeros(c(  round(mfvp, 3) , round(console[1:j,5], 3) ))[1], 
        sp, strsplit(format(Sys.time(), "%a %b %d %X %Y"), " ")[[1]][4]
      )))
    }
    if(j == nGenerations){nonStop = FALSE}
    if(nrow(pop@gv) > 0){
      if(totalVarG < tolVarG){nonStop = FALSE; message("Variance across traits exhausted. Early stop.")}
    }else{
      nonStop = FALSE; message("All individuals discarded. Consider changing some parameter values e.g., mutRate
                           or nQTLperInd (initial number of QTLs) to avoid all 
                           solutions to go beyond the bounds.")
    }
  }# end of for each generation
  ################################
  ################################
  ################################
  # 7) retrieve output of best combinations
  indivPerformance <- do.call(rbind, indivPerformance)
  best <- do.call(c, best)
  # transform the Pop 
  popEvola <- as(best,"evolaPop")
  popEvola@score <- averagePerformance[1:j,,drop=FALSE]
  popEvola@pointMut <- nrow(pointMut)
  popEvola@constCheckUB <- constCheckUB
  popEvola@constCheckLB <- constCheckLB
  popEvola@traits <- traits
  
  ###################
  # Although multiple traits are enabled it is assumed that same QTLs are behind all the traits, differing only in their average allelic effects.
  # we need to recalculate fitness because we have been scaling across generations
  Q <- pullQtlGeno(popEvola, simParam = SP, trait = iTrait)/2 #?/2
  Q <- as(as(as( Q,  "dMatrix"), "generalMatrix"), "CsparseMatrix") # as(Q, Class = "dgCMatrix")
  rownames(Q) <- popEvola@id
  a <- do.call(cbind, lapply(SP$traits, function(x){x@addEff}))
  popEvola@gv <- as.matrix(Q%*% a)
  fitnessValuePop<- do.call("fitnessf", args=list(Y=popEvola@gv, b=b,  Q=Q,
                                                  a=a, D=D, lambda=lambda,
                                                  ... ), quote = TRUE)

  if(!is.matrix(fitnessValuePop)){
    fitnessValuePop <- Matrix::Matrix(fitnessValuePop,ncol=1)
  };  rownames(fitnessValuePop) <- popEvola@id
  popEvola@fitness <- as.vector(fitnessValuePop)
  indivPerformance[,"fitness"] <-  as.vector(fitnessValuePop)
  popEvola@indivPerformance <- if(is.null(indivPerformance)){data.frame()}else{indivPerformance} 
  
  res <- list(pop=popEvola, simParam=SP, call=mc, fitness=fitnessValuePop, console=console)
  class(res) <- "evolaFitMod"
  return(res)
}

