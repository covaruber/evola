evolafit <- function(formula, dt, 
                     constraintsUB, constraintsLB, b,
                     nCrosses=50, nProgeny=20,nGenerations=20, 
                     recombGens=1, nChr=1, mutRate=0,
                     nQTLperInd=NULL, D=NULL, lambda=NULL,
                     propSelBetween=1,propSelWithin=0.5,
                     fitnessf=NULL, verbose=TRUE, dateWarning=TRUE,
                     selectTop=TRUE, tolVarG=1e-6, keepBest=FALSE, 
                     initPop=NULL, simParam = NULL, ...){
  
  my.date <- "2025-04-01"
  your.date <- Sys.Date()
  ## if your month is greater than my month you are outdated
  if(dateWarning & verbose){
    if (your.date > my.date) {
      warning("Version out of date. Please update evola to the newest version using:\ninstall.packages('evola') in a new session\n Use the 'dateWarning' argument to disable the warning message.", call. = FALSE)
    }
  }
  if(propSelBetween==0 | propSelWithin==0){stop("Please ensure that parameters propSelWithin and propSelBetween are different than zero.", call. = FALSE)}
  if(propSelBetween>0 & nCrosses==0){stop("If you apply selection between families you need to set nCrosses to a value > 0.", call. = FALSE)}
  if(missing(formula)){stop("Please provide the formula to know traits and classifiers.", call. = FALSE)}
  mc <- match.call() # create a call
  # add the fitness function (current options are qa=-1 to ln>=0 )
  if(is.null(fitnessf)){
    fitnessf <- ocsFun
  };
  
  elements <- strsplit(as.character(formula), split = "[+]")#[[1]]
  elements <- lapply(elements[-c(1)], function(x){all.vars(as.formula(paste("~",x)))})
  traits <- elements[[1]]
  if(!all(traits%in%colnames(dt))){stop("Specified traits are not traits in the dataset. Please correct.", call. = FALSE)}
  classifiers <- elements[[2]]
  
  if(is.null(initPop)){
    # check that the user has provided a single value for each QTL
    checkNQtls <- table(dt[,classifiers])
    if(length(which(checkNQtls > 1)) > 0){stop("You cannot provide more than one alpha value per QTL. Make sure that your x variable has only one value.", call. = FALSE)}
    if(missing(constraintsUB)){constraintsUB <- rep(Inf,length(traits))}
    if(length(constraintsUB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
    if(missing(constraintsLB)){constraintsLB <- rep(-Inf,length(traits))}
    if(length(constraintsLB) != length(traits)){stop(paste0("Constraints need to have the same length than traits (",length(traits),")"), call. = FALSE)}
    if(missing(b)){b <- rep(1,length(traits))}
    if(length(b) != length(traits)){stop(paste0("Weights need to have the same length than traits (",length(traits),")"), call. = FALSE)}
    if(is.null(lambda)){lambda <- 0}
    if(is.null(D)){D <- Matrix::Diagonal(nrow(dt))}
    if(is.null(nQTLperInd)){nQTLperInd <- nrow(dt)/5}
    nMutations = round(mutRate * nrow(dt)) # number of mutations per individual per generation
    # 1) initialize the population with customized haplotypes to ensure a single QTL per individual
    haplo = matrix(0, nrow=nrow(dt)*2, ncol = nrow(dt)) # rbind( diag(nrow(dt)), diag(nrow(dt)) )
    for (i in 1:nrow(haplo)) {
      haplo[i,sample(1:ncol(haplo), nQTLperInd )] <- 1
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
    founderPop = importHaplo(haplo=haplo, 
                             genMap=genMap,
                             ploidy=2L,
                             ped=ped)
    # founderPop = quickHaplo(nInd=Ne,nChr=1,segSites=nrow(dt), inbred = TRUE)
    SP = SimParam$new(founderPop)
    # 2) add the traits (columns from user) to take the values (rows) as marker effects
    for(iTrait in 1:length(traits)){
      SP$importTrait(markerNames =unlist(lapply(SP$genMap,names)), addEff = dt[,traits[iTrait]]/2) # over 2 because QTL data is diplodized
    }
    alpha = do.call(cbind,lapply(SP$traits, function(x){x@addEff}))
    # 3) set the population
    pop = newPop(founderPop, simParam = SP)
    if(nCrosses > 0){
      pop = randCross(pop, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)
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
  averagePerformance <- matrix(NA, nrow=nGenerations,ncol=5) # to store results
  colnames(averagePerformance) <- c("Average.qa","Best.qa","Average.qDq","nQTL.mu", "deltaC.mu")
  rownames(averagePerformance) <- paste("Generation",seq(nrow(averagePerformance)))
  # 4) Starting the Generational process
  ################################
  ################################
  ################################
  ## FOR GENERATION
  spacing9=paste(rep(" ",10), collapse = "");spacing99=paste(rep(" ",9), collapse = "");spacing999=paste(rep(" ",8), collapse = "")
  j =0 # in 1:nGenerations
  nonStop=TRUE
  pedBest <- list()
  best <- pop[0]; pedBest <- data.frame(matrix(NA,nrow=0, ncol=4)); colnames(pedBest) <- c("id","mother","father","gen")
  while(nonStop) { # for each generation we breed # j=1
    j=j+1
    if(j > 1){
      # group relationship
      qtDq <- Matrix::diag(Q%*%Matrix::tcrossprod(D,Q))
      # calculate base coancestry Ct
      m <- Matrix::Matrix(1,nrow=1,ncol=ncol(D))
      mtDm <- as.vector((m%*%Matrix::tcrossprod(D,m))/(4*(ncol(D)^2)))
      # rate of coancestry qtDq/4p^2 - mtAm/4n^2
      deltaC <- ( (qtDq/(4*(apply(Q/2,1,sum)^2))) - mtDm)/(1-mtDm)
      # if there is variation in min and max values in qtDq standardize
      if((max(qtDq)-min(qtDq)) > 0){ 
        qtDq = (qtDq-min(qtDq))/(max(qtDq)-min(qtDq)) # standardized xAx
      }
      ## apply selection between and within
      qtDq.lam = qtDq * lambda#[iTrait]
      names(qtDq.lam) <- pop@id
      structure = table(paste(pop@mother, pop@father))
      nc = length(structure)
      np = floor(median(structure))
      if( propSelBetween < 1){
        suppressWarnings( popF <- selectFam(pop=pop,nFam = round(nc*propSelBetween), trait = fitnessf, 
                                            b=b,d=qtDq.lam[pop@id],  Q=Q[pop@id,], use = "pheno", simParam = SP,
                                            selectTop=selectTop,...
        ), classes = "warning")
      }else{popF = pop}
      if( propSelWithin < 1 ){
        suppressWarnings( popW <- selectWithinFam(pop = pop, nInd = round(np*propSelWithin), 
                                                  trait = fitnessf,  b=b,d=qtDq.lam[pop@id],  Q=Q[pop@id,], use = "pheno", simParam = SP,
                                                  selectTop=selectTop, ...
        ), classes = "warning")
      }else{popW=pop}
      selected <- intersect(popF@id,popW@id)
      pop <- pop[which(pop@id %in% selected)]
      # solutions selected for tracing
      if(keepBest){
        best <- c(best, pop)
        pedBest = rbind(pedBest, data.frame(id=pop@id, mother=pop@mother, father=pop@father, gen=j) )
      }
      ## create new progeny
      for(k in 1:recombGens){
        if(nCrosses > 0){
          pop <- randCross(pop=pop, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)
        }
      }
      pop <- makeDH(pop=pop, nDH = 1, simParam = SP)
      ## compute constrained traits
      if(all(SP$varG>0)){
        pop = setPheno(pop,h2=rep(.98,length(which(variances>0))), simParam = SP, traits = which(variances > 0) )
      }else{
        pop@pheno <- apply(pop@pheno,2,function(xx){rnorm(length(xx))}) # rnorm(length(pop@pheno))
      }
    }
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
    # extract solutions for the trait 1 because all traits have the same QTLs
    Q <- pullQtlGeno(pop, simParam = SP, trait = 1)  
    Q <- as(Q, Class = "dgCMatrix")
    rownames(Q) <- pop@id
    constCheckUB <- constCheckLB <- matrix(1, nrow=nrow(Q), ncol=length(traits))
    qaFinal <- list()
    
    ################################
    ################################
    ## FOR EACH TRAIT
    for(iTrait in 1:length(traits)){ # iTrait=1
      qaOr <- Q %*% SP$traits[[iTrait]]@addEff # solutions * alpha for the iTrait
      qaFinal[[iTrait]] <- qaOr
      if((max(qaOr)-min(qaOr)) > 0){ # if there is variation
        qa = (qaOr-min(qaOr))/(max(qaOr)-min(qaOr)) # standardized qa
      }
      # calculate the genetic value of solutions using the objective functions
      pop@pheno[,iTrait] <- qa[,1] # as.vector( do.call(fitnessf[[iTrait]], list(dt=dt, qa=qa, qtDq.lam=qtDq.lam, pop=pop, Q=Q, alpha=alpha,...)) ) # qa - (lambda[iTrait] * qtDq ) # breeding value + coancestry
      
      # check the contraints and trace them back
      constCheckUB[,iTrait] <- ifelse( (qaOr[,1] > constraintsUB[iTrait])  , 0 , 1) # ifelse(c1+c2 < 2, 0, 1)
      constCheckLB[,iTrait] <- ifelse( (qaOr[,1] < constraintsLB[iTrait]) , 0 , 1)
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
      # message(paste(length(didntMetConst),"individuals discarded for breaking the lower bounds.", nInd(pop)-length(didntMetConst), "left." ))
      for(iTrait in 1:length(traits)){
        pop <- pop[setdiff(1:nInd(pop),didntMetConst)]
      }
    }
    #  remove individuals that break the constraints LB
    metConstCheckL <- apply(constCheckLB,1,sum) 
    didntMetConstL <- which(metConstCheckL < length(traits))
    # impute with mean value the ones that do not met the constraints
    if(length(didntMetConstL)>0){
      # message(paste(length(didntMetConst),"individuals discarded for breaking the lower bounds.", nInd(pop)-length(didntMetConst), "left." ))
      for(iTrait in 1:length(traits)){
        pop <- pop[setdiff(1:nInd(pop),didntMetConstL)]
      }
    }
    #store the performance of the jth generation for plot functions
    score <- do.call(cbind, qaFinal) %*% b
    if(j > 1){
      if(length(as.vector(score)) == length(as.vector(deltaC))){
        indivPerformance[[j]] <- data.frame(score=as.vector(score),deltaC= as.vector(deltaC) , qtDq= as.vector(qtDq), generation=j, nQTL=apply(Q/2,1,sum)) # save individual solution performance
      }
      averagePerformance[j,] <- c( mean(score,na.rm=TRUE), max(score,na.rm=TRUE) , mean(qtDq,na.rm=TRUE),  mean(apply(Q/2,1,sum),na.rm=TRUE), mean(deltaC,na.rm=TRUE) ) # save summaries of performance
    }
    if(j == nGenerations){nonStop = FALSE}
    
    if(nrow(pop@gv) > 0){
      totalVarG = sum(diag(varG(pop = pop)))
      if(totalVarG < tolVarG){nonStop = FALSE; message("Variance across traits exhausted. Early stop.")}
    }else{
      totalVarG = 0
      nonStop = FALSE; message("All individuals discarded. Consider changing some parameter values (e.g., mutRate).")
    }
    if(verbose){
      if(j==1){message("generation  constrainedUB  constrainedLB  total          varG")}
      nup=length(didntMetConst)
      nlb=length(didntMetConstL)
      nin=nInd(pop)
      message(paste("   ", j, ifelse(j <10,spacing9, spacing99), nup, ifelse(nup <10,spacing9, ifelse(nup <100,spacing99, spacing999)), 
                    nlb, ifelse(nlb <10,spacing9, ifelse(nlb <100,spacing99, spacing999)), 
                    nin, ifelse(nin <10,spacing9, ifelse(nin <100,spacing99, spacing999)), 
                    round(totalVarG,3)
      ))
    }
  }# end of for each generation
  ################################
  ################################
  ################################
  # 7) retrieve output of best combinations
  Q <- pullQtlGeno(pop, simParam = SP, trait=1); Q <- Q/2
  Qb <- pullQtlGeno(best, simParam = SP, trait=1); Qb <- Qb/2
  colnames(Q) <- colnames(Qb) <- apply(data.frame(dt[,classifiers]),1,function(x){paste(x,collapse = "_")})
  indivPerformance <- do.call(rbind, indivPerformance)
  # transform the Pop 
  popEvolaMod <- as(pop,"evolaMod")
  popEvolaMod@Q <- Q
  popEvolaMod@score <- averagePerformance[1:j,]
  popEvolaMod@pointMut <- nrow(pointMut)
  popEvolaMod@indivPerformance <- indivPerformance
  popEvolaMod@constCheckUB <- constCheckUB
  popEvolaMod@constCheckLB <- constCheckLB
  popEvolaMod@traits <- traits
  # transform the Pop
  bestEvolaMod <- as(best,"evolaMod")
  bestEvolaMod@Q <- Qb
  bestEvolaMod@pedTrack <- pedBest
  
  res <- list(pop=popEvolaMod, simParam=SP, call=mc, popBest=bestEvolaMod)
  return(res)
}

