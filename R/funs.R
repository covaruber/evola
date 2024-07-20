evola <- function(formula, dt, 
                 constraintsUB, constraintsLB, traitWeight,
                 nCrosses=50, nProgeny=40,nGenerations=30, recombGens=1,
                 nQTLperInd=1, A=NULL, lambda=1,
                 propSelBetween=1,propSelWithin=0.5){
  
  if(missing(formula)){stop("Please provide the formula to know traits and classifiers.", call. = FALSE)}
  elements <- strsplit(as.character(formula), split = "[+]")#[[1]]
  elements <- lapply(elements[-c(1)], function(x){all.vars(as.formula(paste("~",x)))})
  traits <- elements[[1]]
  classifiers <- elements[[2]]
  
  if(length(constraintsUB) != length(traits)){stop("Constraints need to have the same length than traits.", call. = FALSE)}
  if(length(constraintsLB) != length(traits)){stop("Constraints need to have the same length than traits.", call. = FALSE)}
  if(length(traitWeight) != length(traits)){stop("Weights need to have the same length than traits.", call. = FALSE)}
  if(length(lambda) != length(traits)){stop("Lambda need to have the same length than traits.", call. = FALSE)}
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
  # 3) set the population
  pop = newPop(founderPop, simParam = SP)
  pop = randCross(pop, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)
  pop = setPheno(pop,h2=.01, simParam = SP) # ignore h2 since we will replace it in line 56
  # ***) creating the frame for the plot
  Performance <- matrix(NA, nrow=nGenerations,ncol=4) # to store results
  colnames(Performance) <- c("Average.xa","Best.xa","Average.xAx","nQTL.mu")
  rownames(Performance) <- paste("Generation",seq(nrow(Performance)))
  # 4) Starting the Generational process
  for (j in 1:nGenerations) { # for each generation we breed # j=1
    print(j)
    if(j > 1){
      ## apply selection between and within
      pop <- selectFam(pop=pop,nFam = round(nCrosses*propSelBetween), trait = selIndex, b=traitWeight, use = "pheno", simParam = SP)
      pop <- selectWithinFam(pop = pop, nInd = round(nProgeny*propSelWithin), trait = selIndex,  b=traitWeight, use = "pheno", simParam = SP)
      ## create new progeny
      for(k in 1:recombGens){
        pop <- randCross(pop=pop, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)
      }
      pop <- makeDH(pop=pop, nDH = 1, simParam = SP)
      ## compute constrained traits
      pop <- setPheno(pop=pop, h2=.01, simParam = SP)  # ignore h2 since we will replace it in line 90
    }
    Q <- pullQtlGeno(pop, simParam = SP, trait = iTrait)  # plot((apply(Q,2,sum)/2)/nrow(Q))  # plot((apply(Q,1,sum)/2)/ncol(Q))
    for(iTrait in 1:length(traits)){ # iTrait=3
      cu <- Q %*% SP$traits[[iTrait]]@addEff
      if((max(cu)-min(cu)) > 0){
        cu = (cu-min(cu))/(max(cu)-min(cu)) # standardized cu
      }
      ctAc <- diag(Q%*%A%*%t(Q))
      if((max(ctAc)-min(ctAc)) > 0){
        ctAc = (ctAc-min(ctAc))/(max(ctAc)-min(ctAc)) # standardized cu
      }
      pop@pheno[,iTrait] <- cu - (lambda[iTrait] * ctAc ) # breeding value + coancestry
      # we only apply constraints to traits that will account for the total merit
      if(traitWeight[iTrait] > 0){ # if the trait will be used for selection
        ## pass each trait through all constraints
        for(iConUB in 1:length(constraintsUB)){ ## upper bound constraints ## for each trait constraint
          cuPrime <- Q %*% SP$traits[[iConUB]]@addEff
          if((max(cuPrime)-min(cuPrime)) > 0){
            ub = ( constraintsUB[iConUB]-min(cuPrime) )/(max(cuPrime)-min(cuPrime))
            cuPrime = (cuPrime-min(cuPrime))/(max(cuPrime)-min(cuPrime)) # standardized cu
          }else{
            ub=constraintsUB[iConUB]
          }
          constr <- cuPrime  #+ (lambda[iConUB] * ctAc )  # breeding value of the trait to be used for constraints
          pop@pheno[,iTrait] <- pop@pheno[,iTrait] * ifelse(constr[,1] > ub, 0 , 1) # adjust the BV of the trait after observing the constraints
        }
        for(iConLB in 1:length(constraintsLB)){ # lower bound constraint
          cuPrime <- Q %*% SP$traits[[iConLB]]@addEff
          if((max(cuPrime)-min(cuPrime))  >0 ){
            lb = (constraintsLB[iConLB]- -min(cuPrime) )/(max(cuPrime)-min(cuPrime)) 
            cuPrime = (cuPrime-min(cuPrime))/(max(cuPrime)-min(cuPrime)) # standardized cu
          }else{
            lb=constraintsLB[iConLB]
          }
          constr <- cuPrime # + (lambda[iConLB] * ctAc )  # breeding value of the trait to be used for constraints
          pop@pheno[,iTrait] <- pop@pheno[,iTrait] * ifelse(constr[,1] < lb, 0 , 1) # adjust the BV of the trait after observing the constraints
        } 
      }
    }
    #store the performance of the jth generation for plots
    xa <- list()
    for(kk in 1:length(traits)){xa[[kk]] <- Q %*% SP$traits[[1]]@addEff}
    score <- do.call(cbind, xa) %*% traitWeight
    Performance[j,] <- c( mean(score), max(score) , mean(diag(Q%*%A%*%t(Q))),  mean(apply(Q/2,1,sum)) )
  }
  # 7) retrieve output of best combinations
  M <- pullQtlGeno(pop, simParam = SP, trait=1); M <- M/2
  colnames(M) <- apply(data.frame(dt[,classifiers]),1,function(x){paste(x,collapse = "_")})
  return(list(M=M, score=Performance, pheno=pop@pheno))
}

# if(plotting){
#   plot(NULL, xlim =c(1,nGenerations), ylim = c(0, max((pop@pheno %*% traitWeight)[,1]) ), ylab = "Value", xlab = "Generation")
#   legend("bottomright", c("Best x'a", "Average x'a", "Average x'Ax", "Mean nQTL"), pch = 20, col = c(1, 2, 3, 4),bty = "n")
# }
# if(plotting){
#   # update plot 
#   if(j %in% seq(1,nGenerations,1)) {
#     Sys.sleep(0.1)
#     lines(Performance[,2])
#     lines(Performance[,1], col = "red")
#     lines(Performance[,3], col = "green")
#     lines(Performance[,4], col = "blue")
#   }
#   Sys.sleep(0)
# } # end of plotting
