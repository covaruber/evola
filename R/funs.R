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
  indivPerformance <- list() # store results by generation
  averagePerformance <- matrix(NA, nrow=nGenerations,ncol=5) # to store results
  colnames(averagePerformance) <- c("Average.xa","Best.xa","Average.xAx","nQTL.mu", "deltaC.mu")
  rownames(averagePerformance) <- paste("Generation",seq(nrow(averagePerformance)))
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
    xtAx <- diag(Q%*%A%*%t(Q))
    m <- matrix(1,nrow=1,ncol=ncol(A))
    mAmt <- as.vector(m%*%A%*%t(m)/(4*(ncol(A)^2)))
    deltaC <- (xtAx - mAmt)/(1-mAmt)
    if((max(xtAx)-min(xtAx)) > 0){ # if there is variation
      xtAx = (xtAx-min(xtAx))/(max(xtAx)-min(xtAx)) # standardized xAx
    }
    for(iTrait in 1:length(traits)){ # iTrait=5
      xa <- Q %*% SP$traits[[iTrait]]@addEff
      if((max(xa)-min(xa)) > 0){ # if there is variation
        xa = (xa-min(xa))/(max(xa)-min(xa)) # standardized xa
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
    for(kk in 1:length(traits)){ # save merit of the different solutions
      xaFinal[[kk]] <- Q %*% SP$traits[[1]]@addEff
    }
    score <- do.call(cbind, xaFinal) %*% traitWeight
    indivPerformance[[j]] <- data.frame(score=score,deltaC=deltaC, xtAx=xtAx, generation=j) # save individual solution performance
    averagePerformance[j,] <- c( mean(score), max(score) , mean(diag(Q%*%A%*%t(Q))),  mean(apply(Q/2,1,sum)), mean(deltaC) ) # save summaries of performance
  }
  # 7) retrieve output of best combinations
  M <- pullQtlGeno(pop, simParam = SP, trait=1); M <- M/2
  colnames(M) <- apply(data.frame(dt[,classifiers]),1,function(x){paste(x,collapse = "_")})
  indivPerformance <- do.call(rbind, indivPerformance)
  return(list(M=M, score=averagePerformance, pheno=pop@pheno, pop=pop, indivPerformance=indivPerformance))
}


pmonitor <- function(object,...){
  x <- object$score#[,"Best.xa"]
  mmin <- min(c(0, x[,"Average.xa"]))
  mmax <- max(c( x[,"Best.xa"], x[,"Average.xAx"] ) )
  plot(x[,1], type="o", ylim=c(mmin,mmax), xlab="Generation", ylab="Value")
  if(ncol(x) > 1){
    for(i in 2:ncol(x)){
      par(new=TRUE)
      plot(x[,i], col=i, ylim=c(mmin,mmax), ylab="",xlab="", type="o",...)
    }
  }
  legend("topright",legend = colnames(x), col=1:(ncol(x)), bty="n", lty=1)
}

pareto <- function(object, scaled=TRUE, ...){
  transp<-  function (col, alpha = 0.5) {
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, 
                                                  c[3]/255, alpha))
    return(res)
  }
  dt <- object$indivPerformance 
  # prepare rate ot coancestry
  dt$deltaC <- dt$deltaC * -1
  # prepare performance
  if(scaled){
    minScore <- min(dt$score)
    maxScore <- max(dt$score)
    dt$score = (dt$score-minScore)/(maxScore-minScore) * 100 # standardized xa
  }
  # prepare summaries of rate of coancestry
  dt2 <- as.data.frame(object$score)
  dt2$deltaC.mu <- dt2$deltaC.mu  * -1 
  # prepare summaries of performance
  if(scaled){
  dt2$Average.xa <- (dt2$Average.xa-minScore)/(maxScore-minScore) * 100 # standardized xa
  }
  colfunc <- colorRampPalette(c("plum1", "plum4"))
  
  layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
  # left plot
  if(!scaled){ylabName="Maximum gain (units)"}else{ylabName="Maximum gain (%)"}
  dt$color <- transp(colfunc(max(dt$generation))[dt$generation], alpha = 0.4)
  with(dt, plot(score~ deltaC, col=color, main="Pareto frontier", pch=20,
                xlab="Rate of coancestry (%)", ylab=ylabName, xaxt="n", ... ))
  axis(1, at=seq(-100,0,10),labels=seq(100,0,-10), col.axis="black")
  grid()
  lines(dt2$deltaC.mu, dt2$Average.xa, col = "blue")
  # right plot
  legend_image <- as.raster(matrix(colfunc(max(dt$generation)), ncol=1))
  plot(c(0,2),c(0, max(dt$generation) ),type = 'n', axes = F,xlab = '', ylab = '', main = 'Generation')
  text(x=1.5, y = seq(min(dt$generation),max(dt$generation),l=5), labels = seq(max(dt$generation),min(dt$generation),l=5) )
  rasterImage(legend_image, 0, 0, 1, max(dt$generation) )
  par(mfrow=c(1,1))
  # library(ggplot2)
  # # Basic scatter plot
  # p <- ggplot(dt, aes(x=deltaC, y=score, col=generation)) 
  # p <- p + geom_point() + xlab("Rate of coancestry (dC)") + ylab("Gain")
  # for(i in 1:(nrow(dt2)-1)){
  #   if(i == (nrow(dt2)-1) ){
  #     p <- p + geom_segment(y = dt2$Average.xa[i], x = dt2$deltaC.mu[i], 
  #                           yend = dt2$Average.xa[i+1], xend = dt2$deltaC.mu[i+1],
  #                           arrow = arrow(length = unit(0.5, "cm"))
  #     )
  #   }else{
  #     p <- p + geom_segment(y = dt2$Average.xa[i], x = dt2$deltaC.mu[i], 
  #                           yend = dt2$Average.xa[i+1], xend = dt2$deltaC.mu[i+1] #,
  #                           #arrow = arrow(length = unit(0.5, "cm")) 
  #     )
  #   }
  #   
  # }
  # p
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