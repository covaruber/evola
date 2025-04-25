pmonitor <- function(object, kind=1, ...){
  x <- object$pop@score#[,"Best.qa"]
 
  if(kind==1){i=1;j=2}else if(kind==2){i=4;j=6}else{i=5; j=5}
  
  x2 <- x[which(x[,i] < Inf & x[,i] > -Inf),i]
  x3 <- x[which(x[,j] < Inf & x[,j] > -Inf),j]
  mmin <- min(c(0, x2))
  mmax <- max(c( x3, x2 ) )
  
  plot(x[,i], type="o", ylim=c(mmin,mmax), xlab="Generation", ylab="Value", col=i)
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  if(ncol(x) > 1){
    # for(i in 2){
      par(new=TRUE)
      plot(x[,j], col=j, ylim=c(mmin,mmax), ylab="",xlab="", type="o",...)
    # }
  }
  legend("topleft",legend = colnames(x)[c(i,j)], col=c(i,j), bty="n", lty=1)
}

pareto <- function(object, scaled=TRUE, pch=20, xlim, ...){
  transp<-  function (col, alpha = 0.5) {
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, 
                                                  c[3]/255, alpha))
    return(res)
  }
  dt <- object$pop@indivPerformance 
  if(length(which(!is.na(dt$deltaC))) == 0){
    stop("To see build this plot you need to run evolafit with argument traceDelta=TRUE", call. = FALSE)
  }
  # prepare rate ot coancestry
  dt$deltaC <- dt$deltaC * -1
  # prepare performance
  if(scaled){
    minFitness <- min(dt$fitness)
    maxFitness <- max(dt$fitness)
    dt$fitness = (dt$fitness-minFitness)/(maxFitness-minFitness) * 100 # standardized qa
  }else{dt$fitness <- dt$fitness }
  
  colfunc <- colorRampPalette(c("plum1", "plum4"))
  
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  
  layout(matrix(1:2,ncol=2), widths = c(2,1),heights = c(1,1))
  # left plot
  if(!scaled){ylabName="Maximum gain (units)"}else{ylabName="Maximum gain (%)"}
  dt$color <- transp(colfunc(max(dt$generation))[dt$generation], alpha = 0.4)
  
  if(missing(xlim)){
    # xlim <- quantile(na.omit(dt$deltaC),c(.005,.995))
    xlim=c(min(na.omit(dt$deltaC)),0)
  }
  with(dt, plot(fitness~ deltaC, col=color, main="Pareto frontier", pch=pch,
                xlab="Rate of coancestry", ylab=ylabName, xlim=xlim, xaxt="n",  ... ))
  axis(1, at=seq(xlim[1],xlim[2],diff(xlim)/5),labels=round(seq(xlim[1]*-1,xlim[2]*-1, (diff(xlim)/5)*-1 ),3), col.axis="black")
  grid()
  
  averages <- aggregate(cbind(deltaC,fitness)~generation, FUN=mean, data=dt)
  lines(averages$deltaC, averages$fitness, col = "blue")
  points(x=averages$deltaC[nrow(averages)], y=averages$fitness[nrow(averages)], col="red", pch=20)
  # right plot
  legend_image <- as.raster(matrix(colfunc(max(dt$generation)), ncol=1))
  plot(c(0,2),c(0, max(dt$generation) ),type = 'n', axes = F,xlab = '', ylab = '', main = 'Generation')
  text(x=1.5, y = seq(min(dt$generation),max(dt$generation),l=5), labels = seq(max(dt$generation),min(dt$generation),l=5) )
  rasterImage(legend_image, 0, 0, 1, max(dt$generation) )
  # par(mfrow=c(1,1))
}