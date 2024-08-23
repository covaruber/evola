pmonitor <- function(object,...){
  x <- object$score#[,"Best.xa"]
  x2 <- x[which(x[,"Average.xa"] < Inf & x[,"Average.xa"] > -Inf),"Average.xa"]
  x3 <- x[which(x[,"Best.xa"] < Inf & x[,"Best.xa"] > -Inf),"Best.xa"]
  mmin <- min(c(0, x2))
  mmax <- max(c( x3, x2 ) )
  plot(x[,1], type="o", ylim=c(mmin,mmax), xlab="Generation", ylab="Value")
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  if(ncol(x) > 1){
    for(i in 2:ncol(x)){
      par(new=TRUE)
      plot(x[,i], col=i, ylim=c(mmin,mmax), ylab="",xlab="", type="o",...)
    }
  }
  legend("topright",legend = colnames(x), col=1:(ncol(x)), bty="n", lty=1)
}

pareto <- function(object, scaled=TRUE, pch=20, xlim, ...){
  transp<-  function (col, alpha = 0.5) {
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, 
                                                  c[3]/255, alpha))
    return(res)
  }
  dt <- object$indivPerformance 
  dt2 <- as.data.frame(object$score)
  # prepare rate ot coancestry
  dt$deltaC <- dt$deltaC * -1
  # prepare performance
  if(scaled){
    minScore <- min(dt$score)
    maxScore <- max(dt$score)
    dt$score = (dt$score-minScore)/(maxScore-minScore) * 100 # standardized xa
  }else{dt$score <- dt$score/dt$nQTL }
  # prepare summaries of rate of coancestry
  dt2$deltaC.mu <- dt2$deltaC.mu  * -1 
  dt2 <- dt2[which(!is.nan(dt2$deltaC.mu)),]
  # prepare summaries of performance
  if(scaled){
    dt2$Average.xa <- (dt2$Average.xa-minScore)/(maxScore-minScore) * 100 # standardized xa
  }else{dt2$Average.xa <- dt2$Average.xa/dt2$nQTL.mu  }
  colfunc <- colorRampPalette(c("plum1", "plum4"))
  
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  
  layout(matrix(1:2,ncol=2), widths = c(2,1),heights = c(1,1))
  # left plot
  if(!scaled){ylabName="Maximum gain (units)"}else{ylabName="Maximum gain (%)"}
  dt$color <- transp(colfunc(max(dt$generation))[dt$generation], alpha = 0.4)
  
  if(missing(xlim)){xlim <- quantile(na.omit(dt$deltaC),c(.005,.995))}
  with(dt, plot(score~ deltaC, col=color, main="Pareto frontier", pch=pch,
                xlab="Rate of coancestry", ylab=ylabName, xlim=xlim, xaxt="n",  ... ))
  axis(1, at=seq(xlim[1],xlim[2],diff(xlim)/5),labels=round(seq(xlim[1]*-1,xlim[2]*-1, (diff(xlim)/5)*-1 ),3), col.axis="black")
  grid()
  lines(dt2$deltaC.mu, dt2$Average.xa, col = "blue")
  points(x=dt2$deltaC.mu[nrow(dt2)], y=dt2$Average.xa[nrow(dt2)], col="red", pch=20)
  # right plot
  legend_image <- as.raster(matrix(colfunc(max(dt$generation)), ncol=1))
  plot(c(0,2),c(0, max(dt$generation) ),type = 'n', axes = F,xlab = '', ylab = '', main = 'Generation')
  text(x=1.5, y = seq(min(dt$generation),max(dt$generation),l=5), labels = seq(max(dt$generation),min(dt$generation),l=5) )
  rasterImage(legend_image, 0, 0, 1, max(dt$generation) )
  # par(mfrow=c(1,1))
}