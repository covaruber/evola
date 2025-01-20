# class for evola population class
setClass(
  "evolaMod",
  contains="Pop",
  slots=c(Q="matrix", 
          score="matrix", pointMut="numeric", indivPerformance="data.frame",
          constCheckUB="matrix", constCheckLB="matrix", traits="character", pedTrack="data.frame")
) #-> evolaMod


