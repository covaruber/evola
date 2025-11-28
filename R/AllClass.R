# class for evola population class
setClass(
  "evolaPop",
  contains="Pop",
  slots=c(score="matrix", pointMut="numeric", indivPerformance="data.frame",
          constCheckUB="matrix", constCheckLB="matrix", traits="character",
          qtl="character",  pedTrack="data.frame",
          fitness="numeric", qtlData="data.frame")
) #-> evolaPop


