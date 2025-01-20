# class for evola population class
setClass(
  "evolaMod",
  contains="Pop",
  slots=c(Q="matrix", Qb="matrix", score="matrix", pointMut="numeric", indivPerformance="data.frame",
          constCheckUB="matrix", constCheckLB="matrix", traits="character", pedBest="data.frame")
) #-> evolaMod


