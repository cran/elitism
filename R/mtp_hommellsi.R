# mtp_hommell
# 2019-03-14
# Jiangtao Gou
#
# alpha = 0.10; pvec.sorted<- c(0.051,0.055,0.060, 0.120); pvec.sorted=c(0, 0.01, 0.08, 0.10, 0.50, 0.70, 0.90); mtp.hommelq(pvec.sorted=pvec.sorted, alpha=alpha); mtp.hommell(pvec.sorted=pvec.sorted, alpha=alpha)
#
# alpha = 0.05; mtp.hommell(pvec.sorted=spvec, alpha=alpha)
#
mtp.hommellsi <- function (pvec.sorted, alpha) {
  pvec.length <- length(pvec.sorted)
  #
  if (pvec.sorted[pvec.length] <= alpha) {
    #
    rej.idx <- pvec.length
    #
    return (list(rej.idx=rej.idx, init.count=1, scnd.count=0,consonance=1))
  }
  #
  LCHrslt <- lowerConvexHull(pvec=pvec.sorted)
  counter1 <- LCHrslt$counter
  #
  LALArslt <- lowestAlphaLevelsA(pvec=pvec.sorted, Ivec=LCHrslt$Ivec)
  counter2 <- LALArslt$counter
  #
  LALBrslt <- lowestAlphaLevelsB(pvec=pvec.sorted, Ivec=LCHrslt$Ivec, Jvec=LALArslt$Jvec, alpha=alpha);
  counter3 <- LALBrslt$counter
  #
  init.count <- counter1 + counter2 + counter3
  #
  i <- LALBrslt$den
  j <- LALBrslt$num
  #print(i)
  #print(j)
  #
  if (i == 0) {
    return (list(rej.idx=0, init.count=init.count, scnd.count=0, consonance=1))
  }
  #
  if (j == 1) {
    rej.idx <- pvec.length-i+1
    return (list(rej.idx=rej.idx, init.count=init.count, scnd.count=0,consonance=1))
  } # End of if (j == 1) 
  #
  pkev$global.count.IS <- 0
  #
  rej.idx <- InsertionSearch(a=pvec.sorted[1:(pvec.length-i+j-1)], value=alpha/(i-1), low=1, high=pvec.length-i+j-1, gc.is.included=TRUE, secondStage=TRUE)
  if (rej.idx > 0) {
    return (list(rej.idx=rej.idx, init.count=init.count, scnd.count=pkev$global.count.IS,consonance=1))
  }
#
  return (list(rej.idx=0, init.count=init.count, scnd.count=pvec.length - (i-j+2) + 1, consonance=0))
  #
}
#
#