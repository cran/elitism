# Quick-01b
#  based on mtp.gtxr0
# 2019-03-01
# Jiangtao Gou
# Example: mtp.quick01b(pvec.sorted=c(0.01, 0.06,0.3, 0.7), alpha=0.61,TRUE)
#
mtp.quick01b <- function (pvec.sorted, alpha, gc.is.included=FALSE) {
  if (gc.is.included) {
    pkev$global.count.FS <- 0
    pkev$global.count.IS <- 0
  }
  #
  pvec.length <- length(pvec.sorted)
  #
  ca <- pvec.length/2/(pvec.length-1)*alpha
  #
  if (pvec.sorted[pvec.length] <= alpha) {
    #
    rej.idx <- pvec.length
    #
    return (list(rej.idx=rej.idx, init.count=1))
  } # End of if 
  #
  if (pvec.sorted[1] > ca) {
    #
    rej.idx <- 0
    #
    return (list(rej.idx=rej.idx, init.count=1))
  } # End of if 
  #
  det.idx <- BinarySearch(a=pvec.sorted[1:(pvec.length-1)], value=ca, low=1, high=pvec.length-1, gc.is.included, secondStage=FALSE)
  #
  #print(det.idx)
  #
  rej.idx <- BinarySearch(a=pvec.sorted[1:det.idx], value=alpha/(pvec.length-det.idx+1), low=1, high=det.idx, gc.is.included, secondStage=TRUE)
  #
  return (list(rej.idx=rej.idx, init.count=pkev$global.count.FS+1))
}


