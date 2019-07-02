# Quick-00b
#  based on mtp.gtxr0
# 2019-03-01
# Jiangtao Gou
# Example: mtp.quick00b(pvec.sorted=rev(c(0.97,0.83,0.79,0.73,0.61,0.59,0.47,0.37,0.31,0.29,0.23,0.19,0.13,0.07,0.03)), alpha=0.30,TRUE)
#
mtp.quick00b <- function (pvec.sorted, alpha, gc.is.included=FALSE) {
  if (gc.is.included) {
    pkev$global.count.FS <- 0
    pkev$global.count.IS <- 0
  }
  #
  pvec.length <- length(pvec.sorted)
  #
  if (pvec.sorted[pvec.length] <= alpha) {
    #
    rej.idx <- pvec.length
    #
    return (list(rej.idx=rej.idx, init.count=1))
  } # End of if 
  #
  if (pvec.sorted[1] > alpha/2) {
    #
    rej.idx <- 0
    #
    return (list(rej.idx=rej.idx, init.count=1))
  } # End of if 
  #
  det.idx <- BinarySearch(a=pvec.sorted[1:(pvec.length-1)], value=alpha/2, low=1, high=pvec.length-1, gc.is.included, secondStage=FALSE)
  #
  rej.idx <- BinarySearch(a=pvec.sorted[1:det.idx], value=alpha/(pvec.length-det.idx+1), low=1, high=det.idx, gc.is.included, secondStage=TRUE)
  #
  return (list(rej.idx=rej.idx, init.count=pkev$global.count.FS+1))
}


