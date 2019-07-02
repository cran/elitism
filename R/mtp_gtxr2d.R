# Biometrika 2014
# 2019-02-20
# Jiangtao Gou
# Example 1: mtp.gtxr2d(pvec.sorted=rev(c(0.97,0.83,0.79,0.73,0.61,0.59,0.47,0.37,0.31,0.29,0.23,0.19,0.13,0.07,0.03)), alpha=0.30,TRUE)
# Example 2: mtp.gtxr2d(pvec.sorted=c(0.1,0.28, 0.4), alpha=0.39,TRUE)
#
mtp.gtxr2d <- function (pvec.sorted, alpha, gc.is.included=FALSE) {
  if (gc.is.included) {
    pkev$global.count.IS <- 0
  }
  pvec.length <- length(pvec.sorted)
  #
  if (pvec.sorted[pvec.length] <= alpha) {
    #
    rej.idx <- pvec.length
    #
    return (list(rej.idx=rej.idx, init.count=1))
  } # End of if 
  #
  for (i in 2:pvec.length) {
    if (pvec.sorted[pvec.length-i+1] <= (i+1)/2/i*alpha) {
      #
      if (i == 2) {
        rej.idx <- InsertionSearch(a=pvec.sorted[1:(pvec.length-i+1)], value=alpha/i, low=1, high=pvec.length-i+1, gc.is.included)
      } else {
        rej.idx <- InsertionSearch(a=pvec.sorted[1:(pvec.length-i+1)], value=alpha/i*(1 + alpha^2/12*(1 - 1/(i-2)^2)), low=1, high=pvec.length-i+1, gc.is.included)
      }
      #
      return (list(rej.idx=rej.idx, init.count=i))
    } # End of if (pvec.sorted[pvec.length-i+1] <= (i+1)/2/i*alpha)
  }
  return (list(rej.idx=0, init.count=pvec.length))
}


