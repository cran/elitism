# Biometrika 2014
# 2019-02-18 
# Jiangtao Gou
# Example: mtp.gtxr0b(pvec.sorted=rev(c(0.97,0.83,0.79,0.73,0.61,0.59,0.47,0.37,0.31,0.29,0.23,0.19,0.13,0.07,0.03)), alpha=0.30,TRUE)
#
mtp.gtxr0 <- function (pvec.sorted, alpha, gc.is.included=FALSE) {
  if (gc.is.included) {
    pkev$global.count.IS <- 0
  }
  #
  pvec.length <- length(pvec.sorted)
  #
  for (i in 1:pvec.length) {
    if (pvec.sorted[pvec.length-i+1] <= (i+1)/2/i*alpha) {
      #
      rej.idx <- InsertionSearch(a=pvec.sorted[1:(pvec.length-i+1)], value=alpha/i, low=1, high=pvec.length-i+1, gc.is.included)
      #
      return (list(rej.idx=rej.idx, init.count=i))
    } # End of if (pvec.sorted[pvec.length-i+1] <= (i+1)/2/i*alpha)
  }
  return (list(rej.idx=0, init.count=pvec.length))
}


