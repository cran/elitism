# Biometrika 2014
# 2019-02-18 
# Jiangtao Gou
# Example: mtp.gtxr0b(pvec.sorted=seq(0.0001,0.1000,by=0.0001), alpha=0.05,TRUE)
#
mtp.gtxr0b <- function (pvec.sorted, alpha, gc.is.included=FALSE) {
  if (gc.is.included) {
    pkev$global.count.IS <- 0
  }
  #
  pvec.length <- length(pvec.sorted)
  #
  for (i in 1:pvec.length) {
    if (pvec.sorted[pvec.length-i+1] <= (i+1)/2/i*alpha) {
      #
      rej.idx <- BinarySearch(a=pvec.sorted[1:(pvec.length-i+1)], value=alpha/i, low=1, high=pvec.length-i+1, gc.is.included)
      #
      return (list(rej.idx=rej.idx, init.count=i))
    } # End of if (pvec.sorted[pvec.length-i+1] <= (i+1)/2/i*alpha)
  }
  return (list(rej.idx=0, init.count=pvec.length))
}


