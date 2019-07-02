# Hochberg
# 2019-02-18 
# Jiangtao Gou
# Example: mtp.hochberg(pvec.sorted=rev(c(0.97,0.83,0.79,0.73,0.61,0.59,0.47,0.37,0.31,0.29,0.23,0.19,0.13,0.07,0.03)), alpha=0.30)
#
mtp.hochberg <- function (pvec.sorted, alpha) {
  pvec.length <- length(pvec.sorted)
  for (i in 1:pvec.length) {
    if (pvec.sorted[pvec.length-i+1] <= alpha/i) {
      #
      rej.idx <- pvec.length-i+1
      #
      return (list(rej.idx=rej.idx, init.count=i))
    } # End of if 
  }
  return (list(rej.idx=0, init.count=pvec.length))
}


