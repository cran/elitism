# Holm
# 2019-02-22
# Jiangtao Gou
# Example 1: mtp.holm(pvec.sorted=rev(c(0.97,0.83,0.79,0.73,0.61,0.59,0.47,0.37,0.31,0.29,0.23,0.19,0.13,0.07,0.03)), alpha=0.80)
# Example 2: mtp.holm(pvec.sorted=c(0.01,0.02,0.03,0.04,0.22), alpha=0.2)
#
mtp.holm <- function (pvec.sorted, alpha) {
  pvec.length <- length(pvec.sorted)
  for (i in 1:pvec.length) {
    if (pvec.sorted[i] > alpha/(pvec.length-i+1)) {
      #
      rej.idx <- i-1
      #
      return (list(rej.idx=rej.idx, init.count=i))
    } # End of if 
  }
  return (list(rej.idx=pvec.length, init.count=pvec.length))
}


