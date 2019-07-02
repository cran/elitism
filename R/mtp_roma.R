# Rom-A
# 2019-03-07
# Jiangtao Gou
# Example: mtp.roma(pvec.sorted=rev(c(0.97,0.83,0.79,0.73,0.61,0.59,0.47,0.37,0.31,0.29,0.23,0.19,0.13,0.07,0.03)), alpha=0.30)
#
mtp.roma <- function (pvec.sorted, alpha) {
  pvec.length <- length(pvec.sorted)
  if (pvec.sorted[pvec.length] <= alpha) {
    #
    rej.idx <- pvec.length
    #
    return (list(rej.idx=rej.idx, init.count=1))
  } # End of if 
  if (pvec.sorted[pvec.length-1] <= alpha/2) {
    #
    rej.idx <- pvec.length-1
    #
    return (list(rej.idx=rej.idx, init.count=2))
  } # End of if 
  #
  for (i in 3:5) {
    if (pvec.sorted[pvec.length-i+1] <= alpha/i*(1 + (i-2)*alpha/2/(i-1))) {
      #
      rej.idx <- pvec.length-i+1
      #
      return (list(rej.idx=rej.idx, init.count=i))
    } # End of if 
  }
  #
  for (i in 6:pvec.length) {
    if (pvec.sorted[pvec.length-i+1] <= -log(1-alpha)/i) {
      #
      rej.idx <- pvec.length-i+1
      #
      return (list(rej.idx=rej.idx, init.count=i))
    } # End of if 
  }
  return (list(rej.idx=0, init.count=pvec.length))
}#
# 


