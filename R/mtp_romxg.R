# Rom-X-Global
# For large simulation
# 2019-02-26 10:55 AM
# Jiangtao Gou
# Example: critconst <- critconst_rom(alpha=0.80, n=15); mtp.romxg(pvec.sorted=rev(c(0.97,0.83,0.79,0.73,0.61,0.59,0.47,0.37,0.31,0.29,0.23,0.19,0.13,0.07,0.03)), alpha=0.80, critconst=critconst)
#
mtp.romxg <- function (pvec.sorted, alpha, critconst) {
  pvec.length <- length(pvec.sorted)
  cv <- alpha*critconst
#
  if (pvec.sorted[pvec.length] <= alpha) {
    #
    rej.idx <- pvec.length
    #
    return (list(rej.idx=rej.idx, init.count=1))
  } # End of if 
  for (i in 2:pvec.length) {
    if (pvec.sorted[pvec.length-i+1] <= cv[i]) {
      #
      rej.idx <- pvec.length-i+1
      #
      return (list(rej.idx=rej.idx, init.count=i))
    } # End of if 
  }
  return (list(rej.idx=0, init.count=pvec.length))
}


