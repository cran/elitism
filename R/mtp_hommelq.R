# mtp_hommelq
# 2019-03-11
# Jiangtao Gou
# Example: mtp.hommelq(pvec.sorted=rev(c(0.97,0.83,0.79,0.73,0.61,0.59,0.47,0.37,0.31,0.29,0.23,0.19,0.13,0.07,0.03)), alpha=0.90)
#
# Example: mtp.hommelq(pvec.sorted=c(0.055,0.060, 0.120), alpha=0.10)
#
# Example: mtp.hommelq(pvec.sorted=c(0.045,0.060, 0.120), alpha=0.10)
#
# Example: mtp.hommelq(pvec.sorted=c(0.045,0.067, 0.120), alpha=0.10)
#
# Example: mtp.hommelq(pvec.sorted=c(0.051,0.055,0.060, 0.120), alpha=0.10)
#
# mtp.hommelq(pvec.sorted=spvec, alpha=0.05)
#
mtp.hommelq <- function (pvec.sorted, alpha) {
  pvec.length <- length(pvec.sorted)
  #
  if (pvec.sorted[pvec.length] <= alpha) {
    #
    rej.idx <- pvec.length
    #
    return (list(rej.idx=rej.idx, init.count=1, scnd.count=0,consonance=1))
  }
  #
  for (i in 2:pvec.length) {
    for (j in (i-1):1) {
    if (pvec.sorted[pvec.length-i+j] <= j*alpha/i) {
      if (j == 1) {
        rej.idx <- pvec.length-i+1
        return (list(rej.idx=rej.idx, init.count=i*(i-1)/2+1, scnd.count=0,consonance=1))
      } # End of if (j == 1) 
      #
      init.count <- i*(i-1)/2+2-j
      #
      for (k in (i-j+2):pvec.length) {
        if (pvec.sorted[pvec.length-k+1] <= alpha/(i-1)) {
          rej.idx <- pvec.length-k+1
          scnd.count <- k - (i-j+1)
          #
          #print(i)
          #print(j)
          #print(k)
          #
          return (list(rej.idx=rej.idx, init.count=init.count, scnd.count=scnd.count,consonance=1))
        }
      } # End of for k
      # break # Very important
      return (list(rej.idx=0, init.count=i*(i-1)/2+2-j, scnd.count=pvec.length - (i-j+2) + 1, consonance=0))
      #
    } # End of if 
    } # End of for j
  } # End of for i
  if (i == pvec.length && j == 1) {
    return (list(rej.idx=0, init.count=(pvec.length-1)*pvec.length/2+1, scnd.count=0, consonance=1))
  } else {
    return ('Error')
  }
}
#
#<https://stackoverflow.com/questions/37571028/breaking-out-of-nested-loops-in-r>
#

