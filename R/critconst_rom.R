# critconst_rom.R
# 2019-02-20/26
# Example 1: critconst_rom(n=2)
# Example 2: ptm.tot <- proc.time(); cc <- critconst_rom(alpha=0.05,n=1e3); myptm.tot <- proc.time() - ptm.tot; myptm.tot
#
critconst_rom <- function (alpha=0.05, n=10) {
  cRX <- rep(1,times=n)
  cRX[1:3] <- c(1, 1/2, 1/3*(1+alpha/4))
  if (n>=4) {
    for (i in 4:n) {
      partA <- (1-alpha^(i-1))/(1-alpha)
      partB <- 0
      for (k in 2:(i-1)) {
        partB <- partB + choose(n=i,k=k-1)*cRX[k]^(i-k+1) *alpha^(i-k)
      }
      cRX[i] <- (partA - partB)/i
    }
  } #if (n>=4) 
  return (cRX[1:n])
}
# 
