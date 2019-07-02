# critconst_quick
# 2019-02-20
# 2019-03-01, based on critconst_quick_v0.R
#
# Example: critconst_quick(alpha=0.05, n=10)
#
critconst_quick <- function (alpha=0.05, n=10) {
  # Objective function
  quick.equation <- function(c, alpha=0.05, n=10) {
    tIer <- alpha^n
    termA <- 0
    for (r in 1:(n-1)) {
      termA <- (1-alpha)^r*alpha^(n-r)*choose(n,r)
      termB <- 0
      for (j in (r+1):n) {
        termB <- termB + choose(n-r, j-r-1) *(1-c)^(j-r-1) *(c^(n-j+1) - (c-1/j)^(n-j+1))
      }
      tIer <- tIer + termA*termB
    }
    return(tIer-alpha)
  }
  #
  r <- uniroot(quick.equation, alpha=alpha, n=n, interval=c(1/2,3/4), tol = 1e-15)
  #
  return (r$root)
}
#


