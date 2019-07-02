# lowestAlphaLevelsB
# 2019-03-14
# 
# pvec=c(0, 0.01, 0.08, 0.10, 0.50, 0.70, 0.90); spvec= pvec; result1 <- lowerConvexHull(pvec=pvec); result2 <- lowestAlphaLevelsA(pvec=pvec, Ivec=result1$Ivec); lowestAlphaLevelsB(pvec=pvec, Ivec=result1$Ivec, Jvec=result2$Jvec, alpha=0.05);
#
lowestAlphaLevelsB <- function (pvec, Ivec, Jvec, alpha=0.05, n=length(pvec)) {
  m <- length(Ivec)
  counter <- 1
  #
  if (m == 1) {
    if (pvec[Ivec[m]] <= alpha) {
      return (list(orderP=Ivec[m],num=1,den=1,counter=1))
    } else {
      return (list(orderP=0,num=0,den=0,counter=1))
    }
  } # End of if m
  #
  if (pvec[Ivec[m]] <= alpha) {
    return (list(orderP=Ivec[m],num=1,den=1,counter=1))
  }
  #
  for (i in (m-1):1) {
    for (j in (Jvec[i+1]-1):Jvec[i]) {
      counter <- counter +1
      if (pvec[Ivec[i]] <= (Ivec[i] - j)/(Ivec[m] -j)*alpha) {
        return (list(orderP=Ivec[i], num=Ivec[i]-j, den=Ivec[m]-j, counter=counter))
      }
    } # End of for j
  } # End of for i
  return (list(orderP=0,num=0,den=0,counter=counter))
}