# lowestAlphaLevels
# 2019-03-13
# 
# pvec <- c(0.01, 0.03, 0.04, 0.20); result <- lowerConvexHull(pvec=pvec); lowestAlphaLevelsA(pvec=pvec, Ivec=result$Ivec)
#
# pvec=c(0, 0.01, 0.08, 0.10, 0.50, 0.70, 0.90); result <- lowerConvexHull(pvec=pvec); lowestAlphaLevelsA(pvec=pvec, Ivec=result$Ivec)
#
# pvec=c(0.5004383, 0.5546876, 0.6801005, 0.9940105, 0.9994619); result <- lowerConvexHull(pvec=pvec); lowestAlphaLevelsA(pvec=pvec, Ivec=result$Ivec)
#
lowestAlphaLevelsA <- function (pvec, Ivec, n=length(pvec)) {
  m <- length(Ivec)
  if (m == 1) {
    return(list(Jvec=0,counter=0))
  }
  initVec <- c(0, Ivec[1:(m-1)]-1)
  Jvec <- initVec
  #print(Jvec)
  counter <- 0
  #
  for (i in m:2) {
    k <- 0
    slopeA <- pvec[Ivec[i]]/(Ivec[i] - initVec[i] - k)
    slopeB <- pvec[Ivec[i-1]]/(Ivec[i-1]  - initVec[i] - k)
    #print('i,k,A,B')
    #print(i)
    #print(k)
    #print(slopeA)
    #print(slopeB)
    moveInd <- slopeA < slopeB
    #
    counter <- counter+1
    #
    while (moveInd) {
      k <- k-1
      slopeA <- pvec[Ivec[i]]/(Ivec[i] - initVec[i] - k)
      slopeB <- pvec[Ivec[i-1]]/(Ivec[i-1]  - initVec[i] - k)
      #print('i,k,A,B')
      #print(i)
      #print(k)
      #print(slopeA)
      #print(slopeB)
      moveInd <- slopeA < slopeB #
      counter <- counter+1
      #print(Jvec)
    } # End of while
    k <- k+1
    Jvec[i] <- initVec[i] + k
    # Update the initial search point to make sure the time complexity is linear.
    initVec[i-1] <- min(initVec[i-1], Jvec[i]-1)
    #
  }
  return(list(Jvec=Jvec,counter=counter))
}
#
#
#
# [1] 0.0002425157 0.0011716709 0.0031837686 0.0045148425 0.0064246729
# [6] 0.0074960615 0.0091432874 0.0416449245 0.0529986855 0.1801078617
# [11] 0.3423716407 0.4163590073 0.4524396561 0.5104196527 0.5198045983
# [16] 0.5259765805 0.5668694328 0.6089291836 0.7007743240 0.7334719171
#
