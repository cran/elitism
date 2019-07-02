# lowerConvexHull
# 2019-03-12
# lowerConvexHull(pvec=c(0, 0.01, 0.08, 0.10, 0.50, 0.70, 0.90))
# lowerConvexHull(pvec=c(0.01, 0.03, 0.04, 0.20))
# N<- 1e4; result.vec <- rep(0,times=N); for (i in 1:N) {n.null <-4; p.null <- runif(n=n.null, min=0, max=1); p.srt <- sort(p.null); result<- lowerConvexHull(pvec=p.srt); result.vec[i] <- result$counter}; print(max(result.vec))
#
lowerConvexHull <- function (pvec, n=length(pvec)) {
  pvec <- pvec[1:n]
  adn <- 2
  excludeTag <- -99
  epvec <- c(1,0,pvec)# expanded pvec
  eIvec <- rep(excludeTag,times=n+adn)
  eIvec[-1+adn] <- -1
  eIvec[0+adn] <- 0
  eIvec[1+adn] <- 1
  eIvec[2+adn] <- 2
  #
  k <- 2
  counter <- 0
  while (eIvec[k+adn] < n+1) {
    #print(k)
    #print(eIvec)
    counter <- counter + 1
    #
    forwardOrBackward <- (epvec[eIvec[k+adn] + adn] - epvec[eIvec[k-1+adn] + adn])/(eIvec[k+adn] - eIvec[k-1+adn]) > (epvec[eIvec[k-1+adn] + adn] - epvec[eIvec[k-2+adn] + adn])/(eIvec[k-1+adn] - eIvec[k-2+adn])
    #
    if (forwardOrBackward) {
      eIvec[k+1+adn] <- eIvec[k+adn] + 1
      k <- k+1
    } else {
      eIvec[k-1+adn] <- eIvec[k+adn]
      eIvec[k+adn] <- excludeTag
      k <- k-1
    } # End of if
  } # End of while
  eIvec[k+adn] <- excludeTag
  Ivec <- eIvec[which(eIvec>0)]
  #
  return (list(Ivec=Ivec,counter=counter))
}