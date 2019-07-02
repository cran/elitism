# strange.pvec <- check.hommelq.hommell(); p.adjust(strange.pvec, method='hommel')
#
check.hommelq.hommell <- function (n=20, N=1000, alpha=0.05, rho=0.5, null.pcnt=0.5, altr.prmtr=2, mtp.type='hommelq', mtp.type.rf='hommell', seed=201903, gc.is.included=TRUE, filename='writeCAT.txt') {
  # Hommel Tag
  hommel.tag <- FALSE
  if (mtp.type=='hommelq' || mtp.type=='hommell') {
    hommel.tag <- TRUE
  }
  # Check rho
  if (rho <= -1/(n-1)) {
    rho <- -1/n
    print('The common correlation should not be less than -1/(n-1).')
  }
  # Start timer
  ptm.tot <- proc.time()
  # Check method
  if (mtp.type == 'romxg') {
    critconst <- critconst_rom(alpha=alpha, n=n)
  }
  #
  if (mtp.type == 'quickxg' || mtp.type == 'quickxbg') {
    critconst <- critconst_quick(alpha=alpha, n=n)
  }
  #
  set.seed(seed=seed)
  #
  n.null <- round(n*null.pcnt)
  n.altr <- n - n.null
  #
  if (rho != 0) {
    mu <- c(rep(altr.prmtr,times=n.altr), rep(0,times=n.null)) 
    Sigma <- (1-rho)*diag(n) + rho*matrix(rep(1,times=n*n),nrow=n)
  }
  #
  indctr.null <- c(rep(0,times=n.altr),rep(1,times=n.null))
  rej.rslt <- rep(0,times=n) # arbitrary
  #
  if (gc.is.included) {
    IS.count.RSLT <- rep(0,times=N)
    pkev$global.count.IS <- 0
  }
  #
  init.count.RSLT <- rep(0,times=N)
  tpos <- rep(0,times=N)
  fpos <- rep(0,times=N)
  tneg <- rep(0,times=N)
  fneg <- rep(0,times=N)
  if (hommel.tag) {
    cnsnnc <- rep(0,times=N)
  }
  #
  for (iN in 1:N) {
    #print('iN'); print(iN)
    #
    if (rho == 0) {
      p.null <- runif(n=n.null, min=0, max=1)
      p.altr <- pnorm(q=rnorm(n=n.altr,mean=altr.prmtr, sd=1), mean=0, sd=1, lower.tail=FALSE)
      p.vec <- c(p.altr, p.null) 
    } else {
      nrn <- MASS::mvrnorm(n=1, mu=mu, Sigma=Sigma, tol = 1e-6)
      p.vec <- pnorm(q=nrn, mean=0, sd=1, lower.tail=FALSE)
    } # End if (rho == 0)
    #
    p.srt <- sort(p.vec, method='quick', index.return=TRUE)
    #
    #if( iN == 438) {
    #  myvec <- p.srt
    #  print(p.srt)
    #}
    ###
    # Select a method
    ###
    mtp.rslt <- switch(mtp.type,
                       # Quick (10)
                       quick00 = mtp.quick00(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick01 = mtp.quick01(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick10 = mtp.quick10(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick11 = mtp.quick11(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quickxg = mtp.quickxg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst, gc.is.included=gc.is.included),
                       quick00b = mtp.quick00b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick01b = mtp.quick01b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick10b = mtp.quick10b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick11b = mtp.quick11b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quickxbg = mtp.quickxbg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst, gc.is.included=gc.is.included),
                       # GTXR (6)
                       gtxr0 = mtp.gtxr0(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr1c = mtp.gtxr1c(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr2d = mtp.gtxr2d(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr0b = mtp.gtxr0b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr1cb = mtp.gtxr1cb(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr2db = mtp.gtxr2db(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       # Rom (3=4-1)
                       rom1 = mtp.rom1(pvec.sorted=p.srt$x, alpha=alpha), 
                       romx = mtp.romx(pvec.sorted=p.srt$x, alpha=alpha), 
                       romxg = mtp.romxg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst),
                       roma = mtp.roma(pvec.sorted=p.srt$x, alpha=alpha),
                       # Hommel, Hochberg and Holm (4)
                       hommelq = mtp.hommelq(pvec.sorted=p.srt$x, alpha=alpha), 
                       hommell = mtp.hommell(pvec.sorted=p.srt$x, alpha=alpha), 
                       hochberg = mtp.hochberg(pvec.sorted=p.srt$x, alpha=alpha), 
                       holm = mtp.holm(pvec.sorted=p.srt$x, alpha=alpha)
    ) # End of switch
    #
    #print(mtp.rslt$rej.idx)
    #print(p.srt$x)
    #
    rej.rslt.srt <- c(rep(1,times=mtp.rslt$rej.idx), rep(0,n-mtp.rslt$rej.idx))
    #
    rej.rslt[p.srt$ix] <- rej.rslt.srt
    rej.rslt.srt.1 <- rej.rslt.srt
    rej.rslt.1 <- rej.rslt
    #
    init.count.RSLT[iN] <- mtp.rslt$init.count
    fpos[iN] <- sum(rej.rslt*indctr.null)
    tneg[iN] <- n.null - fpos[iN]
    tpos[iN] <- sum(rej.rslt*(1-indctr.null))
    fneg[iN] <- n.altr - tpos[iN]
    if (gc.is.included) {
      IS.count.RSLT[iN] <- pkev$global.count.IS
      pkev$global.count.IS <- 0
    }
    if (hommel.tag) {
      IS.count.RSLT[iN] <- mtp.rslt$scnd.count
      cnsnnc[iN] <- mtp.rslt$consonance
    }
    #
    #
    #
    ###
    # Select another method
    ###
    mtp.rslt <- switch(mtp.type.rf,
                       # Quick (10)
                       quick00 = mtp.quick00(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick01 = mtp.quick01(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick10 = mtp.quick10(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick11 = mtp.quick11(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quickxg = mtp.quickxg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst, gc.is.included=gc.is.included),
                       quick00b = mtp.quick00b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick01b = mtp.quick01b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick10b = mtp.quick10b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quick11b = mtp.quick11b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       quickxbg = mtp.quickxbg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst, gc.is.included=gc.is.included),
                       # GTXR (6)
                       gtxr0 = mtp.gtxr0(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr1c = mtp.gtxr1c(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr2d = mtp.gtxr2d(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr0b = mtp.gtxr0b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr1cb = mtp.gtxr1cb(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       gtxr2db = mtp.gtxr2db(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included), 
                       # Rom (3=4-1)
                       rom1 = mtp.rom1(pvec.sorted=p.srt$x, alpha=alpha), 
                       romx = mtp.romx(pvec.sorted=p.srt$x, alpha=alpha), 
                       romxg = mtp.romxg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst),
                       roma = mtp.roma(pvec.sorted=p.srt$x, alpha=alpha),
                       # Hommel, Hochberg and Holm (4)
                       hommelq = mtp.hommelq(pvec.sorted=p.srt$x, alpha=alpha), 
                       hommell = mtp.hommell(pvec.sorted=p.srt$x, alpha=alpha), 
                       hochberg = mtp.hochberg(pvec.sorted=p.srt$x, alpha=alpha), 
                       holm = mtp.holm(pvec.sorted=p.srt$x, alpha=alpha)
    ) # End of switch
    #
    #print(mtp.rslt$rej.idx)
    #print(p.srt$x)
    #
    rej.rslt.srt <- c(rep(1,times=mtp.rslt$rej.idx), rep(0,n-mtp.rslt$rej.idx))
    #
    rej.rslt[p.srt$ix] <- rej.rslt.srt
    rej.rslt.srt.2 <- rej.rslt.srt
    rej.rslt.2 <- rej.rslt
    #
    init.count.RSLT[iN] <- mtp.rslt$init.count
    fpos[iN] <- sum(rej.rslt*indctr.null)
    tneg[iN] <- n.null - fpos[iN]
    tpos[iN] <- sum(rej.rslt*(1-indctr.null))
    fneg[iN] <- n.altr - tpos[iN]
    if (gc.is.included) {
      IS.count.RSLT[iN] <- pkev$global.count.IS
      pkev$global.count.IS <- 0
    }
    if (hommel.tag) {
      IS.count.RSLT[iN] <- mtp.rslt$scnd.count
      cnsnnc[iN] <- mtp.rslt$consonance
    }
    #
    #
    #print(iN)
    #if (iN == 65) {
    #  print(iN)
    #  print(p.srt$x)
    #  print(rej.rslt.srt.1)
    #  print(rej.rslt.srt.2)
    #}
    if (any(rej.rslt.1 != rej.rslt.2)) {
      print('rej.iN'); print(iN)
      #print(p.srt$x)
      #print(rej.rslt.srt.1)
      #print(rej.rslt.srt.2)
      #print(p.adjust(p.srt$x,method='hommel'))
      #strange.pvec <- p.srt$x
    }
  } # End of for (iN in 1:N)
  #return (strange.pvec)
  return (0)
  #
}
###
# c(0.5004383, 0.5546876, 0.6801005, 0.9940105, 0.9994619)
