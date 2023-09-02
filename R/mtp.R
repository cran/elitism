#' @name mtp
#' @title Multiple Test Procedures
#'
#' @description Given a set of p-values and the level of significance, returns decisions whether the corresponding hypotheses should be rejected or not, including the hybrid Hochberg-Hommel procedure (Gou et al., 2014) and Quick procedure (Gou and Zhang, 2022).
#'
#' @param p vector of p-values.
#' @param alpha the level of significance.
#' @param method multiplicity correction method, including the Holm procedure ("holm"), the Hochberg procedure ("hochberg", "chochberg"), the Hommel procedure ("hommel", "hommelq", "hommell", "hommellsi",  "hommellsb"), the Rom procedure ("rom", "rom1", "roma", "romx"), the Gou-Tamhane-Xi-Rom procedure ("gtxr", "gtxr0i", "gtxr1ci", "gtxr2di", "gtxr0b", "gtxr1cb", "gtxr2db"), and the Quick procedure ("quick", "quick00i", "quick01i", "quick10i", "quick11i", "quickxi", "quick00b", "quick01b", "quick10b", "quick11b", "quickxb").
#' @param n number of p-values.
#' @return a list, including a binary vector of rejections, the total number of comparisons, and an indicator of consonance.
#' @author Jiangtao Gou
#' @references
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics} \bold{6}, 65-70.
#'
#' Hochberg, Y. and Tamhane, A. C. (1987). Multiple Comparison Procedures. John Wiley and Sons, New York.
#'
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. \emph{Biometrika} \bold{75}, 800-802.
#'
#' Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. \emph{Biometrika} \bold{75}, 383-386.
#'
#' Rom, D. M. (1990). A sequentially rejective test procedure based on a modified Bonferroni inequality. \emph{Biometrika} \bold{77}, 663-665.
#'
#' Wright, S. P. (1992). Adjusted p-values for simultaneous inference. \emph{Biometrics} \bold{48}, 1005-1013.
#'
#' Gou, J., Tamhane, A. C., Xi, D., and Rom, D. (2014). A class of improved hybrid Hochberg-Hommel type step-up multiple test procedures. \emph{Biometrika} \bold{101}, 899-911.
#'
#' Gou, J., and Tamhane, A. C. (2014). On generalized Simes critical constants. \emph{Biometrical Journal} \bold{56}, 1035-1054.
#'
#' Gou, J., and Tamhane, A. C. (2018). Hochberg procedure under negative dependence. \emph{Statistica Sinica} \bold{28}, 339-362.
#'
#' Tamhane, A. C., and Gou, J. (2018). Advances in p-value based multiple test procedures. \emph{Journal of Biopharmaceutical Statistics} \bold{28}, 10-27.
#'
#' Meijer, R. J., Krebs, T. J. P., and Goeman, J. J. (2019). Hommel's procedure in linear time. \emph{Biometrical Journal} \bold{61}, 73-82.
#'
#' Tamhane, A. C., and Gou, J. (2022). Chapter 2 Multiple test procedures based on p-values. In X. Cui, T. Dickhaus, Y. Ding, and J. C. Hsu (Eds.), \emph{Handbook of multiple comparisons} (Vol. 45, pp. 11-34).
#'
#' Gou, J.(2022). Quick multiple test procedures and p-value adjustments, \emph{Statistics in Biopharmaceutical Research} \bold{14}, 636-650.
#'
#' @details
#' Given a set of p-values, returns a binary vector of decisions, where 1 stands for rejection, and 0 stands for acceptance. There are six families of procedures.
#' \enumerate{
#'  \item Holm procedure (1 procedure)
#'  \enumerate{
#'   \item \emph{holm}, the Holm (1979) step-down method.
#'   }
#'  \item Hochberg procedure (2 procedures)
#'  \enumerate{
#'   \item \emph{hochberg}, the Hochberg (1988) step-up method.
#'   \item \emph{chochberg}, the conservative Hochberg method developed by Gou and Tamhane (2018).
#'   }
#'  \item Hommel procedure (5 procedures)
#'  \enumerate{
#'   \item \emph{hommel}, the Hommel (1988) step-up method, linear time algorithm with standard binary search, equivalent to \emph{hommellsb}.
#'   \item \emph{hommelq}, the Hommel (1988) step-up method, quadratic time algorithm.
#'   \item \emph{hommell}, the Hommel (1988) step-up method, linear time algorithm by Meijer, Krebs and Goeman (2019).
#'   \item \emph{hommellsb}, the Hommel (1988) step-up method, linear time algorithm by Meijer, Krebs and Goeman (2019), with standard binary search enhancement.
#'   \item \emph{hommellsi}, the Hommel (1988) step-up method, linear time algorithm by Meijer, Krebs and Goeman (2019), with interpolation search enhancement.
#'   }
#'   \item Rom procedure (4 procedure)
#'  \enumerate{
#'   \item \emph{rom}, the Rom (1990) step-up method, equivalent to \emph{romx}.
#'   \item \emph{rom1}, the Rom-1 method proposed by Gou and Zhang (2020).
#'   \item \emph{roma}, the Rom-1A method proposed by Gou and Zhang (2020).
#'   \item \emph{romx}, the Rom (1990) step-up method.
#'   }
#'  \item Gou-Tamhane-Xi-Rom procedure (7 procedure)
#'  \enumerate{
#'   \item \emph{gtxr}, the zeroth order hybrid Hommel-Hochberg procedure, proposed by Gou et al. (2014), equivalent to \emph{gtxr0b}.
#'   \item \emph{gtxr0b}, the zeroth order GTXR procedure, with standard binary search enhancement.
#'   \item \emph{gtxr1cb}, the GTXR procedure with refined c critical constants, with standard binary search enhancement.
#'   \item \emph{gtxr2db}, the GTXR procedure with refined d critical constants, with standard binary search enhancement.
#'   \item \emph{gtxr0i}, the zeroth order GTXR procedure, with interpolation search enhancement.
#'   \item \emph{gtxr1ci}, the GTXR procedure with refined c critical constants, with interpolation search enhancement.
#'   \item \emph{gtxr2di}, the GTXR procedure with refined d critical constants, with interpolation search enhancement.
#'   }
#'  \item Quick procedure (11 procedure)
#'  \enumerate{
#'   \item \emph{quick}, the Quick method, proposed by Gou and Zhang (2020), equivalent to \emph{quick00b}.
#'   \item \emph{quick00b}, the zeroth order Quick procedure, proposed by Gou and Zhang (2020), with standard binary search enhancement.
#'   \item \emph{quick01b}, the Quick procedure with refined d critical constants, with standard binary search enhancement.
#'   \item \emph{quick10b}, the Quick procedure with refined c critical constants, with standard binary search enhancement.
#'   \item \emph{quick11b}, the Quick procedure with refined c and d critical constants, with standard binary search enhancement.
#'   \item \emph{quickxb}, the exact Quick procedure with refined c critical constants, with standard binary search enhancement.
#'   \item \emph{quick00i}, the zeroth order Quick procedure, proposed by Gou and Zhang (2020), with interpolation search enhancement.
#'   \item \emph{quick01i}, the Quick procedure with refined d critical constants, with interpolation search enhancement.
#'   \item \emph{quick10i}, the Quick procedure with refined c critical constants, with interpolation search enhancement.
#'   \item \emph{quick11i}, the Quick procedure with refined c and d critical constants, with interpolation search enhancement.
#'   \item \emph{quickxi}, the exact Quick procedure with refined c critical constants, with interpolation search enhancement.
#'   }
#' }
#' @examples
#' library(elitism)
#'  pvalues.raw <- c(0.002,0.007,0.005,0.024,0.022,0.009,0.007,0.036,0.060,0.035)
#'  pkev <- new.env(); pkev$global.count.IS <- 0; pkev$global.count.FS <- 0;
#'  decision.hoch <- mtp(pvalues.raw, alpha = 0.025, method = "hochberg")
#'  pkev <- new.env(); pkev$global.count.IS <- 0; pkev$global.count.FS <- 0;
#'  decision.quick <- mtp(pvalues.raw, alpha = 0.025, method = "quick")
#'  pkev <- new.env(); pkev$global.count.IS <- 0; pkev$global.count.FS <- 0;
#'  decision.gtxr <- mtp(pvalues.raw, alpha = 0.025, method = "gtxr")
#' @seealso \code{elitism::p.adjust}
#' @export
#' @import stats
#'
# utils::globalVariables(c("global.count.IS", "global.count.FS"))
#
#global.count.IS <<- 0
#global.count.FS <<- 0
#gc.is.included <<- TRUE
#
#pkev <- new.env()
#
mtp <- function (p, alpha = 0.05, method = "gtxr", n = length(p)) {
  # pkev <- new.env()
  #pkev$global.count.IS <- 0
  #pkev$global.count.FS <- 0
  #
  p.vec <- p
  mtp.type <- method
  gc.is.included <- TRUE
# Hommel Tag
hommel.tag <- FALSE
if (mtp.type=='hommelq' || mtp.type=='hommell' || mtp.type=='hommellsi'|| mtp.type=='hommellsb') {
  hommel.tag <- TRUE
}
# Check method
if (mtp.type == 'romx' || mtp.type == 'rom') {
  critconst <- critconst_rom(alpha=alpha, n=n)
}
#
if (mtp.type == 'quickxi' || mtp.type == 'quickxb') {
  critconst <- critconst_quick(alpha=alpha, n=n)
}
#
rej.rslt <- rep(0,times=n) # arbitrary
#
if (gc.is.included) {
  IS.count.RSLT <- 0
  pkev$global.count.IS <- 0
}
#
init.count.RSLT <- 0
#
if (hommel.tag) {
  cnsnnc <- 0
}
#
###
#
#
p.srt <- sort(p.vec, method='quick', index.return=TRUE)
###
# Select a method
###
mtp.rslt <- switch(mtp.type,
                   # Quick (1+10)
                   quick = mtp.quick00b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   quick00i = mtp.quick00(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   quick01i = mtp.quick01(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   quick10i = mtp.quick10(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   quick11i = mtp.quick11(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   quickxi = mtp.quickxg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst, gc.is.included=gc.is.included),
                   quick00b = mtp.quick00b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   quick01b = mtp.quick01b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   quick10b = mtp.quick10b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   quick11b = mtp.quick11b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   quickxb = mtp.quickxbg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst, gc.is.included=gc.is.included),
                   # GTXR (1+6)
                   gtxr = mtp.gtxr0b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   gtxr0i = mtp.gtxr0(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   gtxr1ci = mtp.gtxr1c(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   gtxr2di = mtp.gtxr2d(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   gtxr0b = mtp.gtxr0b(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   gtxr1cb = mtp.gtxr1cb(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   gtxr2db = mtp.gtxr2db(pvec.sorted=p.srt$x, alpha=alpha, gc.is.included=gc.is.included),
                   # Rom (1+3)
                   rom = mtp.romxg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst),
                   rom1 = mtp.rom1(pvec.sorted=p.srt$x, alpha=alpha),
                   romx = mtp.romxg(pvec.sorted=p.srt$x, alpha=alpha, critconst=critconst),
                   roma = mtp.roma(pvec.sorted=p.srt$x, alpha=alpha),
                   # Hommel (1+4)
                   hommel = mtp.hommellsb(pvec.sorted=p.srt$x, alpha=alpha),
                   hommelq = mtp.hommelq(pvec.sorted=p.srt$x, alpha=alpha),
                   hommell = mtp.hommell(pvec.sorted=p.srt$x, alpha=alpha),
                   hommellsi = mtp.hommellsi(pvec.sorted=p.srt$x, alpha=alpha),
                   hommellsb = mtp.hommellsb(pvec.sorted=p.srt$x, alpha=alpha),
                   # Hochberg and Holm (3)
                   hochberg = mtp.hochberg(pvec.sorted=p.srt$x, alpha=alpha),
                   chochberg = mtp.chochberg(pvec.sorted=p.srt$x, alpha=alpha),
                   holm = mtp.holm(pvec.sorted=p.srt$x, alpha=alpha)
) # End of switch
#
#print(mtp.rslt$rej.idx)
#print(p.srt$x)
#
rej.rslt.srt <- c(rep(1,times=mtp.rslt$rej.idx), rep(0,n-mtp.rslt$rej.idx))
#
rej.rslt[p.srt$ix] <- rej.rslt.srt
#
init.count.RSLT <- mtp.rslt$init.count
#
if (gc.is.included) {
  IS.count.RSLT <- pkev$global.count.IS
  pkev$global.count.IS <- 0
}
if (hommel.tag) {
  IS.count.RSLT <- mtp.rslt$scnd.count
  cnsnnc <- mtp.rslt$consonance
}
#
count1 <- init.count.RSLT
if (gc.is.included || hommel.tag) {
  count2 <- IS.count.RSLT
  count <- count1 + count2
} else {
  count <- count1
}
if (hommel.tag) {
  consonance <- cnsnnc
} else {
  consonance <- 1
}
rslt <- list("rejection" = rej.rslt, "steps" = count, "consonance" = consonance)
}

