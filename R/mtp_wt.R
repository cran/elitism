#' @name mtp.wt
#' @title Multiple Test Procedures with known true significances and non-significances 
#'
#' @description Given a set of p-values, a set of indicators whether the corresponding hypothesis is true significance or true null, and the level of significance, returns a summary of false/true positive/negative, including the hybrid Hochberg-Hommel procedure (Gou et al., 2014) and Quick procedure (Gou and Zhang, 2020).
#'
#' @param p vector of p-values. 
#' @param indctr.sig vector of indicators, 1 stands for true significance and 0 stands for true null.
#' @param alpha the level of significance.
#' @param method multiplicity correction method, including the Holm procedure ("holm"), the Hochberg procedure ("hochberg", "chochberg"), the Hommel procedure ("hommel", "hommelq", "hommell", "hommellsi",  "hommellsb"), the Rom procedure ("rom", "rom1", "roma", "romx"), the Gou-Tamhane-Xi-Rom procedure ("gtxr", "gtxr0i", "gtxr1ci", "gtxr2di", "gtxr0b", "gtxr1cb", "gtxr2db"), and the Quick procedure ("quick", "quick00i", "quick01i", "quick10i", "quick11i", "quickxi", "quick00b", "quick01b", "quick10b", "quick11b", "quickxb"). 
#' @param n number of p-values.
#' @return a list, including five integers and a binary indicator: number of false positives, number of true negatives, number of true positives, number of false negatives, the number of total comparisons, and an indicator of consonance. 
#' @author Jiangtao Gou
#' @author Fengqing (Zoe) Zhang
#' @references
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure. \emph{Scandinavian Journal of Statistics} \bold{6}, 65-70. <http://www.jstor.org/stable/4615733>
#' 
#' Hochberg, Y. and Tamhane, A. C. (1987). Multiple Comparison Procedures. John Wiley and Sons, New York. <http://dx.doi.org/10.1002/9780470316672>
#' 
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. \emph{Biometrika} \bold{75}, 800-802. <http://dx.doi.org/10.2307/2336325>
#' 
#' Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. \emph{Biometrika} \bold{75}, 383-386. <http://dx.doi.org/10.1093/biomet/75.2.383>
#' 
#' Rom, D. M. (1990). A sequentially rejective test procedure based on a modified Bonferroni inequality. \emph{Biometrika} \bold{77}, 663-665. <https://dx.doi.org/10.1093/biomet/77.3.663>
#' 
#' Wright, S. P. (1992). Adjusted p-values for simultaneous inference. \emph{Biometrics} \bold{48}, 1005-1013. <http://www.jstor.org/stable/2532694>
#' 
#' Gou, J., Tamhane, A. C., Xi, D., and Rom, D. (2014). A class of improved hybrid Hochberg-Hommel type step-up multiple test procedures. \emph{Biometrika} \bold{101}, 899-911. <https://dx.doi.org/10.1093/biomet/asu032>
#' 
#' Gou, J., and Tamhane, A. C. (2014). On generalized Simes critical constants. \emph{Biometrical Journal} \bold{56}, 1035-1054. <http://dx.doi.org/10.1002/bimj.201300258>
#' 
#' Gou, J., and Tamhane, A. C. (2018). Hochberg procedure under negative dependence. \emph{Statistica Sinica} \bold{28}, 339-362. <http://www3.stat.sinica.edu.tw/statistica/J28N1/J28N116/J28N116.html>
#' 
#' Tamhane, A. C., and Gou, J. (2018). Advances in p-value based multiple test procedures. \emph{Journal of Biopharmaceutical Statistics} \bold{28}, 10-27. <https://dx.doi.org/10.1080/10543406.2017.1378666>
#' 
#' Meijer, R. J., Krebs, T. J. P., and Goeman, J. J. (2019). Hommel's procedure in linear time. \emph{Biometrical Journal} \bold{61}, 73-82. <https://dx.doi.org/110.1002/bimj.201700316>
#' 
#' Gou, J., and Zhang, F. (2020). Quick multiple test procedures and p-value adjustments. Technical report.
#' 
#' @details
#' Given a set of p-values with a binary vector of true signficances, where 1 stands for true significances, and 0 stands for true nulls. There are six families of procedures. 
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
#'  indctr.sig <- c(1, 0, 0, 0, 1,  1, 1, 1, 0, 0)
#'  pkev <- new.env(); pkev$global.count.IS <- 0; pkev$global.count.FS <- 0;
#'  summary.hoch <- mtp.wt(pvalues.raw, indctr.sig, alpha = 0.025, method = "hochberg")
#'  pkev <- new.env(); pkev$global.count.IS <- 0; pkev$global.count.FS <- 0;
#'  summary.quick <- mtp.wt(pvalues.raw, indctr.sig, alpha = 0.025, method = "quick")
#'  pkev <- new.env(); pkev$global.count.IS <- 0; pkev$global.count.FS <- 0;
#'  summary.gtxr <- mtp.wt(pvalues.raw, indctr.sig, alpha = 0.025, method = "gtxr")
#' @seealso \code{elitism::mtp}
#' @export
#' @import stats
#'
#'
#'
#'
#'
#'
#'
#'
#'
#utils::globalVariables(c("global.count.IS", "global.count.FS"))
#
#global.count.IS <<- 0
#global.count.FS <<- 0
#gc.is.included <<- TRUE
#
# pkev <- new.env(); pkev$global.count.IS <- 0; pkev$global.count.FS <- 0
#
mtp.wt <- function (p, indctr.sig, alpha = 0.05, method = "gtxr", n = length(p)) {
  # pkev <- new.env()
  #pkev$global.count.IS <- 0
  #pkev$global.count.FS <- 0
  #
  indctr.null <- 1 - indctr.sig
  n.null <- sum(indctr.null)
  n.altr <- sum(indctr.sig)
  #rslt <- elitism::mtp(p, alpha, method, n)
  rslt <- mtp(p, alpha, method, n)
  rej.rslt <- rslt$rejection
  fpos<- sum(rej.rslt*indctr.null)
  tneg<- n.null - fpos
  tpos <- sum(rej.rslt*(1-indctr.null))
  fneg <- n.altr - tpos
  #
  rsltwt <- list("falsePos" = fpos, "trueNeg" = tneg, "truePos" = tpos, "falseNeg" = fneg,"steps" = rslt$steps, "consonance" = rslt$consonance)
}