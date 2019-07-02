#' @name p.adjust
#' @title Adjust P-values for Multiple Test Procedures
#'
#' @description Given a set of p-values, returns adjusted p-values, including the hybrid Hochberg-Hommel procedure (Gou et al., 2014) and Quick procedure (Gou and Zhang, 2020).
#'
#'
#'
#' @param p vector of p-values.
#' @param method multiplicity correction method, "gtxr" is the hybrid Hochberg-Hommel method, "quick" is the Quick method. Other methods include:"holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none" from the standard R function p.adjust.
#' @param n number of p-values.
#' @return a vector of corrected p-values.
#' @author Jiangtao Gou
#' @references
#' Gou, J., Tamhane, A. C., Xi, D., and Rom, D. (2014). A class of improved hybrid Hochberg-Hommel type step-up multiple test procedures. \emph{Biometrika} \bold{101}, 899-911. <https://dx.doi.org/10.1093/biomet/asu032>
#' 
#' Gou, J., and Zhang, F. (2020). Quick multiple test procedures and p-value adjustments. Technical report.
#' 
#' @details
#' Given a set of p-values, returns p-values adjusted using one of several methods. The default method is "gtxr". Another option is "quick". Other adjustment methods have been included in function p.adjust in R package stats.
#' @examples
#' library(elitism)
#'  pvalues.raw <- c(0.002,0.007,0.005,0.024,0.022,0.009,0.007,0.036,0.060,0.035)
#'  p.adj.hoch <- elitism::p.adjust(pvalues.raw, method = "hochberg")
#'  p.adj.quick <- elitism::p.adjust(pvalues.raw, method = "quick")
#'  p.adj.gtxr <- elitism::p.adjust(pvalues.raw, method = "gtxr")
#' @seealso \code{stats::p.adjust}
#' @export
#' @import stats
#'
# utils::globalVariables(c("global.count.IS", "global.count.FS"))
#
#global.count.IS <<- 0
#global.count.FS <<- 0
#gc.is.included <<- TRUE
#
p.adjust <- function (p, method = "gtxr", n = length(p)) {
  if (method == "holm" || method == "hochberg" || method == "hommel" || method == "bonferroni" || method == "BH" || method ==  "BY" || method == "fdr" || method ==  "none") {
    result <- stats::p.adjust(p,method,n);
    return(result);
  } else if (method == "gtxr") {
    cvec <- (2:(n+1))/(2*(1:n));
    dvec <- 1/(1:n);
    p.sorted.adjust <- rep(0,n);
    result <- rep(0,n);
    p.sorted <- sort(p,decreasing = TRUE, index.return = TRUE);
    for (i in 1:n) {
      p.sorted.adjust[i] <- min(pmax(p.sorted$x[1:i]/cvec[1:i], p.sorted$x[i]/dvec[1:i]));
    }
    result[p.sorted$ix] <- p.sorted.adjust[1:n];
    return(result);
  }  else if (method == "quick") {
    cvec <- c(1, rep(1/2,n-1));
    dvec <- 1/(1:n);
    p.sorted.adjust <- rep(0,n);
    result <- rep(0,n);
    p.sorted <- sort(p,decreasing = TRUE, index.return = TRUE);
    for (i in 1:n) {
      p.sorted.adjust[i] <- min(pmax(p.sorted$x[1:i]/cvec[1:i], p.sorted$x[i]/dvec[1:i]));
    }
    result[p.sorted$ix] <- p.sorted.adjust[1:n];
    return(result);
  } else {
    print("Methods include: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none, gtxr, quick.");
    return(p);
  }
}
#
# End of Fucntion p.adjust
#

##########
#
# Function Extra
#
# Fucntion p.adjust
#
# Jiangtao Gou
# 2016-09-24
# 2019-06-25 add 'quick'
# Given a set of p-values, returns p-values adjusted using one of several methods. The default method is "gtxr".
# Example: p.adjust.gtxr(c(0.002,0.007,0.005,0.024,0.022,0.009,0.007,0.036,0.060,0.035), method = "gtxr")
# Ref: Gou J, Tamhane AC, Xi D, Rom D. A class of improved hybrid Hochberg-Hommel type step-up multiple test procedures. Biometrika. 2014; 101:899-911.
#
