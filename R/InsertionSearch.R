# Insertion Search
# 2019-02-18 
# 2019-03-01 revised
# Jiangtao Gou
# 
# 2019-03-01: add an additional global counter. Input secondStage = TRUE is the default setting for the second stage of GTXR, Hommel, Quick procedure. When secondStage = False, we use the other global counter for the Quick procedure. 
#
# Example1 : InsertionSearch(a=rev(c(0.97,0.83,0.79,0.73,0.61,0.59,0.47,0.37,0.31,0.29,0.23,0.19,0.13,0.07,0.03)), value=0.40, low=1, high=15)
# Example 2: InsertionSearch(a=c(0.01, 0.07, 0.13, 0.13, 0.13, 0.13, 0.13, 0.43, 0.67, 0.97), value=0.43, low=1, high=10)
# Example 3: InsertionSearch(a=c(0.01, 0.07, 0.13, 0.13, 0.13, 0.13, 0.13, 0.43, 0.67, 0.97), value=0.13, low=1, high=10)
#
InsertionSearch <- function (a, value, low, high, gc.is.included=FALSE, secondStage=TRUE) {
  #
  if (gc.is.included) {
    if (secondStage) {
      pkev$global.count.IS <- pkev$global.count.IS + 1
    } else {
      pkev$global.count.FS <- pkev$global.count.FS + 1
    } # End of if (secondStage)
  } # End of if (gc.is.included)
  #
  if (a[low] > value) {
    return (low-1)
  } else if (a[high] <= value) {
    return (high)
  }
  mid <- low + (value-a[low])/(a[high]-a[low])*(high-low)
  if (mid < low) {
    return (low-1)
  } else if (mid >= high) {
    return (high)
  } else {
    mid <- round(mid)
  }
  if(a[mid]>value) {
    return (InsertionSearch(a,value,low,mid-1, gc.is.included,secondStage))
  } else if(a[mid] <= value) {
    return (InsertionSearch(a,value,mid+1,high, gc.is.included,secondStage))
  }
}
#
# Ref: 
# <https://www.jianshu.com/p/d78dd2b32580>
# <https://d.cosx.org/d/12252-12252>
# <https://www.zhihu.com/question/20387324>
#