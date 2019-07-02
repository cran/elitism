# Binary Search 
# 2019-03-05
# Jiangtao Gou
#
# Example 3: BinarySearch(a=c(0.01, 0.07, 0.13, 0.13, 0.13, 0.13, 0.13, 0.43, 0.67, 0.97), value=0.13, low=1, high=10)
#
# BinarySearch(a=c(9.410985e-05, 8.030631e-04, 3.215159e-03, 1.193564e-02), value=0.025, low=1, high=4)

#
BinarySearch <- function(a, value, low, high, gc.is.included=FALSE, secondStage=TRUE) {
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
  #
  mid = round(low+(high-low)/2)
  #
  #
  if(a[mid] > value) {
    return (BinarySearch(a,value,low,mid-1, gc.is.included,secondStage))
  } else if(a[mid] <= value) {
    return (BinarySearch(a,value,mid+1,high, gc.is.included,secondStage))
  }
  #
}


#<https://www.jianshu.com/p/d78dd2b32580>
