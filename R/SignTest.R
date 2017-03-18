#' Sign test
#' @description Test for consistent differences between pairs of observations
#' @param x - first variable
#' @param y - second variable
#' @param method - In "One-tail" method, the p-value is computed by Pr(x>y).
#' In "Two-tail" method, the p-value is computed by Pr(x>y)+Pr(x<y).
#' @return p_value - p value in the test
#' @export

SignTest<-function(x, y, method = "One-tail"){

  n<-length(x)
  diff<-x-y
  positive<-length(diff[diff>0])

  # one-tail test
  if(method == "One-tail"){
    p_value<-binom.test(positive,n)$p.value/2
  }

  # two-tail test\
  if(method == "Two-tail"){
    p_value<-binom.test(positive,n)$p.value
  }

  # output
  return(list(p_value = p_value))

}
