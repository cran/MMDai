#' Chi-Square test
#' @description Test the independence of two variables
#' @param JT - joint distribution table used
#' @param v1 - first variable
#' @param v2 - second variable
#' @return chi_square - chi_square test statistic
#' @return p_value - p value in the test
#' @export

ChiSquareTest<-function(JT,v1,v2){
  # dimensional parameters
  p<-ncol(JT)-1      # number of variables
  d<-rep(0,2)        # number of categories
  d[1]<-length(unique(JT[,v1]))
  d[2]<-length(unique(JT[,v2]))

  # observation table
  obs<-matrix(0,nrow = d[1],ncol = d[2])
  for(d1 in 1:d[1]){
    for(d2 in 1:d[2]){
      obs[d1,d2]<-sum(JT[JT[,v1]==d1 & JT[,v2]==d2,p+1])
    }
  }

  # expectation table
  exp<-matrix(0,nrow = d[1],ncol = d[2])
  for(d1 in 1:d[1]){
    for(d2 in 1:d[2]){
      exp[d1,d2]<-sum(obs[d1,])*sum(obs[,d2])
    }
  }


  # chi-square
  obs<-obs*100
  exp<-exp*100
  chi_square<-sum((obs-exp)^2/exp)
  df<-(d[v1]-1)*(d[v2]-1)
  p_value<-1-pchisq(chi_square,df)

  # output
  return(list(chi_square = chi_square,p_value = p_value))



}
