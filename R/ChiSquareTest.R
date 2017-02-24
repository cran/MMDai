#' Chi-Square test
#' @description Test the independence of two variables
#' @param theta - vector of probability for each component
#' @param psi - specific probability for each variable in each component
#' @param v1 - first variable
#' @param v2 - second variable
#' @return chi_square - chi_square test statistic
#' @return p_value - p value in the test
#' @export

ChiSquareTest<-function(theta, psi, v1 = 1, v2 =2){
  # dimensional parameters
  p<-length(psi)    # number of variables
  k<-length(theta) # number of components

  d<-rep(0,p)   # category for each variable
  for(j in 1:p){
    d[j]<-ncol(psi[[j]])
  }

  # observation table
  obs<-matrix(0,nrow = d[v1],ncol = d[v2])
  for(c1 in 1:d[v1]){
    for(c2 in 1:d[v2]){
      for(h in 1:k){
        obs[c1,c2]<-obs[c1,c2]+theta[h]*psi[[v1]][h,c1]*psi[[v2]][h,c2]
      }
    }
  }

  # expectation table
  exp<-matrix(0,nrow = d[v1],ncol = d[v2])
  for(c1 in 1:d[v1]){
    for(c2 in 1:d[v2]){
      exp[c1,c2]<-sum(obs[c1,])*sum(obs[,c2])
    }
  }


  # chi-square
  obs<-obs*100
  exp<-exp*100
  chi_square<-sum((obs-exp)^2/exp)
  df<-(d[v1]-1)*(d[v2]-1)
  p_value<-1-pchisq(chi_square,df)

  # output
  return(list(chi_square = chi_square, p_value = p_value))

}
