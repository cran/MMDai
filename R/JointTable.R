#' Generate joint distribution table
#' @description This function is used to generate general joint distribution table
#' @param p - number of variables
#' @param d - number of categories
#' @param alpha - prior parameter
#' @return JT - generated joint distribution table
#' @export

JointTable<-function(p, d, alpha = 0.5){
  # empty joint table
  JT<-matrix(0,nrow = prod(d),ncol = p+1)

  # joint input
  for(j in 1:p){
    if(j==p)
      times<-1
    else
      times<-prod(d[(j+1):p])

    JT[,j]<-rep(c(rep(1:d[j],rep(prod(d[1:j-1]),d[j]))),times)
  }

  # joint probability
  probability<-rdirichlet(1,rep(alpha,prod(d)))
  JT[,p+1]<-probability

  # output
  return(JT)
}
