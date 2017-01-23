#' Joint distribution table estimation
#' This function is used to esitmate joint distribution table based on theta and psi
#' @param theta - vector of probability for each component
#' @param psi - specific probability for each variable in each component
#' @return table - joint distribution table
#' @export

JointTableEst<-function(theta,psi){
  # dimension parameter
  k<-length(theta) # number of components
  p<-length(psi)   # number of variables

  d<-rep(0,p)      # number of categories
  for(j in 1:p){
    d[j]<-ncol(psi[[j]])
  }

  # empty table
  table<-matrix(0,nrow = prod(d),ncol = p+1)

  # joint input
  for(j in 1:p){
    if(j==p)
      times<-1
    else
      times<-prod(d[(j+1):p])

    table[,j]<-rep(c(rep(1:d[j],rep(prod(d[1:j-1]),d[j]))),times)
  }

  # joint probability
  for(dd in 1:prod(d)){
    for(h in 1:k){
      # probability in each component
      temp<-theta[h]
      for(j in 1:p){
        temp<-temp*psi[[j]][h,table[dd,j]]
      }

      table[dd,p+1]<-table[dd,p+1]+temp

    }
  }


  # output
  return(table)
}
