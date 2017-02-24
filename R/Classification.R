#' Classification
#' @description This function is used to classify each observation.
#' @param data - Imputed dataset
#' @param theta - vector of probability for each component
#' @param psi - specific probability for each variable in each component
#' @return class - class for each observation.
#' @import stats
#' @export

Classification<-function(data, theta, psi){

  # dimensional parameters
  p<-ncol(data)    # number of variables
  n<-nrow(data)    # number of observations
  k<-length(theta) # number of components


  # classification
  class<-rep(0,n)

  for(i in 1:n){
    prob<-rep(0,k)

    for(h in 1:k){
      prob[h]<-theta[h]
      for(j in 1:p){
        prob[h]<-prob[h]*psi[[j]][h,data[i,j]]
      }
    }

    class[i]<-which.max(prob)

  }

  # output
  return(class)

}


