#' Imputation
#' @description This function is used to perform multiple imputation for missing data given the estimates of theta and psi.
#' @param data - incomplete dataset
#' @param theta - a vector that sum to 1, denotes probability of latent class.
#' @param psi - an array with dimension c(k,p,d), specific probability for each variables in each component.
#' @param m - number of trials
#' @param method - methods for imputation, including Sampling and Expectation.
#' @note  k - number of components
#' @note  n - number of observations
#' @note  p - number of variables
#' @note  d - number of choices. Here d are identical for all variables.
#' @return CompleteData - dataset has been imputated.
#' @import stats
#' @examples
#' data("IncompleteData")
#' theta <- IncompleteData$theta
#' psi <- IncompleteData$psi
#' Imputation(data = IncompleteData$data,theta,psi,m = 10,method = "Sampling")
#' @export


Imputation<-function(data,theta,psi,m,method = "Sampling"){
  # dimensional parameters
  n<-dim(data)[1]    # number of samples
  p<-dim(data)[2]    # number of variables
  d<-dim(data)[3]    # number of categories
  k<-length(theta)   # number of components

  # multiple imputation
  CompleteData<-data

  # Sampling method
  if(method == "Sampling"){
    for(i in 1:n){
      for(j in 1:p){
        if(is.na(data[i,j,1])){
          h<-rmultinom(1,1,theta)
          CompleteData[i,j,]<-rmultinom(1,m,psi[h,j,])
        }
      }
    }
  }

  # Expectation method
  if(method == "Expectation"){
    jointpsi<-matrix(0,nrow = p, ncol = d)
    for(h in 1:k){
      jointpsi<-jointpsi+theta[h]*psi[h,,]
    }
    for(i in 1:n){
      for(j in 1:p){
        if(is.na(data[i,j,1])){
          CompleteData[i,j,]<-rmultinom(1,m,jointpsi[j,])
        }
      }
    }
  }

  # output
  return(list(CompleteData = CompleteData))

}
