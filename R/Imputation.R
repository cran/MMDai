#' Imputation
#' @description This function is used to perform multiple imputation for missing data given the joint distribution.
#' @param data - incomplete dataset
#' @param theta - vector of probability for each component
#' @param psi - specific probability for each variable in each component
#' @param method - methods for imputation, including "Sampling" and "MaxProb".
#' In "Sampling" method, sample missing values from the conditional distribution randomly.
#' In "MaxProb" method, impute missing values with maximal probability from the conditional distribution.
#' @return ImputedData - dataset has been imputated.
#' @import stats
#' @export


Imputation<-function(data, theta, psi, method = "MaxProb"){
  # dimensional parameters
  p<-ncol(data)    # number of variables
  n<-nrow(data)    # number of observations
  k<-length(theta) # number of components

  d<-rep(0,p)   # category for each variable
  for(j in 1:p){
    d[j]<-ncol(psi[[j]])
  }


  # multiple imputation
  ImputedData<-data

  # MaxProb method
  if(method == "MaxProb"){
    # imputation
    for(i in 1:n){
      if(any(is.na(data[i,]))){
        miss<-which(is.na(data[i,]))
        obs<-which(!is.na(data[i,]))

        for(j in miss){
          CondProb<-rep(0,d[j])
          for(c in 1:d[j]){
            for(h in 1:k){
              prob<-theta[h]*psi[[j]][h,c]
              for(jj in obs){
                prob<-prob*psi[[jj]][h,data[i,jj]]
              }
              CondProb[c]<-CondProb[c]+prob
            }
          }

          ImputedData[i,j]<-which.max(CondProb)

        }

      }
    }
  }


  # Sampling method
  if(method == "Sampling"){
    # imputation
    for(i in 1:n){
      if(any(is.na(data[i,]))){
        miss<-which(is.na(data[i,]))
        obs<-which(!is.na(data[i,]))

        for(j in miss){
          CondProb<-rep(0,d[j])
          for(c in 1:d[j]){
            for(h in 1:k){
              prob<-theta[h]*psi[[j]][h,c]
              for(jj in obs){
                prob<-prob*psi[[jj]][h,data[i,jj]]
              }
              CondProb[c]<-CondProb[c]+prob
            }
          }

          ImputedData[i,j]<-which(rmultinom(1,1,CondProb)==1)

        }

      }
    }
  }

  # output
  return(ImputedData)

}
