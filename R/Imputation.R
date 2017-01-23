#' Imputation
#' @description This function is used to perform multiple imputation for missing data given the joint distribution.
#' @param data - incomplete dataset
#' @param JT - joint distribution used
#' @param method - methods for imputation, including "Sampling" and "MaxProb".
#' In "Sampling" method, sample missing values from the conditional distribution randomly.
#' In "MaxProb" method, impute missing values with maximal probability from the conditional distribution.
#' @return ImputedData - dataset has been imputated.
#' @import stats
#' @export


Imputation<-function(data, JT, method = "MaxProb"){
  # dimensional parameters
  p<-ncol(data)    # number of variables
  n<-nrow(data)    # number of observations

  # multiple imputation
  ImputedData<-data

  # MaxProb method
  if(method == "MaxProb"){
    # imputation
    for(i in 1:n){
      if(any(is.na(data[i,]))){
        CondJT<-JT
        for(j in 1:p){
          if(!is.na(data[i,j])){
            CondJT<-CondJT[CondJT[,j]==data[i,j],]
          }
        }

        imp<-CondJT[which.max(CondJT[,p+1]),1:p]

        ImputedData[i,]<-imp
      }
    }
  }


  # Sampling method
  if(method == "Sampling"){
    # imputation
    for(i in 1:n){
      if(any(is.na(data[i,]))){
        CondJT<-JT
        for(j in 1:p){
          if(!is.na(data[i,j])){
            CondJT<-CondJT[CondJT[,j]==data[i,j],]
          }
        }

        index<-which(rmultinom(1,1,CondJT[,p+1])==1)

        ImputedData[i,]<-CondJT[index,1:p]
      }
    }
  }

  # output
  return(ImputedData)

}
