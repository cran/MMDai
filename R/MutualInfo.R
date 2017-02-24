#' mutual information
#' @description This function is used to calculate the mutual information (MI) and maximal information coefficient (MIC) between two variables.
#' @param theta - vector of probability for each component
#' @param psi - specific probability for each variable in each component
#' @param v1 - first variable
#' @param v2 - second variable
#' @return MI - mutual information between two variables
#' @return MIC - maximal information coefficient between two variables
#' @export

MutualInfo<-function(theta, psi, v1, v2){
  # dimensional parameters
  p<-length(psi)    # number of variables
  k<-length(theta) # number of components

  d<-rep(0,p)   # category for each variable
  for(j in 1:p){
    d[j]<-ncol(psi[[j]])
  }

  # Mutual Information
  Hv1<-0
  for(c1 in 1:d[v1]){
    prob<-0
    for(h in 1:k){
      prob<-prob+theta[h]*psi[[v1]][h,c1]
    }

    Hv1<-Hv1-prob*log(prob)
  }

  Hv2<-0
  for(c2 in 1:d[v2]){
    prob<-0
    for(h in 1:k){
      prob<-prob+theta[h]*psi[[v2]][h,c2]
    }

    Hv2<-Hv2-prob*log(prob)
  }

  Hv1v2<-0
  for(c1 in 1:d[v1]){
    for(c2 in 1:d[v2]){
      prob<-0
      for(h in 1:k){
        prob<-prob+theta[h]*psi[[v1]][h,c1]*psi[[v2]][h,c2]
      }

      Hv1v2<-Hv1v2-prob*log(prob)
    }
  }

  # mutual information
  MI<-Hv1+Hv2-Hv1v2

  # maximal information coefficient
  MIC<-MI/(min(Hv1,Hv2))

  # output
  return(list(MI = MI, MIC = MIC))
}
