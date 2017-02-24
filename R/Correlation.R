#' Correlation estimation
#' @description Calculate Pearson's correlation between two specific variables
#' @param theta - vector of probability for each component
#' @param psi - specific probability for each variable in each component
#' @param v1 - first variable
#' @param v2 - second variable
#' @return Cor - correlation
#' @export

Correlation<-function(theta, psi, v1 = 1, v2 = 2){
  # dimensional parameters
  p<-length(psi)    # number of variables
  k<-length(theta) # number of components

  d<-rep(0,p)   # category for each variable
  for(j in 1:p){
    d[j]<-ncol(psi[[j]])
  }

  # correlation matrix
  Ev1v2<-0
  for(c1 in 1:d[v1]){
    for(c2 in 1:d[v2]){
      prob<-0
      for(h in 1:k){
        prob<-prob+theta[h]*psi[[v1]][h,c1]*psi[[v2]][h,c2]
      }
      Ev1v2<-Ev1v2+prob*c1*c2
    }
  }

  Ev1<-0
  Ev1v1<-0
  for(c1 in 1:d[v1]){
    prob<-0
    for(h in 1:k){
      prob<-prob+theta[h]*psi[[v1]][h,c1]
    }
    Ev1<-Ev1+prob*c1
    Ev1v1<-Ev1v1+prob*c1*c1
  }

  Ev2<-0
  Ev2v2<-0
  for(c2 in 1:d[v2]){
    prob<-0
    for(h in 1:k){
      prob<-prob+theta[h]*psi[[v2]][h,c2]
    }
    Ev2<-Ev2+prob*c2
    Ev2v2<-Ev2v2+prob*c2*c2
  }

  Varv1<-Ev1v1-(Ev1)^2
  Varv2<-Ev2v2-(Ev2)^2


  if(v1==v2)
    Cor<-(Ev1v1-Ev1*Ev1)/sqrt(Varv1*Varv1)
  else
    Cor<-(Ev1v2-Ev1*Ev2)/sqrt(Varv1*Varv2)

  return(Cor)

}
