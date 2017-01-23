#' mutual information
#' @description This function is used to calculate the mutual information (MI) and maximal information coefficient (MIC) between two variables.
#' @param JT - joint distribution table used
#' @param v1 - first variable
#' @param v2 - second variable
#' @return MI - mutual information between two variables
#' @return MIC - maximal information coefficient between two variables
#' @export

MutualInfo<-function(JT, v1, v2){

  # dimension parameter
  p<-ncol(JT)-1    # number of variables

  d<-rep(0,p)       # number of categories
  for(j in 1:p){
    d[j]<-length(unique(JT[,j]))
  }


  # Mutual Information
  Hv1<-0
  for(c1 in 1:d[v1]){
    pv1<-sum(JT[JT[,v1]==c1,p+1])
    if(pv1==0)
      next
    else
      Hv1<-Hv1-pv1*log(pv1)
  }

  Hv2<-0
  for(c2 in 1:d[v2]){
    pv2<-sum(JT[JT[,v2]==c2,p+1])
    if(pv2==0)
      next
    else
      Hv2<-Hv2-pv2*log(pv2)
  }

  Hv1v2<-0
  for(c1 in 1:d[v1]){
    for(c2 in 1:d[v2]){
      pv1v2<-sum(JT[JT[,v1]==c1 & JT[,v2]==c2,p+1])

      if(pv1v2==0)
        next
      else
        Hv1v2<-Hv1v2-pv1v2*log(pv1v2)
    }
  }

  # mutual information
  MI<-Hv1+Hv2-Hv1v2

  # maximal information coefficient
  MIC<-MI/(min(Hv1,Hv2))

  # output
  return(list(MI = MI, MIC = MIC))
}
