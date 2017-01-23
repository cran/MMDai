#' Correlation estimation
#' @description Calculate correlation matrix among variables
#' @param JT - joint distribution table used
#' @return cor - correlation matrix
#' @export

Correlation<-function(JT){
  # dimensional parameters
  p<-ncol(JT)-1     # number of variables

  d<-rep(0,p)       # number of categories
  for(j in 1:p){
    d[j]<-length(unique(JT[,j]))
  }

  # correlation matrix
  cor<-matrix(0,nrow = p,ncol = p)
  for(v1 in 1:p){
    for(v2 in 1:p){
      if(v1==v2){
        Ev1v2<-0
        for(c1 in 1:d[v1]){
          Ev1v2<-Ev1v2+c1*c1*sum(JT[JT[,v1]==c1,p+1])
        }
      }else{
        Ev1v2<-0
        for(c1 in 1:d[v1]){
          for(c2 in 1:d[v2]){
            Ev1v2<-Ev1v2+c1*c2*sum(JT[JT[,v1]==c1 & JT[,v2]==c2,p+1])
          }
        }
      }

      Ev1<-0
      Ev1v1<-0
      for(c1 in 1:d[v1]){
        Ev1v1<-Ev1v1+c1*c1*sum(JT[JT[,v1]==c1,p+1])
        Ev1<-Ev1+c1*sum(JT[JT[,v1]==c1,p+1])
      }

      Ev2v2<-0
      Ev2<-0
      for(c2 in 1:d[v2]){
        Ev2v2<-Ev2v2+c2*c2*sum(JT[JT[,v2]==c2,p+1])
        Ev2<-Ev2+c2*sum(JT[JT[,v2]==c2,p+1])
      }


      # variance
      Varv1<-Ev1v1-Ev1^2
      Varv2<-Ev2v2-Ev2^2

      # correlation
      cor[v1,v2]<-(Ev1v2-Ev1*Ev2)/sqrt(Varv1*Varv2)

    }
  }

  return(cor)

}
