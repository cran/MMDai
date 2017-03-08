#' Sample joint distribution table
#' @description This function is used to calcualte sample joint distribution table given samples.
#' @param data - samples are observed
#' @param d - number of categories
#' @return SJT - sample joint distribution table
#' @export

SampleJointTable<-function(data, d){
  # dimension parameters
  n<-nrow(data)
  p<-ncol(data)

  # empty joint table
  SJT<-JointTable(p, d)
  SJT[,p+1]<-0

  # observation count
  for(i in 1:n){
    index<-rep(T,nrow(SJT))
    for(j in 1:p){
      index<-index & SJT[,j]==data[i,j]
    }
    SJT[index,p+1]<-SJT[index,p+1]+1
  }

  SJT[,p+1]<-SJT[,p+1]/n

  # output
  return(SJT)
}
