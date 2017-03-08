#' Bhattacharyya distance
#' @param JT1 - first joint distribution table
#' @param JT2 - second joint distribution table

#' @export

DBhattacharyya<-function(JT1, JT2){

  # dimension parameter
  p<-ncol(JT1)-1

  # joint probability
  prob1<-JT1[,p+1]
  prob2<-JT2[,p+1]

  # Bhattacharyya coefficient
  BC<-sum(sqrt(prob1*prob2))

  # Bhattacharyya distance
  DB<--log(BC)

  # output
  return(DB)
}
