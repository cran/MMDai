#' Matrix To Factor
#' @description This function is used to transform data from matrix into factor.
#' @param data_m - data in matrix formation
#' @param d - number of categoriesin each variable
#' @return data_f - data in factor formation
#' @export
MatrixToFactor<-function(data_m,d){

  # dimensional parameters
  p<-ncol(data_m)    # number of variables

  # transformation
  data_f<-data.frame(data_m)
  for(j in 1:p){
    data_f[,j]<-as.factor(data_f[,j])
  }

  # output
  return(data_f)

}
