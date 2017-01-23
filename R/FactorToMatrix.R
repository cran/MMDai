#' Factor To Matrix
#' @description This function is used to transform data from factor into matrix.
#' @param data_f - data in factor formation
#' @param d - number of choices in each variable
#' @return data_m - data in matrix formation
#' @export
FactorToMatrix<-function(data_f,d){

  # dimensional parameters
  p<-ncol(data_f)    # number of variables
  n<-nrow(data_f)    # number of observations

  # transformation
  data_m<-matrix(0,nrow = n,ncol = p)
  for(j in 1:p){
    data_m[,j]<-as.numeric(data_f[,j])
  }

  # output
  return(data_m)

}
