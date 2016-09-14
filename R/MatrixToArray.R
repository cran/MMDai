#' Matrix To Array
#' @description If the number of trials m = 1, the observed dataset may be a matrix format. This function transforms observed dataset into array format, which is used in DPMCM package
#' @param data_m - a matrix format of observed dataset
#' @param d - number of choices
#' @note  n - number of observations
#' @note  p - number of variables
#' @return data_a - an array format of observed dataset
#' @export
MatrixToArray<-function(data_m,d){

  # dimensional parameters
  n<-nrow(data_m)    # number of observations
  p<-ncol(data_m)    # number of variables

  # transformation
  data_a<-array(0,dim = c(n,p,d))
  for(i in 1:n){
    for(j in 1:p){
      data_a[i,j,data_m[i,j]]<-1
    }
  }

  # output
  return(data_a)

}
