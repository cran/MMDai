#' Generate random dataset
#' @description This function is used to generate random datasets. It helps to evaluate the performance of this package.
#' @param n - number of samples
#' @param p - number of variables
#' @param d - a vector which denotes the number of categories for each variable. It could be distinct among variables.
#' @param method - the method to generate random datasets, including "MixMultinomial" and "General".
#' In "MixMultinomial" method, the dataset follows finite mixture of multinomial distribution
#' In "General" method, the dataset follows general joint distribution.
#' @param k - parameter used in "MixMultinomial" method
#' @param theta - parameter used in "MixMultinomial" method
#' @param psi - parameter used in "MixMultinomial" method
#' @param JT - parameter used in "General" method
#' @return data - generated random dataset, a matrix with n rows and p columns.
#' @importFrom DirichletReg rdirichlet
#' @examples
#' # dimension parameters
#' n<-200; p<-5; d<-rep(2,p);
#' # generate complete data
#' Complete<-GenerateData(n, p, d, k = 3, method = "MixMultinomial")
#' @export
GenerateData<-function(n, p, d, method = "MixMultinomial", k = 3, theta = rdirichlet(1,rep(10,k)), psi = InitialPsi(p,d,k), JT = JointTable(p,d)){
  # initial dataset
  data<-matrix(0,nrow = n,ncol = p)

  if(method == "MixMultinomial"){
    # generate complete data
    for(i in 1:n){
      z<-which(rmultinom(1,1,theta)==1)
      for(j in 1:p){
        data[i,j]<-which(rmultinom(1,1,psi[[j]][z,])==1)
      }
    }
  }

  if(method == "General"){
    # generate complete data
    for(i in 1:n){
      index<-which(rmultinom(1,1,JT[,p+1])==1)
      data[i,]<-JT[index,1:p]
    }
  }

  # output
  return(data)

}
