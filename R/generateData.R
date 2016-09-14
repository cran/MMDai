#' Generate random mixture multinomial dataset
#' @description This function is used to generate random mixture multinomial dataset. It helps to evaluate our estimatior.
#' @param n - number of observations
#' @param p - number of variables
#' @param d - number of categories
#' @param m - number of trials
#' @param k - number of components
#' @param theta - a vector that sum to 1, denotes probability of latent components. Default value is Dirichlet random number with equal weight.
#' @param psi - an array with dimension c(k,p,d), specific probability for each variables in each component. Default value is Dirichlet random number with equal weight.
#' @param miss - rate of missing data in all observations. Missing mechanism is Missing Completely At Random (MCAR). The value of miss is from 0 to 1. Default value is 0. Missing data is denoted by NA.
#' Here assume all trials in data[i,j,] are missing or observed at the same time.
#' @return data - random dataset generated.
#' @return theta - values of theta are used to generate random dataset.
#' @return psi - values of psi are used to generate random dataset.
#' @importFrom DirichletReg rdirichlet
#' @examples
#' n<-100; k<-2; d<-2; m<-10; p<-2;
#' completedata<-generateData(n,p,d,m,k)
#' incompletedata<-generateData(n,p,d,m,k,miss = 0.25)
#' @export
generateData<-function(n,p,d,m,k,theta = rdirichlet(1,rep(1,k)), psi = array(rdirichlet(k*p,rep(1,d)),dim = c(k,p,d)),miss = 0){
  # generate complete data
  data<-array(0,dim = c(n,p,d))
  z<-rep(0,n)  # latent component
  for(i in 1:n){
    z[i]<-which(rmultinom(1, 1, theta)==1)
    for(j in 1:p){
      data[i,j,]<-rmultinom(1, m, psi[z[i],j,])
    }
  }

  # missing completely at random
  missing<-sample(1:(n*p),size = floor(n*p*miss),replace = F)
  for(c in 1:d){
    data_temp<-data[,,c]
    data_temp[missing]<-NA
    data[,,c]<-data_temp
  }

  # output
  return(list(data = data,theta = theta, psi = psi))

}
