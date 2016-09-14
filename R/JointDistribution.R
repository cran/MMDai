#' Joint distribution
#' @description This function is used to compute joint distribution given the estimates of theta and psi.
#' @param theta - a vector that sum to 1, denotes probability of latent class.
#' @param psi - an array with dimension c(k,p,d), specific probability for each variables in each component.
#' @param input - a p-dimensional vector denote specific variable, e.g., assume p=2, input = (1,2) means the joint probability of V_1=1 and V_2=2.
#' @note  k - number of components
#' @note  n - number of observations
#' @note  p - number of variables
#' @note  d - number of choices. Here d are identical for all variables.
#' @importFrom DirichletReg rdirichlet
#' @examples
#' k<-2; p<-2; d<-2;
#' theta<-c(0.3,0.7)
#' psi <- array(rdirichlet(k*p,rep(1,d)),dim = c(k,p,d))
#' JointDistribution(theta,psi,input = c(1,1))
#' @export


JointDistribution<-function(theta,psi,input){
  # dimensional parameters
  k<-dim(psi)[1]    # number of components
  p<-dim(psi)[2]    # number of variables
  d<-dim(psi)[3]    # number of choices

  # joint probabilility for input
  JointProb<-0
  for(h in 1:k){
    # probability in each component
    prob_temp<-1
    for(j in 1:p){
      prob_temp<-prob_temp*psi[h,j,input[j]]
    }

    JointProb<-JointProb+theta[h]*prob_temp
  }

  return(list(JointProb = JointProb))

}
