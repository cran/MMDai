#' Finite k components and multiple m trials
#' @description This function is applied when number of components k is known and number of trials m > 2k -1. Missing data is allowed in dataset.
#' This function estimates parameters of latent joint distribution which observations are generated from.
#' @param data - an array with dimension c(n,p,d). data[i,j,c] denotes the number of trials occurs in i-th observation, j-th variable and c-th categories.
#' @param k - number of components is known
#' @param T - number of iterations in Gibbs sampler, default value is 2000. T should be an even number for 'burn-in'. The estimates are computed by the second-half iterations.
#' @note  m - number of trials, i.e., m = sum(data[i,j,]). Here m are identical for all i and j.
#' m > 2k - 1 is a condition that guarantees identifiability (Elmore and Wang, 2003).
#' @note  n - number of observations
#' @note  p - number of variables
#' @note  d - number of categories. Here d are identical for all variables.
#' @note  r - matrix with n rows and p columns. If data[i,j,] is missing, r[i,j] = 0. Otherwise r[i,j] = 1. data[i,j,] must be missing or observed for all trials.
#' @note  q - matrix with k rows and p columns, denotes missing probability for each variable in different components.
#' @return theta - a vector that sum to 1, denotes probability of latent class.
#' @return psi - an array with dimension c(k,p,d), specific probability for each variables in each component.
#' @import stats
#' @importFrom DirichletReg rdirichlet
#' @examples
#' k <- 2
#' ## IncompleteData example
#' data("IncompleteData")
#' FinMix(data = IncompleteData$data,k)
#' @references
#' [1] Elmore, Ryan, and Shaoli Wang. "ldentifiability and Estimation in Finite Mixture Models with Multinomial Components."
#' @export

FinMix<-function(data,k,T = 2000){

  # dimensional parameters
  n<-dim(data)[1]    # number of observations
  p<-dim(data)[2]    # number of variables
  d<-dim(data)[3]    # number of categories

  # missing data position
  r<-matrix(0,nrow = n,ncol = p)
  pos<-which(!is.na(data[,,1]))
  for(c in 1:d){
    if(all(pos==which(!is.na(data[,,c]))))
      next
    else
      return(cat("Error: data[i,j,] must be missing or observed for all trials.\n"))
  }
  r[pos]<-1


  # Bayesian inference
  # hyper parameter
  delta<-rep(1,k)    # theta ~ Dirichlet(delta)
  beta<-rep(1,d)     # psi ~ Dirichlet(beta)
  epsilon1<-1        # q ~ beta(epsilon1,epsilon2)
  epsilon0<-1

  # initial values
  theta<-rdirichlet(1,delta)

  psi<-array(0,dim = c(k,p,d))
  for(h in 1:k){
    psi[h,,]<-rdirichlet(p,beta)
  }

  q<-matrix(0, nrow = k, ncol = p)
  for(h in 1:k){
    q[h,]<-rbeta(p,epsilon1,epsilon0)
  }

  # gibbs sampler
  # parameter track
  theta_track<-matrix(0,nrow = T,ncol = k)
  psi_track<-array(0,dim = c(T,k,p,d))
  q_track<-array(0,dim = c(T,k,p))
  for(t in 1:T){
    # update z
    z<-rep(0,n)   # latent class
    for(i in 1:n){
      f<-rep(0,k)
      for(h in 1:k){
        f[h]<-theta[h]
        for(j in 1:p){
          if(r[i,j]==1){
            f[h]<-f[h]*q[h,j]
            for(c in 1:d){
              f[h]<-f[h]*psi[h,j,c]^data[i,j,c]
            }
          }else{
            f[h]<-f[h]*(1-q[h,j])
          }
        }
      }

      f<-f/sum(f)
      z[i]<-which(rmultinom(1,1,f)==1)

    }

    # update theta
    delta_temp<-rep(0,k)
    for(i in 1:n){
      delta_temp[z[i]]<-delta_temp[z[i]]+1
    }

    theta<-rdirichlet(1,delta+delta_temp)


    # update q
    q<-matrix(0, nrow = k, ncol = p)
    for(h in 1:k){
      for(j in 1:p){
        q[h,j]<-rbeta(1,epsilon1+sum(r[z==h,j]),epsilon0+length(r[z==h,j])-sum(r[z==h,j]))
      }
    }

    # update psi
    beta_temp<-array(0,dim = c(k,p,d))
    for(i in 1:n){
      for(j in 1:p){
        if(r[i,j]==1){
          for(c in 1:d){
            beta_temp[z[i],j,c]<-beta_temp[z[i],j,c]+data[i,j,c]
          }
        }
      }
    }

    beta_temp<-beta_temp+beta

    for(h in 1:k){
      for(j in 1:p){
        psi[h,j,]<-rdirichlet(1,beta_temp[h,j,])
      }
    }

    # reorder to avoid label switching
    order<-rank(theta)

    theta<-theta[order]
    q<-q[order,]
    psi<-psi[order,,]

    # parameter track
    theta_track[t,]<-theta
    psi_track[t,,,]<-psi
    q_track[t,,]<-q

  }

  # parameter estimate
  theta_est<-rep(0,k)
  psi_est<-array(0,dim = c(k,p,d))
  theta_sd<-rep(0,k)
  psi_sd<-array(0,dim = c(k,p,d))

  for(t in (T/2+1):T){
    theta_est<-theta_est+theta_track[t,]
    psi_est<-psi_est+psi_track[t,,,]
  }

  theta_est<-theta_est/(T/2)
  psi_est<-psi_est/(T/2)

  for(t in (T/2+1):T){
    theta_sd<-theta_sd+(theta_track[t,]-theta_est)^2
    psi_sd<-psi_sd+(psi_track[t,,,]-psi_est)^2
  }

  theta_sd<-theta_sd/(T/2)
  psi_sd<-psi_sd/(T/2)

  # output
  return(list(theta = theta_est,psi = psi_est,theta_sd = theta_sd,psi_sd = psi_sd))


}
