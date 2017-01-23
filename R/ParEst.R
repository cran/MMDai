#' Estimate theta and psi in multinomial mixture model
#' @description This function is used to estimate theta and psi in multinomial mixture model given the number of components k.
#' @param data - data in matrix formation with n rows and p columns
#' @param d - number of categories for each variable
#' @param k - number of components
#' @param T - number of iterations in Gibbs sampler, default value is 1000. T should be an even number for 'burn-in'.
#' @return theta - vector of probability for each component
#' @return psi - specific probability for each variable in each component
#' @import stats
#' @importFrom DirichletReg rdirichlet
#' @export

ParEst<-function(data, d, k, T = 1000){

  # dimensional parameters
  p<-ncol(data)    # number of variables
  n<-nrow(data)    # number of observations

  # missing data
  d<-d+1           # including missing category
  for(j in 1:p){
    data[is.na(data[,j]),j]<-d[j]
  }

  # Bayesian inference
  # initial theta
  theta<-rdirichlet(1,rep(1,k))

  # initial psi
  psi<-list()
  for(j in 1:p){
    psi[[j]]<-rdirichlet(k,rep(1,d[j]))
  }

  # initial z
  z<-rep(0,n)
  for(i in 1:n){
    z[i]<-which(rmultinom(1,1,theta)==1)
  }

  # gibbs sampler
  # initial log posterior probability
  pp<-0
  for(i in 1:n){
    for(j in 1:p){
      pp<-pp+log(psi[[j]][z[i],data[i,j]])
    }
  }

  # posterior probability track
  pp_track<-rep(pp,T)

  # initial estimator
  theta_est<-theta
  psi_est<-psi
  for(t in 2:T){
    # update z
    z<-rep(0,n)   # latent class
    for(i in 1:n){
      f<-rep(0,k)
      for(h in 1:k){
        f[h]<-theta[h]
        for(j in 1:p){
          f[h]<-f[h]*(psi[[j]][h,data[i,j]])
        }
      }

      f<-f/sum(f)
      z[i]<-which(rmultinom(1,1,f)==1)
    }

    # update theta
    delta<-rep(0,k)
    for(h in 1:k){
      delta[h]<-length(z[z==h])
    }

    theta<-rdirichlet(1,delta+rep(1,k))

    # update psi
    for(j in 1:p){
      for(h in 1:k){
        subdata<-data[z==h,j]
        beta<-rep(1,d[j])
        for(c in 1:d[j]){
          beta[c]<-beta[c]+length(subdata[subdata==c])
        }
        psi[[j]][h,]<-rdirichlet(1,beta)
      }
    }

    # reorder to avoid label switching
    if(k>1){
      order<-rank(theta)

      theta<-theta[order]
      for(j in 1:p){
        psi[[j]]<-psi[[j]][order,]
      }
    }

    # posterior probability
    pp<-0
    for(i in 1:n){
      for(j in 1:p){
        pp<-pp+log(psi[[j]][z[i],data[i,j]])
      }
    }

    pp_track[t]<-pp

    # parameter estimation
    if(pp==max(pp_track)){
      theta_est<-theta
      psi_est<-psi
    }
  }

  # renormalize
  for(j in 1:p){
    for(h in 1:k){
      psi_est[[j]][h,1:(d[j]-1)]<-psi_est[[j]][h,1:(d[j]-1)]/(1-psi_est[[j]][h,d[j]])
    }

    psi_est[[j]]<-psi_est[[j]][,1:(d[j]-1)]
  }

  # avoid reduce to vector
  if(k==1){
    for(j in 1:p){
      psi_est[[j]]<-matrix(psi_est[[j]],nrow = 1)
    }
  }

  # output
  return(list(theta = theta_est,psi = psi_est))

}
