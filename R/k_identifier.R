#' Identify the number of components k
#' @description This function is applied when number of components k is unknown and number of trials m > 2k -1. Missing data is allowed in dataset.
#' This function identifies the number of components. Then if needed, parameters could be estimated by function FinMix().
#' @param data - an array with dimension c(n,p,d). data[i,j,c] denotes the number of trials occurs in i-th observation, j-th variable and c-th categories.
#' @param T - number of iterations in Gibbs sampler, default value is 1000. T should be an even number for 'burn-in'. The estimates are computed by the second-half iterations.
#' @note  m - number of trials, i.e., m = sum(data[i,j,]). Here m are identical for all i and j.
#' m > 2k - 1 is a condition that guarantees identifiability (Elmore and Wang, 2003).
#' @note  n - number of observations
#' @note  p - number of variables
#' @note  d - number of categories. Here d are identical for all variables.
#' @return k - number of components
#' @import stats
#' @importFrom DirichletReg rdirichlet
#' @examples
#' ## IncompleteData example
#' data("IncompleteData")
#' k_identifier(data = IncompleteData$data)
#' @references
#' [1] Elmore, Ryan, and Shaoli Wang. "ldentifiability and Estimation in Finite Mixture Models with Multinomial Components."
#' @export

k_identifier<-function(data,T = 1000){

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
  beta<-rep(1,d)     # psi ~ Dirichlet(beta)
  epsilon1<-1        # q ~ beta(epsilon1,epsilon2)
  epsilon0<-1
  a<-0.25            # alpha ~ gamma(a,b)
  b<-0.25

  # initial values
  # initial alpha
  alpha<-rgamma(1,shape = a,rate = b)

  # initial V, theta, k
  V<-c()
  theta<-c()
  h<-1
  V[1]<-rbeta(1,1,alpha)    # V ~ beta(1,alpha)
  theta[1]<-V[1]
  while(sum(theta)<0.9999){
    h<-h+1
    V[h]<-rbeta(1,1,alpha)
    theta[h]<-V[h]*prod(1-V[1:h-1])
  }
  k<-h

  # inital z
  z<-rep(0,n)
  for(i in 1:n){
    z[i]<-which(rmultinom(1,1,theta)==1)
  }

  # initial psi
  psi<-array(0,dim = c(k,p,d))
  for(h in 1:k){
    psi[h,,]<-rdirichlet(p,beta)
  }

  # initial q
  q<-matrix(0, nrow = k, ncol = p)
  for(h in 1:k){
    q[h,]<-rbeta(p,epsilon1,epsilon0)
  }

  # gibbs sampler
  # parameter track
  k_track<-rep(0,T)
  for(t in 1:T){
    # update u
    u<-rep(0,n)
    for(i in 1:n){
      u[i]<-runif(1,0,theta[z[i]])
    }


    # update V,theta,k
    h<-0
    theta<-c()
    while(sum(theta)<1-min(u) & sum(theta)<0.9999){
      h<-h+1
      # lower bound
      lb_num<-u[z==h]
      if(length(lb_num)==0)
        lb_num<-0

      lb_num<-max(lb_num)

      if(h==1){
        lb_dom<-1
      }else{
        lb_dom<-prod(1-V[1:h-1])
      }

      lb<-lb_num/lb_dom

      # upper bound
      pos<-which(z>h)
      ub<-c()
      for(ii in pos){
        ub_num<-u[ii]*(1-V[h])
        ub_dom<-V[z[ii]]*prod(1-V[1:z[ii]-1])

        ub<-c(ub,ub_num/ub_dom)
      }

      if(length(ub)==0)
        ub<-0

      ub<-1-max(ub)

      # inverse sampling
      rnd<-runif(1,0,1)
      V[h]<-1-((1-lb)^alpha-rnd*((1-lb)^alpha-(1-ub)^alpha))^(1/alpha)

      # update theta
      if(h==1){
        theta[h]<-V[h]
      }else{
        theta[h]<-V[h]*prod(1-V[1:h-1])
      }
    }

    k<-h
    V<-V[1:k]

    # update q
    q<-matrix(0, nrow = k, ncol = p)
    for(h in 1:k){
      for(j in 1:p){
        q[h,j]<-rbeta(1,epsilon1+sum(r[z==h,j]),epsilon0+length(r[z==h,j])-sum(r[z==h,j]))
      }
    }

    # update psi
    beta_temp<-array(0,dim = c(k,p,d))
    for(h in 1:k){
      for(j in 1:p){
        for(c in 1:d){
          beta_temp[h,j,c]<-sum(data[z==h,j,c],na.rm = TRUE)
        }
      }
    }

    beta_temp<-beta_temp+beta

    psi<-array(0,dim = c(k,p,d))
    for(h in 1:k){
      for(j in 1:p){
        psi[h,j,]<-rdirichlet(1,beta_temp[h,j,])
      }
    }

    # update z
    z<-rep(0,n)
    for(i in 1:n){
      f<-rep(0,k)
      for(h in 1:k){
        if(u[i]<theta[h]){
          f[h]<-1
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
        }else
          f[h]<-0
      }

      f<-f/sum(f)
      z[i]<-which(rmultinom(1,1,f)==1)

    }

    # update alpha
    V_temp<-V
    V_temp[V==1]<-0.99999
    alpha<-rgamma(1,a+k,b-sum(log(1-V_temp)))

    # parameter track
    k_track[t]<-k

  }

  # parameter estimate
  k_est<-as.numeric(names(which.max(table(k_track[(T/2+1):T]))))

  # output
  return(list(k_track = k_track,k_est = k_est))

}
