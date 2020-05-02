#' Estimate theta and psi in multinomial mixture model
#' @description This function is generate random sample from Dirichlet distribution
#' @param n - sample size
#' @param alpha - parameters in Dirichlet distribution
#' @return out - generated data
#' @import stats
#' @examples
#' # dimension parameters
#' rdirichlet(n=10,alpha=c(1,1,1))
#' @export


# random sample from Dirichelet distribution
rdirichlet <- function(n=1,alpha=c(1,1)) {
  # Simulations from the Dirichlet Distribution, according to the method of Winipedia
  len_alpha <- length(alpha)
  out <- matrix(0,n,len_alpha)
  gams <- matrix(0,n,len_alpha)
  for (i in 1:len_alpha) {
    gams[,i] <- matrix(rgamma(n,alpha[i]),n,1)
  }
  gamtotal <- matrix(rowSums(gams),n,len_alpha)
  out <- gams/gamtotal
  return(out)
}


