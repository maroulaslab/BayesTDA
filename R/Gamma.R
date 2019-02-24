#' Gamma
#' @name Gamma
#' @param x: a 2 by 1 vector where the posterior intensity is computed. Cosequently, x is a coordinate in the wedge where the tilted PD is defined.
#' @param Dy: a list of n vectors(2 by 1) representing points observed in a tilted persistence diagram of a fixed homological feature.
#' @param alpha: 0<=alpha<1. The probability of a feature in the prior will be detected in the observation. 
#' @param prob.prior: The prior cardinality pmf is defined as a binomial and prob.prior is the probability term (\eqn{p_x} in Eqn.(2)).
#' @param weight.prior: a N by 1 vector of mixture weights (\eqn{c^X_{j}} of Eqn (1)) for the prior density estimation.
#' @param mean.prior: a list of N vectors(2 by 1) each represets mean of the prior density (\eqn{\mu^X_{j}} of Eqn (1))
#' @param sigma.prior: a N by 1 vector of positive constants, \eqn{\sigma^X_{j}} of Eqn (1).
#' @param sigma.y: a positive constant. Variance coefficient (\eqn{\sigma}) of the likelihood density \eqn{l(y|x)} defined in the description above. This represents the degree of faith on the observed PDs representing the prior.
#' @param weights.unexpected: a M by 1 vector of mixture weights for the unexpected features. i.e., (\eqn{c^Y_{i}} of Eqn (2) above.
#' @param mean.unexpected: a list of M vectors (2 by 1),each represets mean of the Gaussian mixture density (\eqn{\mu^Y_{i}} of Eqn (2)) for the unexpected features.
#' @param sigma.unexpected: a M by 1 vector of positive constants, \eqn{\sigma^Y_{i}} of Eqn (2).
#' @param Nmax: The maximum number of points on which the posterior cardinality will be truncated, i.e., {p_n} is computed for n=0 to n=Nmax. Also, this Nmax will be used to define prior cardinality too. 
#' So, it should be large enough so that all Binomials involved in the computation make sense
#' @param b: 0 or 1
#' @keywords internal
#' @return numeric
#' @export

Gamma = function(Dy,alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.y,weights.unexpected,mean.unexpected,sigma.unexpected,Nmax,b){
  
  if( alpha ==  1)
    stop('alpha can not be 1')
  else
    
    
  K = length(Dy) 
  prob.unexpected =1/Nmax
  pdys = unlist(lapply(0:K,dbinom,size = K,prob=prob.unexpected)) 
  qy = lapply(Dy,function(y){mapply(function(mu,sig1,sig2){dmvnorm(y,mean=mu,
                                                                   sigma=(sig1+sig2)*diag(2))},mean.prior,sigma.prior,MoreArgs = list(sig2 = sigma.y))}) 
  tilde_c = function(means,sigmas,weights){pmvnorm(lower = c(0,0),upper = c(Inf,Inf),mean=means,
                                                   sigma=sigmas*diag(2))[1]*weights} 
  tilde_c_Dx = (1-alpha)*sum(mapply(tilde_c,mean.prior,sigma.prior,weight.prior))
  f_Dys = unlist(lapply(Dy,Wedge_Gaussian_Mixture,weights = weights.unexpected,means = mean.unexpected,sigmas = sigma.unexpected))
  input_el_sym = mapply(function(q,fDys){alpha*t(weight.prior)%*%q*(1/fDys)},qy,f_Dys)
  el_sym = c(1,unlist(lapply(1:K,el_symmetric,values = input_el_sym)))
  
  perm_coeff = function(n,j){
    if(n>=j){
      return(factorial(n)/factorial(n-j))
    }else{return(0)}
  } 
  
  bottom = integer(Nmax+1)
  for(tau in 0:Nmax){
    bot_sum = integer(min(K,tau)+1)
    for(k in 0:min(K,tau)){
      bot_sum[k+1] =  factorial(K-k)*perm_coeff(tau,k+b)*pdys[K-k+1]*
        (tilde_c_Dx^(tau-k-b))*el_sym[k+1]
    }
    bottom[tau+1] = dbinom(tau,Nmax,prob.prior)*sum(bot_sum)
  }
  return(sum(bottom))
}