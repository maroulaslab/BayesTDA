#' Posterior cardinality for persistence diagrams modeled as IID cluster Point Process
#' @description This function use the Bayesian inference obtained by characterizing persistence diagrams (PDs) as IID cluster point process. 
#' Posterior cardinality can be computed from a prior distribution and a set of observed persistence diagrams. 
#' The inputs consists of two intensity functions estimated by using \code{Wedge_Gaussian_Mixture} function: one is the prior intensity and another one is denoted as unexpected.
#' And two cardinality functions defined as Binomials.
#' The prior distribution will be depending on the knowledge we have about the underlying truth. The observed PDs can exhibit two features.
#' One which is believed to be associated to underlying truth and one is not anticipated by underlying truth or produce due to the noise. We call the later one unexpected features.
#' Usually these features are considered to occur near the birth axis as they can be formed due to noise.
#' @name postCardiid
#' @usage postCardiid(n,Dy,alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.Dyo,prob.unexpected,weights.unexpected,mean.unexpected,sigma.unexpected,Nmax)
#' @param n: positive integer. The posterior cardinalities are computed at n and 0<= n<= Nmax.
#' @param Dy: list of two element vectors representing points observed in a persistence diagram. Here we consider a fixed homological feature.
#' The coordinates (birth, death)/(birth, persistence) need to be a list of two element vectors. 
#' @param alpha: 0<=alpha<1. The probability of a feature in the prior will be detected in the observation. 
#' @param prob.prior: The prior cardinality is defined as a binomial and prob.prior is the probability term, i.e., prior cardinality ~ B(Nmax,prob.prior). 
#' @param weight.prior: a list of mixture weights for the prior density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate prior density.
#' @param mean.prior: a list of two element vector, means of the prior density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate prior density.
#' @param sigma.prior: a list of positive constants, sigmas in covariance matrices of the prior density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate prior density.  
#' @param sigma.Dyo: positive constant. variance coefficient of the likelihood density. This represents the degree of faith on the observed PDs representing the underlying truth.
#' @param prob.unexpected: The unexpected cardinality is defined as a binomial and prob.unexpected is the probability term, i.e., unexpected cardinality ~ B(Nmax,unexpected). 
#' @param weights.unexpected: a list of mixture weights for the unexpected density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate unexpected density.
#' @param mean.unexpected: a list of two element vector, means of the unexpected density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate unexpected density.
#' @param sigma.unexpected: a list of positive constants, sigmas in covariance matrices of the unexpected density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate unexpected density.  
#' @param Nmax: The maximum number of points at which the posterior cardinality will be truncated, i.e., {p_n} is computed for n=0 to n=Nmax. Also, this Nmax will be used to define prior and unexpected cardinality. 
#' So, it should be large enough so that all Binomials involved in the computation make sense
#' @return The function \code{postCardiid} returns posterior cardinality given prior and set of observed PDs using Bayesian framework, where PDs are characterized by IID cluster point process.
#' @details Required packages are \code{mvtnorm}, \code{polynom}, \code{purrr} and \code{TDA}
#' @references Bayesian Inference for Persistent Homology, V Maroulas, F Nasrin, C Oballe, \url{https://arxiv.org/abs/1901.02034}
#' @example 
#' # sample data created from a unit circle to define prior
#' # sample data created from a unit circle to define prior
#' set.seed(88)
#' t = seq(from=0, to = 1, by = 0.01)
#' x = cos(2*pi*t)
#' y = sin(2*pi*t)
#' coord.mat = cbind(x,y)
#' 
#' # persistence diagram using ripsDiag of TDA package.  
#' #library(TDA)
#' circle.samp = coord.mat[sample(1:nrow(coord.mat),50),]
#' pd.circle = ripsDiag(circle.samp,maxscale=3,maxdimension = 1)
#' circle.diag = matrix(pd.circle$diagram,ncol=3)
#' circle.holes = matrix(circle.diag[which(circle.diag[,1]==1),],ncol=3)
#' # convert the persistence diagram to a tilted persistence diagram.
#' circle.tilted = cbind(circle.holes[,2],circle.holes[,3]-circle.holes[,2])
#' circle.df = data.frame(Birth = circle.tilted[,1], Persistence = circle.tilted[,2]) ## create the dataframe consisting of pd coordinates.

#' # these parameters are required to define the object unexpected
#' unex.mean = list(c(0.5,0))
#' unex.noise = 1
#' unex.weight = 1
#' unex.card = length(unex.mean)
#' # these parameters are required to define the object prior
#' inf.mean = list(c(0.5,1.2))
#' inf.noise = 0.2
#' inf.weight = 1
#' inf.card = length(inf.mean)
#' ## next we define the likelihood from the observed pds. Here we consider a dataset from the same unit circle that is perturbed by a gaussian noise.
#' set.seed(7)
#' sy = 0.01 #the variance of the likelihood density function
#' circle.noise = 0.001 ## the variance coefficient of noise in this dataset 
#' noise = rmvnorm(nrow(coord.mat),mean = c(0,0), sigma = circle.noise*diag(2))
#' coord.mat = coord.mat + noise ## gives us the dataset which we will consider as observation
#' 
#' 
#' ## following section creates observed PD
#' circle.samp = coord.mat[sample(1:nrow(coord.mat),50),]
#' pd.circle = ripsDiag(circle.samp,maxscale=3,maxdimension = 1)
#' circle.diag = matrix(pd.circle$diagram,ncol=3)
#' circle.holes = matrix(circle.diag[which(circle.diag[,1]==1),],ncol=3)
#' circle.tilted = cbind(circle.holes[,2],circle.holes[,3]-circle.holes[,2])
#' circle.df = data.frame(Birth = circle.tilted[,1], Persistence = circle.tilted[,2]) 
#' dy = lapply(1:nrow(circle.df),function(x){unlist(c(circle.df[x,][1],circle.df[x,][2]))}) ## creates the parameter Dy
#' 
#' ##plot the observed data
#' plot(x,y, type='l',col = 'red', lwd = 1.5, asp = 1)
#' points(circle.samp[,1],circle.samp[,2],pch = 20, col = 'blue')
#' 
#' ## next we compute the posterior
#' alpha = 0.99 #probablity of a feature in prior to be detected in observations
#' Nmax = 7
#' sig_Dyo = 0.01
#' prob_Dx = inf.card/Nmax
#'    
#'    ## plot the prior cardinality
#'    card.inf.prior = as.data.frame(cbind(0:Nmax,unlist(lapply(0:Nmax,dbinom,size = Nmax,prob=prob_Dx))))
#'    g.card.inf.prior = ggplot(data = card.inf.prior,aes(V1,V2))+geom_bar(stat="identity")+theme_classic()+xlab("number of points")+ylab("Prior Cardinality")+ggtitle("Informative")
#'    print(g.card.inf.prior)
#'    
#'    ##compute the posterior cardinality
#'    card.inf = as.data.frame(cbind(0:Nmax,unlist(lapply(0:Nmax,postCardiid,Dy = dy,alpha = alpha,prob.prior = prob_Dx,weight.prior = inf.weight,
#'                                                     mean.prior=inf.mean,sigma.prior=inf.noise,sigma.Dyo = sig_Dyo,weights.unexpected=unex.weight,mean.unexpected=unex.mean,sigma.unexpected=unex.noise,Nmax=Nmax))))
#'                                                     
#'    ## plot the cardinality distribution as a pmf
#'    g.card.inf = ggplot(data = card.inf,aes(V1,V2))+geom_bar(stat="identity")+theme_classic()+xlab("number of points")+ylab("Posterior Cardinality")
#'    print(g.card.inf)
#' @export



postCardiid = function(n,Dy,alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.Dyo,weights.unexpected,mean.unexpected,sigma.unexpected,Nmax){
  
  K = length(Dy) 
  prob.unexpected =1/Nmax
  pdys = unlist(lapply(0:K,dbinom,size = K,prob=prob.unexpected)) 
  qy = lapply(Dy,function(y){mapply(function(mu,sig1,sig2){dmvnorm(y,mean=mu,
                                                                   sigma=(sig1+sig2)*diag(2))},mean.prior,sigma.prior,MoreArgs = list(sig2 = sigma.Dyo))}) 
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
  
  top_sum = integer(K+1)
  for(k in 0:K){
    top_sum[k+1] =  factorial(K-k)*perm_coeff(n,k)*pdys[K-k+1]*
      (tilde_c_Dx^(n-k))*el_sym[k+1]
  }
  top = dbinom(n,Nmax,prob.prior)*sum(top_sum)
  
  bottom = integer(Nmax+1)
  for(tau in 0:Nmax){
    bot_sum = integer(min(K,tau)+1)
    for(k in 0:min(K,tau)){
      bot_sum[k+1] =  factorial(K-k)*perm_coeff(tau,k)*pdys[K-k+1]*
        (tilde_c_Dx^(tau-k))*el_sym[k+1]
    }
    bottom[tau+1] = dbinom(tau,Nmax,prob.prior)*sum(bot_sum)
  }
  return(top/sum(bottom))
}



