#' Posterior intensity for persistence diagrams modeled as i.i.d. cluster point process
#' @description This function uses the Bayesian inference obtained by characterizing tilted persistence diagrams (PDs) as i.i.d. cluster point processes. 
#' A tilted persistence diagram is produced by the mapping \eqn{T(b,d)=(b,d-b)}. A i.i.d. cluster point process is a generalization of Poisson point process and characterized 
#' by the cardinality and intensity  simultaneously. Similarly as Poisson point process the intensity reflects the expected number of elements and serves as an analog to the first order moment for a random variable. 
#' Cardinality takes the form of probability mass function(pmf) of the number of points in the point process. Hence i.i.d point process characterization of any object needs to integrate the intensity and cardinality densities. 
#' To that end, we define the prior object of the Bayesian model as i.i.d cluster point process with intensity as Gaussian mixture:
#' \deqn{\lambda_X(x) = \sum_{j = 1}^{N}c^X_{j}N*(x;\mu^X_{j},\sigma^X_{j}I)---------------------(1)}
#' and cardinality pmf as \deqn{p_X(n) = B(N_0,n)p_x^n (1-p_x)^{N_0-n}---------------------(2)}
#' where  N* is the restricted Gaussian on the wedge W where the tilted PD is defined and B(N_0,n) is the binomial coefficient. 
#' The weights \eqn{c^X_j}, means \eqn{\mu^X_j} and variance \eqn{\sigma^X_j} of the Gaussian mixture and the probability \eqn{p_x} depend on the prior knowledge about a PD.  Posterior intensity is computed from a prior intensity and 
#' a set of observed persistence diagrams. The observed PDs can exhibit two features: expected and unexpected features. Features which are believed to be associated to prior are called expected feature and others are unexpected features. 
#' The unexpected features of the observed PDs are also defined similarly as prior having the intensity as Gaussian mixture form:
#' \deqn{\lambda_Y(y) = \sum_{i = 1}^{M}c^Y_{i}N*(y;\mu^Y_{i},\sigma^Y_{i}I)---------------------(3)}
#' and cardinality pmf as \deqn{p_Y(n) = B(M_0,n)p_y^n (1-p_y)^{M_0-n}.---------------------(4)}
#' Their respective weights \eqn{c^Y_i}, means \eqn{\mu^Y_i} and variance \eqn{\sigma^Y_i} and probabilty \eqn{p_y} need to be predefined. 
#' The expected features formed the likelihood density \eqn{l(y|x)=N*(y;x,\sigma I)}, where the variance \eqn{\sigma}is a preselected parameter that reflect the level of faith on the observed PDs.
#' @name postIntensityiid
#' @usage postIntensityiid(x,Dy,alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.y,prob.unexpected,weights.unexpected,mean.unexpected,sigma.unexpected,Nmax)
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
#' @return The function \code{postIntensityiid} returns posterior intensity given prior and set of observed PDs using Bayesian framework, where PDs are characterized by IID cluster point process.
#' @details Required packages are \code{mvtnorm}, \code{polynom}, \code{purrr} and \code{TDA}
#' @references Bayesian Inference for Persistent Homology, V Maroulas, F Nasrin, C Oballe, \url{https://arxiv.org/abs/1901.02034}
#' @example # sample data created from a unit circle to define prior
#' set.seed(88)
#' t = seq(from=0, to = 1, by = 0.01)
#' x = cos(2*pi*t)
#' y = sin(2*pi*t)
#' coord.mat = cbind(x,y)
#' # persistence diagram using ripsDiag of TDA package.  
#' #library(TDA)
#' circle.samp = coord.mat[sample(1:nrow(coord.mat),50),]
#' pd.circle = ripsDiag(circle.samp,maxscale=3,maxdimension = 1)
#' circle.diag = matrix(pd.circle$diagram,ncol=3)
#' circle.holes = matrix(circle.diag[which(circle.diag[,1]==1),],ncol=3)
#' # convert the persistence diagram to a tilted persistence diagram.
#' circle.tilted = cbind(circle.holes[,2],circle.holes[,3]-circle.holes[,2])
#' circle.df = data.frame(Birth = circle.tilted[,1], Persistence = circle.tilted[,2]) ## create the dataframe consisting of pd coordinates.
#' 
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
#' ## plot the prior intensity on a grid. Required package mvtnorm and ggplot2
#' #library(mvtnorm)
#' #library(ggplot2)
#' values.x = seq(from = 0, to = 3, by = 0.1)
#' values.y = seq(from = 0, to = 4, by = 0.1)
#' grid = expand.grid(values.x,values.y)
#' inf.dens = apply(grid,1,Wedge_Gaussian_Mixture,inf.weight,inf.mean,inf.noise)
#' max_inf.dens = max(inf.dens)
#' normalized_inf.dens = inf.dens/max_inf.dens
#' inf.df = data.frame(Birth = grid[,1],Persistence = grid[,2],Intensity = normalized_inf.dens)
#' g.inf = ggplot(data = inf.df, aes(x = Birth, y = Persistence)) + geom_tile(aes(fill=Intensity)) + coord_equal()
#' g.inf = g.inf + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow",midpoint = 0.5, limit = c(0,1))
#' g.inf = g.inf + labs(title='Informative Prior',x='Birth',y='Persistence',fill='')  + theme_grey(base_size=16) +theme(plot.title = element_text(hjust = 0.5))
#' print(g.inf)
#' 
#' ## next we define the likelihood from the observed pds. Here we consider a dataset from the same unit circle that is perturbed by a Gaussian noise.
#' set.seed(7)
#' sy = 0.01 #the variance of the likelihood density function
#' circle.noise = 0.001 ## the variance coefficient of noise in this dataset 
#' noise = rmvnorm(nrow(coord.mat),mean = c(0,0), sigma = circle.noise*diag(2))
#' coord.mat = coord.mat + noise ## gives us the dataset which we will consider as observation
#' 
#'  ## following section creates observed PD
#' circle.samp = coord.mat[sample(1:nrow(coord.mat),50),]
#' pd.circle = ripsDiag(circle.samp,maxscale=3,maxdimension = 1)
#' circle.diag = matrix(pd.circle$diagram,ncol=3)
#' circle.holes = matrix(circle.diag[which(circle.diag[,1]==1),],ncol=3)
#' circle.tilted = cbind(circle.holes[,2],circle.holes[,3]-circle.holes[,2])
#' circle.df = data.frame(Birth = circle.tilted[,1], Persistence = circle.tilted[,2]) 
#' dy = lapply(1:nrow(circle.df),function(x){unlist(c(circle.df[x,][1],circle.df[x,][2]))}) ## creates the parameter Dy
#' plot(x,y, type='l',col = 'red', lwd = 1.5, asp = 1)
#' points(circle.samp[,1],circle.samp[,2],pch = 20, col = 'blue')
#'  ## next we compute the posterior over a grid
#'  alpha = 0.99 #probablity of a feature in prior to be detected in observations
#'  Nmax = 7
#'  sig_Dyo = 0.01
#'  prob_Dx = inf.card/Nmax
#' 
#' intens.inf = apply(grid,1,postIntensityiid,dy,alpha,prob_Dx,inf.weight,
#'                    inf.mean,inf.noise,sig_Dyo,unex.weight,unex.mean,unex.noise,Nmax)
#' 
#' max_intens.inf = max(intens.inf)
#' normalized_intens.inf = intens.inf/max_intens.inf
#' post.df = data.frame(Birth = grid[,1],Persistence = grid[,2], Intensity = normalized_intens.inf)
#' g.post.inf = ggplot(data = post.df, aes(x = Birth, y = Persistence)) + geom_tile(aes(fill=Intensity)) + coord_equal()+
#'   scale_fill_gradient2(low = "blue", high = "red", mid = "yellow",midpoint = 0.5, limit = c(0,1))+
#'   geom_point(data = circle.df, aes(x = Birth, y= Persistence), size=2, color = "Green")+
#'   labs(title='Posterior Using Inf.',x='Birth',y='Persistence', fill = "") + theme_grey(base_size=16) +  theme(plot.title = element_text(hjust = 0.5))
#' print(g.post.inf)
#' @export

postIntensityiid = function(x,Dy,alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.y,weights.unexpected,mean.unexpected,sigma.unexpected,Nmax){
  
  if( alpha ==  1)
    stop('alpha can not be 1')
  else
    
    
    
  first_term = (1-alpha)*Wedge_Gaussian_Mixture(x,weight.prior,mean.prior,sigma.prior)*
    Gamma(Dy,alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.y,weights.unexpected,mean.unexpected,sigma.unexpected,Nmax,1)*
    (1/Gamma(Dy,alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.y,weights.unexpected,mean.unexpected,sigma.unexpected,Nmax,0))
  qy = lapply(Dy,function(y){mapply(function(mu,sig1,sig2){dmvnorm(y,mean=mu,
                                                                   sigma=(sig1+sig2)*diag(2))},mean.prior,sigma.prior,MoreArgs = list(sig2 = sigma.y))}) 
  K = length(Dy) 
  prob.unexpected = 1/Nmax
  g = Gamma(Dy,alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.y,weights.unexpected,mean.unexpected,sigma.unexpected,Nmax,0)
  to_sum = matrix(0,nrow = K,length(weight.prior)) 
  for(j in 1:K){
    for(i in 1:length(weight.prior)){
      to_sum[j,i]=Gamma(Dy[-j],alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.y,weights.unexpected,mean.unexpected,sigma.unexpected,Nmax,1)*
        alpha*qy[[j]][i]*weight.prior[i]*(1/Wedge_Gaussian_Mixture(Dy[[j]],weights = weights.unexpected,means = mean.unexpected,sigmas = sigma.unexpected))*
        (1/g)*
        dmvnorm(x,mean = (sigma.prior[i]*Dy[[j]]+sigma.y*mean.prior[[i]])/(sigma.prior[[i]]+sigma.y),
                sigma = (sigma.y*sigma.prior[[i]]/(sigma.prior[[i]]+sigma.y))*diag(2))
    }
  }
  return(first_term+sum(to_sum))
}



