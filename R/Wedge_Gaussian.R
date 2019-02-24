#' Restricted Gaussian density
#' The function \code{Wedge_Gaussian} computes Gaussian density restricted to wedge W.
#' @name Wedge_Gaussian
#' @usage Wedge_Gaussian(z, mean, sigma)
#' @param z: two element vector. These are the points to evaluate density at.
#' @param mean: two element vector, mean of the restricted Gaussian
#' @param sigma: positive constant. Covariance matrix of the restricted Gaussian has the form sigma*I_2, where I_2 is the 2x2 identity matrix.
#' @details The function \code{Wedge_Gaussian} computes Gaussian density restricted to wedge W. 
#' For a death vs. birth persistence diagram wedge is defined as the two dimensional coordinate (b,d) such that d>b>0. 
#' For a tilted representation, i.e., persistence vs. birth persistence diagram wedge is defined as T(b,d) = (b, d-b).
#' @return The function \code{Wedge_Gaussian} returns numerical value of the  Gaussian density restricted to wedge W.
#' @description The function \code{Wedge_Gaussian} computes Gaussian density restricted to wedge W. Require \code{mvtnorm} package for computing multivariate normal pdf and cdf.
#' @examples # input the mean and constant sigma of covariance matrix
#' w <- Wedge_Gaussian(z = c(1,2), mean = c(1,1), sigma = 0.1)
#' @keywords internal
#' @export

Wedge_Gaussian = function(z,mean,sigma){
  num = dmvnorm(z,mean = mean, sigma = sigma*diag(2))
  return(num)
}


