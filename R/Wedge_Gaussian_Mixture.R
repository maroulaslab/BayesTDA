#' Restricted Gaussian Mixture density
#' @description Computes mixed Gaussian density restricted to wedge W. Require \code{mvtnorm} package.
#' @usage Wedge_Gaussian_Mixture(x, weights, means, sigmas)
#' @name Wedge_Gaussian_Mixture
#' @param x: two element vector. These are the points to evaluate density at.
#' @param weights: a vector of the mixture weights
#' @param means: a list of two element vector where the vectors are the means of the restricted Gaussian
#' @param sigmas: a vector of positive constants, sigmas in covariance matrices of the restricted Gaussian
#' @return The function \code{Wedge_Gaussian_Mixture} returns numerical value of the mixed Gaussian density restricted to wedge W.
#' @examples # input the mean and constant sigma of covariance matrix
#' w <- Wedge_Gaussian_Mixture(x = c(1,2), weights = 1, means = list(c(1,1)), sigmas= 0.1)
#' @keywords internal
#' @export

Wedge_Gaussian_Mixture = function(x, weights, means, sigmas){
  to_sum = weights*mapply(Wedge_Gaussian,means,sigmas,MoreArgs = list(z = x))
  return(sum(to_sum))
}



