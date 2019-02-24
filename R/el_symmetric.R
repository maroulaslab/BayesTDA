#'  Elementary Symmetric Function
#'  @name el_symmetric
#'  @usage el_symmetric(values,k)
#'  @param values: vector of inputs
#'  @param k: index for the function i.e., k in the formula above
#'  @return numerical value of the function
#'  @description Evaluate elementary symmetric function, i.e., e(v,k) = e({x_1,x_2,x_3},1) = (x_1+x_2+x_3), or e({x_1,x_2,x_3},2) = (x_1*x_2+x_2*x_3+x_3*x_1)
#'  or,e({x_1,x_2,x_3},3) = (x_1*x_2*x_3).
#'  @example el_symmetric(c(1,2,3),3)
#'  @details Required package is polynom
#'  @export 
#'  @keyword internal
#'  
el_symmetric = function(values,k){
  if(length(values) == 0){
    return(1)
  } else{
    M = length(values)
    polynomial = as.numeric(poly.calc(values))
    return(((-1)^k)*polynomial[M+1]*polynomial[M-k+1])
  }
}

