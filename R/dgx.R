#' pmf Discrete Gaussian Exponential Distribution
#' @export
pmfDGX <- function(x, mu, sig){
  aValue <- .getAConstant(mu, sig)
  ((aValue/x)*exp(-((log(x) - mu)^2)/(2*sig^2)))
}
