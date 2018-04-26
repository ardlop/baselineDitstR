### Negative Binomia Zero Truncated
#' @export
pmfNBzt <- function(x, gam, p){
  (gamma(gam + x) * p^(x) * (1 - p)^gam)/(factorial(x) * gamma(gam) * (1 - (1 - p)^gam))
}
