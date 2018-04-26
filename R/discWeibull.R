# Discrete Weibull truncated
#' @export
pmfDiscreteWeibull <- function(p, v, k){
  ((p^(k^(v)) - p^((k+1)^(v)))/p)
}

