#' @export
dzipf <- function(x, alpha){
  (x^(-alpha))/VGAM::zeta(alpha)
}
