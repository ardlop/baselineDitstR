.zeta_Distribution <- function(alpha, nSize, freq, values){
  -( -alpha*sum(freq * log(values)) - nSize * log(VGAM::zeta(alpha)))
}


#' @export
zipfFit <- function(data, init_alpha, level = 0.95, ...){
  Call <- match.call()

  if(!is.numeric(init_alpha)){
    stop('Wrong intial value for the parameter.')
  }

  tryCatch(
    {
      estResults <- stats::optim(par = init_alpha, .zeta_Distribution,
                                 nSize = sum(data[,2]), freq = data[,2],
                                 values = data[,1], hessian = TRUE, method='Brent',
                                 lower = 1, upper = 100000)
      estAlpha <- as.numeric(estResults$par[1])
      paramSD <- sqrt(diag(solve(estResults$hessian)))
      #He puesto 0 para evitar tocar la definicion de la funcion, aquí no hay extra param.
      paramsCI <- .getConfidenceIntervals(paramSD, estAlpha,0,level)

      structure(class = "zipf_R", list(alphaHat = estAlpha,
                                        alphaSD = paramSD[1],
                                        alphaCI = c(paramsCI[1,1], paramsCI[1,2]),
                                        logLikelihood = -estResults$value,
                                        hessian = estResults$hessian,
                                        call = Call))
    },
    error=function(cond) {
      print(cond)
      return(NA)
    })
}

#' @rdname zipfFit
#' @export
residuals.zipf_R<- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  fitted.values <- fitted(object)
  residual.values <- as.numeric(dataMatrix[, 2]) - fitted.values
  return(residual.values)
}

#' @rdname zipfFit
#' @export
fitted.zipf_R <- function(object, ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  N <- sum(as.numeric(dataMatrix[, 2]))
  fitted.values <- N*sapply(as.numeric(as.character(dataMatrix[,1])), dzipf, alpha = object[['alphaHat']])
  return(fitted.values)
}

#' @rdname zipfFit
#' @export
coef.zipf_R <- function(object, ...){
  estimation <- matrix(nrow = 1, ncol = 4)
  estimation[1, ] <- c(object[['alphaHat']], object[['alphaSD']], object[['alphaCI']][1], object[['alphaCI']][2])
  #estimation[2, ] <- c(object[['pHat']], object[['pSD']], object[['pCI']][1], object[['pCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("alpha")
  estimation
}

#' @rdname zipfFit
#' @export
plot.zipf_R <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  graphics::plot(dataMatrix[,1], dataMatrix[,2], log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting Zipf Distribution", ...)

  graphics::lines(as.numeric(dataMatrix[,1]), fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'Zipf Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname zipfFit
#' @export
print.zipf_R <- function(x, ...){
  cat('Call:\n')
  print(x[['call']])
  cat('\n')
  cat('Initial Values:\n')
  cat(sprintf('alpha: %s\n', format(eval(x[['call']]$init_alpha), digits = 3)))
  #cat(sprintf('p: %s\n', format(eval(x[['call']]$init_p), digits = 3)))
  cat('\n')
  cat('Coefficients:\n')
  print(coef(x))
  cat('\n')
  cat('Metrics:\n')
  cat(sprintf('Log-likelihood: %s\n', logLik(x)))
  cat(sprintf('AIC: %s\n', AIC(x)))
  cat(sprintf('BIC: %s\n', BIC(x)))
}

#' @rdname zipfFit
#' @export
summary.zipf_R <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname zipfFit
#' @export
logLik.zipf_R <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname zipfFit
#' @export
AIC.zipf_R <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 1)
  return(aic)
}

#' @rdname zipfFit
#' @export
BIC.zipf_R <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  bic <- .get_BIC(object[['logLikelihood']], 1, sum(as.numeric(dataMatrix[, 2])))
  return(bic)
}
