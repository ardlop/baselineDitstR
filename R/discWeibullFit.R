.mleDisWeib <- function(params, values, freq){
  p <- params[1]
  v <- params[2]

  probs <- sapply(values, function(i, p, v){
    pmfDiscreteWeibull(p, v, i)
  }, p = p, v = v)
  -(sum(log(probs) * freq))
}

#' @export
discWeibullZTFit <- function(data, init_p, init_v, level = 0.95, ...){
  Call <- match.call()

  if(!is.numeric(init_p) || !is.numeric(init_v)){
    stop('Wrong intial values for the parameters.')
  }

  tryCatch(
    {
      estResults <- stats::optim(par = c(init_p, init_v), .mleDisWeib,
                                 values = data[,1], freq = data[,2],
                                 hessian = TRUE)
      #.paramEstimationBase(data, c(init_alpha, init_beta), .mloglikelihood, ...)
      estp <- as.numeric(estResults$par[1])
      estv <- as.numeric(estResults$par[2])
      paramSD <- sqrt(diag(solve(estResults$hessian)))
      paramsCI <- .getConfidenceIntervals(paramSD, estp, estv, level)

      structure(class = "dwZTR", list(pHat = estp,
                                          vHat = estv,
                                          pSD = paramSD[1],
                                          vSD = paramSD[2],
                                          pCI = c(paramsCI[1,1],paramsCI[1,2]),
                                          vCI = c(paramsCI[2,1],paramsCI[2,2]),
                                          logLikelihood = -estResults$value,
                                          hessian = estResults$hessian,
                                          call = Call))
    },
    error=function(cond) {
      print(cond)
      return(NA)
    })
}

#' @rdname discWeibullFit
#' @export
residuals.dwZTR<- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  fitted.values <- fitted(object)
  residual.values <- as.numeric(dataMatrix[, 2]) - fitted.values
  return(residual.values)
}

#' @rdname discWeibullFit
#' @export
fitted.dwZTR <- function(object, ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  N <- sum(as.numeric(dataMatrix[, 2]))
  fitted.values <- N*sapply(as.numeric(as.character(dataMatrix[,1])), pmfDiscreteWeibull, p = object[['pHat']],
                            v = object[['vHat']])
  return(fitted.values)
}

#' @rdname discWeibullFit
#' @export
coef.dwZTR <- function(object, ...){
  estimation <- matrix(nrow = 2, ncol = 4)
  estimation[1, ] <- c(object[['pHat']], object[['pSD']], object[['pCI']][1], object[['pCI']][2])
  estimation[2, ] <- c(object[['vHat']], object[['vSD']], object[['vCI']][1], object[['vCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("p", "v")
  estimation
}

#' @rdname discWeibullFit
#' @export
plot.dwZTR <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  graphics::plot(dataMatrix[,1], dataMatrix[,2], log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting DW ZT Distribution", ...)

  graphics::lines(as.numeric(dataMatrix[,1]), fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'DW ZT Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname discWeibullFit
#' @export
print.dwZTR <- function(x, ...){
  cat('Call:\n')
  print(x[['call']])
  cat('\n')
  cat('Initial Values:\n')
  cat(sprintf('p: %s\n', format(eval(x[['call']]$init_p), digits = 3)))
  cat(sprintf('v: %s\n', format(eval(x[['call']]$init_v), digits = 3)))
  cat('\n')
  cat('Coefficients:\n')
  print(coef(x))
  cat('\n')
  cat('Metrics:\n')
  cat(sprintf('Log-likelihood: %s\n', logLik(x)))
  cat(sprintf('AIC: %s\n', AIC(x)))
  cat(sprintf('BIC: %s\n', BIC(x)))
}

#' @rdname discWeibullFit
#' @export
summary.dwZTR <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname discWeibullFit
#' @export
logLik.dwZTR <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname discWeibullFit
#' @export
AIC.dwZTR <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 2)
  return(aic)
}

#' @rdname discWeibullFit
#' @export
BIC.dwZTR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  bic <- .get_BIC(object[['logLikelihood']], 2, sum(as.numeric(dataMatrix[, 2])))
  return(bic)
}
