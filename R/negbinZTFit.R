
.logLikNegBinzt <- function(gam, p, values, freq, N){
  val1 <- sum(lgamma(gam + values) * freq)
  val2 <- log(p) * sum(values * freq)
  val3 <- sum(lfactorial(values) * freq)
  val4 <- N * (gam*log(1-p)-log(1-((1-p)^(gam))) - lgamma(gam))
  (val1 + val2 - val3 + val4)
}



.mleNBzt <- function(params, values, freq, N){
  gam <- params[1]
  p <- params[2]
  print(c(gam, p))
  -(.logLikNegBinzt(gam, p, values, freq, N))
}


#' Negative binomial truncated
#' MLE de la Binomial negativa truncada en cero.
#'
#' @export
negbinZTFit <- function(data, init_gamma, init_p, level = 0.95, ...){
  Call <- match.call()

  if(!is.numeric(init_gamma) || !is.numeric(init_p)){
    stop('Wrong intial values for the parameters.')
  }

  tryCatch(
    {
      estResults <- stats::optim(par = c(init_gamma, init_p), .mleNBzt,
                                 values = data[,1], freq = data[,2],
                                 N = sum(data[,2]), hessian = TRUE)
        #.paramEstimationBase(data, c(init_alpha, init_beta), .mloglikelihood, ...)
      estGamma <- as.numeric(estResults$par[1])
      estp <- as.numeric(estResults$par[2])
      paramSD <- sqrt(diag(solve(estResults$hessian)))
      paramsCI <- .getConfidenceIntervals(paramSD, estGamma, estp, level)

      structure(class = "negbinZTR", list(gammaHat = estGamma,
                                         pHat = estp,
                                         gammaSD = paramSD[1],
                                         pSD = paramSD[2],
                                         gammaCI = c(paramsCI[1,1],paramsCI[1,2]),
                                         pCI = c(paramsCI[2,1],paramsCI[2,2]),
                                         logLikelihood = -estResults$value,
                                         hessian = estResults$hessian,
                                         call = Call))
    },
    error=function(cond) {
      print(cond)
      return(NA)
    })
}

#' @rdname negbinZTFit
#' @export
residuals.negbinZTR<- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  fitted.values <- fitted(object)
  residual.values <- as.numeric(dataMatrix[, 2]) - fitted.values
  return(residual.values)
}

#' @rdname negbinZTFit
#' @export
fitted.negbinZTR <- function(object, ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  N <- sum(as.numeric(dataMatrix[, 2]))
  fitted.values <- N*sapply(as.numeric(as.character(dataMatrix[,1])), pmfNBzt, gam = object[['gammaHat']],
                            p = object[['pHat']])
  return(fitted.values)
}

#' @rdname negbinZTFit
#' @export
coef.negbinZTR <- function(object, ...){
  estimation <- matrix(nrow = 2, ncol = 4)
  estimation[1, ] <- c(object[['gammaHat']], object[['gammaSD']], object[['gammaCI']][1], object[['gammaCI']][2])
  estimation[2, ] <- c(object[['pHat']], object[['pSD']], object[['pCI']][1], object[['pCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("gamma", "p")
  estimation
}

#' @rdname negbinZTFit
#' @export
plot.negbinZTR <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  graphics::plot(dataMatrix[,1], dataMatrix[,2], log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting Neg. Bin ZT Distribution", ...)

  graphics::lines(as.numeric(dataMatrix[,1]), fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'Neg Bin ZT Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname negbinZTFit
#' @export
print.negbinZTR <- function(x, ...){
  cat('Call:\n')
  print(x[['call']])
  cat('\n')
  cat('Initial Values:\n')
  cat(sprintf('gamma: %s\n', format(eval(x[['call']]$init_gamma), digits = 3)))
  cat(sprintf('p: %s\n', format(eval(x[['call']]$init_p), digits = 3)))
  cat('\n')
  cat('Coefficients:\n')
  print(coef(x))
  cat('\n')
  cat('Metrics:\n')
  cat(sprintf('Log-likelihood: %s\n', logLik(x)))
  cat(sprintf('AIC: %s\n', AIC(x)))
  cat(sprintf('BIC: %s\n', BIC(x)))
}

#' @rdname negbinZTFit
#' @export
summary.negbinZTR <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname negbinZTFit
#' @export
logLik.negbinZTR <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname negbinZTFit
#' @export
AIC.negbinZTR <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 2)
  return(aic)
}

#' @rdname negbinZTFit
#' @export
BIC.negbinZTR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  bic <- .get_BIC(object[['logLikelihood']], 2, sum(as.numeric(dataMatrix[, 2])))
  return(bic)
}
