# dgx distribution
.getAConstant <- function(mu, sig, infSize = 1000){
  probs <- sapply(1:infSize, function(i, mu, sig){
    (1/i)*exp(-((log(i) - mu)^2)/(2*sig^2))
  }, mu = mu, sig = sig)
  (sum(probs)^(-1))
}

.logLikDGX <- function(mu, sig, values, freq, N){
  aValue <- .getAConstant(mu, sig)
  val1 <- sum((log(values) + ((log(values) - mu)^2/(2*sig^2)))*freq)
  (N * log(aValue) - val1)
}

.mleDGX <- function(params, values, freq, N){
  mu <- params[1]
  sig <- params[2]
  # print(c(mu, sig))
  -(.logLikDGX(mu, sig, values, freq, N))
}


#' mle DGX
#' @export
dgxFit <- function(data, init_mu, init_sig, level = 0.95, ...){
  Call <- match.call()

  if(!is.numeric(init_mu) || !is.numeric(init_sig)){
    stop('Wrong intial values for the parameters.')
  }

  tryCatch(
    {
      estResults <- stats::optim(par = c(init_mu, init_sig), .mleDGX,
                                 values = data[,1], freq = data[,2],
                                 N = sum(data[,2]), hessian = TRUE)
      estmu <- as.numeric(estResults$par[1])
      estsig <- as.numeric(estResults$par[2])
      paramSD <- sqrt(diag(solve(estResults$hessian)))
      paramsCI <- .getConfidenceIntervals(paramSD, estmu, estsig, level)

      structure(class = "dgxR", list(muHat = estmu,
                                      sigHat = estsig,
                                      muSD = paramSD[1],
                                      sigSD = paramSD[2],
                                      muCI = c(paramsCI[1,1],paramsCI[1,2]),
                                      sigCI = c(paramsCI[2,1],paramsCI[2,2]),
                                      logLikelihood = -estResults$value,
                                      hessian = estResults$hessian,
                                      call = Call))
    },
    error=function(cond) {
      print(cond)
      return(NA)
    })
}

#' @rdname dgxFit
#' @export
residuals.dgxR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  fitted.values <- fitted(object)
  residual.values <- as.numeric(dataMatrix[, 2]) - fitted.values
  return(residual.values)
}

#' @rdname dgxFit
#' @export
fitted.dgxR <- function(object, ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  N <- sum(as.numeric(dataMatrix[, 2]))
  fitted.values <- N*sapply(as.numeric(as.character(dataMatrix[,1])), pmfDGX, mu = object[['muHat']],
                            sig = object[['sigHat']])
  return(fitted.values)
}

#' @rdname dgxFit
#' @export
coef.dgxR <- function(object, ...){
  estimation <- matrix(nrow = 2, ncol = 4)
  estimation[1, ] <- c(object[['muHat']], object[['muSD']], object[['muCI']][1], object[['muCI']][2])
  estimation[2, ] <- c(object[['sigHat']], object[['sigSD']], object[['sigCI']][1], object[['sigCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("mu", "sig")
  estimation
}

#' @rdname dgxFit
#' @export
plot.dgxR <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  graphics::plot(dataMatrix[,1], dataMatrix[,2], log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting DGX Distribution", ...)

  graphics::lines(as.numeric(dataMatrix[,1]), fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'DGX Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname dgxFit
#' @export
print.dgxR <- function(x, ...){
  cat('Call:\n')
  print(x[['call']])
  cat('\n')
  cat('Initial Values:\n')
  cat(sprintf('mu: %s\n', format(eval(x[['call']]$init_mu), digits = 3)))
  cat(sprintf('sig: %s\n', format(eval(x[['call']]$init_sig), digits = 3)))
  cat('\n')
  cat('Coefficients:\n')
  print(coef(x))
  cat('\n')
  cat('Metrics:\n')
  cat(sprintf('Log-likelihood: %s\n', logLik(x)))
  cat(sprintf('AIC: %s\n', AIC(x)))
  cat(sprintf('BIC: %s\n', BIC(x)))
}

#' @rdname dgxFit
#' @export
summary.dgxR <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname dgxFit
#' @export
logLik.dgxR <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname dgxFit
#' @export
AIC.dgxR <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 2)
  return(aic)
}

#' @rdname dgxFit
#' @export
BIC.dgxR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  bic <- .get_BIC(object[['logLikelihood']], 2, sum(as.numeric(dataMatrix[, 2])))
  return(bic)
}












