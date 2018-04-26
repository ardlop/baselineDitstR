# Metrics --------------------

.get_AIC <- function(loglike, K) {
  -2*loglike + 2*K
}

.get_BIC <- function(loglike, K, N) {
  -2*loglike + K*log(N)
}

# Utils ----------------

.getConfidenceIntervals <- function(paramSD, alpha, beta, level){
  result <- matrix(nrow=2, ncol=2)
  levelCoef <- round(stats::qnorm(1-((1-level)/2)), 2)
  offset <- levelCoef * paramSD
  result[1, ] <- c(alpha - offset[1], alpha + offset[1])
  result[2, ] <- c(beta - offset[2], beta + offset[2])
  colnames(result) <- c('Inf. CI', 'Sup. CI')
  rownames(result) <- c('alpha', 'beta')
  return(result)
}
