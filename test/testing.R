library(baselineDistR)
dataRaw <- read.csv('/home/ariel/Documents/ac/trunk/Ziva_test/atCKN_deg-all.txt', header= T, sep='\t')
data <- as.data.frame(table(dataRaw[,3]))
# data[,1] <- as.numeric(levels(data[,1]))[data[,1]]
data[,1] <- as.numeric(data[,1])

# -------------------- Testing NB Zero Truncate Distribution -------------
nb_zt <- negbinZTFit(data, 0.5, 0.5)
coef(nb_zt)
logLik(nb_zt)
AIC(nb_zt)
BIC(nb_zt)
print(nb_zt)
fitted(nb_zt)
plot(nb_zt)
residuals(nb_zt)

# --------------- Testing Discrete Weibull Zero truncated distribution ---------------
dw_zt <- discWeibullZTFit(data, init_p = 0.5, init_v = 0.5)
coef(dw_zt)
logLik(dw_zt)
AIC(dw_zt)
BIC(dw_zt)
print(dw_zt)
fitted(dw_zt)
plot(dw_zt)
residuals(dw_zt)

# ---------------- Testing DGX distribution --------------------
dgx_zt <- dgxFit(data, init_mu = 1, init_sig = 5)
coef(dgx_zt)
logLik(dgx_zt)
AIC(dgx_zt)
BIC(dgx_zt)
print(dgx_zt)
fitted(dgx_zt)
plot(dgx_zt)
residuals(dgx_zt)

# ------------------- Testing Zipf Distribution  ----------------------
zipf_fit <- zipfFit(data, init_alpha = 1.5)
coef(zipf_fit)
logLik(zipf_fit)
AIC(zipf_fit)
BIC(zipf_fit)
print(zipf_fit)
fitted(zipf_fit)
plot(zipf_fit)
residuals(zipf_fit)


