source('estimation.R')

Lambda0 <- seq(0.2, 0.7, by = 0.025)
c <- 0.75
S <- 10
M <- 10
dat <- readRDS('data.rds')
for (i in 8:length(dat))
  saveRDS(lapply(dat[[i]], function(dat) {
    results <- estimate(dat$Y, dat$D, dat$X, Lambda0, c, S, M)
    AIC <- sapply(results, `[[`, 'AIC')
    BIC <- sapply(results, `[[`, 'BIC')
    GIC <- sapply(results, `[[`, 'GIC')
    result.AIC <- results[[which.min(AIC)]]
    result.BIC <- results[[which.min(BIC)]]
    result.GIC <- results[[which.min(GIC)]]
    print(c(result.AIC$lambda0, result.BIC$lambda0, result.GIC$lambda0))
    print(round(result.AIC$beta, 3))
    print(round(result.BIC$beta, 3))
    print(round(result.GIC$beta, 3))
    list(results = results, AIC = AIC, BIC = BIC, GIC = GIC, result.AIC = result.AIC, result.BIC = result.BIC, result.GIC = result.GIC)
  }), paste0('simulation_', i, '.rds'))
