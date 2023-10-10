source('estimation.R')

lambda0 <- seq(0.3, 0.7, by = 0.05)
c <- 1
dat <- readRDS('data.rds')
for (i in seq_along(dat)) {
  results <- lapply(dat[[i]], function(dat) {
    results <- lapply(lambda0, function(lambda0) {
      print(lambda0)
      print(system.time(result <- estimate(dat$Y, dat$D, dat$X, lambda0, c)))
      print(round(result$beta, 3))
      result
    })
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
  })
  saveRDS(results, paste0('simulation_', i, '.rds'))
}
