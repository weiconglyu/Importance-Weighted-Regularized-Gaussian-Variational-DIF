library(torch)
library(matrixStats)

ics <- c('AIC', 'BIC', 'GIC')

results <- lapply(1:8, function(r) {
  j <- 20 * c(0.2, 0.6)[(r - 1) %% 2 + 1]
  reps <- lapply(readRDS(paste0('simulation_', r, '.rds')), function(rep) {
    lapply(ics, function(ic) {
      result <- rep[[paste0('result.', ic)]]
      beta <- result$beta[-1, ]
      gamma <- result$gamma[-1, , ]
      lambda0 <- result$lambda0
      dif <- (beta != 0) | (rowSums(gamma != 0, dims = 2) != 0)
      d <- dif[, 1:j]
      power <- c(mean(colAnys(d)), mean(d[1, ]), mean(d[2, ]))
      d <- dif[, -(1:j)]
      typeI <- c(mean(colAnys(d)), mean(d[1, ]), mean(d[2, ]))
      prop <- rbind(power, typeI)
      colnames(prop) <- c('Omnibus', 'Low', 'High')
      list(prop = prop, lambda0 = lambda0)
    })
  })
  prop <- torch_stack(lapply(reps, function(rep) {
    torch_stack(lapply(rep, `[[`, 'prop'))
  }))
  props <- lapply(torch_cat(torch_std_mean(prop, 1))$unbind(), function(prop) {
    prop <- as.array(prop)
    rownames(prop) <- c('Power', 'TypeI')
    colnames(prop) <- c('Omnibus', 'Low', 'High')
    prop
  })
  lambda0 <- do.call(rbind, lapply(reps, function(rep) {
    unlist(lapply(rep, `[[`, 'lambda0'))
  }))
  lambda0s <- as.data.frame(rbind(colMeans(lambda0), colMins(lambda0), colMaxs(lambda0)))
  result <- mapply(rbind, props, lambda0 = cbind(lambda0s, lambda0s), SIMPLIFY = F)
  mean <- setNames(result[4:6], ics)
  sd <- setNames(result[1:3], ics)
  list(mean = mean, sd = sd)
})

results[[1]]$mean
results[[2]]$mean
results[[3]]$mean
results[[4]]$mean
results[[5]]$mean
results[[6]]$mean
results[[7]]$mean
results[[8]]$mean
