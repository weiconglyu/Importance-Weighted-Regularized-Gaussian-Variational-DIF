library(torch)
library(matrixStats)

raw <- lapply(1:8, function(r) {
  N <- c(1500, 3000)[((r - 1) %% 4 > 1) + 1]
  J <- 20 * c(0.2, 0.6)[(r - 1) %% 2 + 1]
  rep <- readRDS(paste0('simulation_', r, '.rds'))
  list(N = N, J = J, rep = rep)
})

c.GIC <- seq(0, 1, by = 0.1)
best <- lapply(raw, function(raw) {
  N <- raw$N
  J <- raw$J
  sapply(raw$rep, function(rep) {
    ic <- sapply(rep$results, function(rep) {
      c(dev = -2 * rep$ll, l0 = sum(rep$gamma != 0) + sum(rep$beta != 0))
    })
    setNames(sapply(c.GIC, function(c) {
      which.min(ic[1, ] + c * ic[2, ] * log(N) * log(log(N)))
    }), c.GIC)
  })
})

c.GIC <- 0.7
result <- lapply(raw, function(raw) {
  N <- raw$N
  J <- raw$J
  rep <- lapply(raw$rep, function(rep) {
    ic <- sapply(rep$results, function(rep) {
      c(dev = -2 * rep$ll, l0 = sum(rep$gamma != 0) + sum(rep$beta != 0))
    })
    best <- c(BIC = which.min(ic[1, ] + ic[2, ] * log(N)), GIC = which.min(ic[1, ] + c.GIC * ic[2, ] * log(N) * log(log(N))))
    lapply(rep$results[best], function(result) {
      beta <- result$beta[-1, ]
      gamma <- result$gamma[-1, , ]
      lambda0 <- result$lambda0
      dif <- (beta != 0) | (rowSums(gamma != 0, dims = 2) != 0)
      d <- dif[, 1:J]
      power <- c(mean(colAnys(d)), mean(d[1, ]), mean(d[2, ]))
      d <- dif[, -(1:J)]
      typeI <- c(mean(colAnys(d)), mean(d[1, ]), mean(d[2, ]))
      prop <- rbind(power, typeI)
      colnames(prop) <- c('Omnibus', 'Low', 'High')
      list(prop = prop, lambda0 = lambda0)
    })
  })
  prop <- torch_stack(lapply(rep, function(rep) {
    torch_stack(lapply(rep, `[[`, 'prop'))
  }))
  prop <- lapply(torch_cat(torch_std_mean(prop, 1))$unbind(), function(prop) {
    prop <- as.array(prop)
    rownames(prop) <- c('Power', 'TypeI')
    colnames(prop) <- c('Omnibus', 'Low', 'High')
    prop
  })
  lambda0 <- do.call(rbind, lapply(rep, function(rep) {
    unlist(lapply(rep, `[[`, 'lambda0'))
  }))
  lambda0 <- as.data.frame(rbind(colMeans(lambda0), colMins(lambda0), colMaxs(lambda0)))
  result <- mapply(rbind, prop, lambda0 = cbind(lambda0, lambda0), SIMPLIFY = F)
  mean <- setNames(result[3:4], c('BIC', 'GIC'))
  sd <- setNames(result[1:2], c('BIC', 'GIC'))
  list(mean = mean, sd = sd)
})
lapply(result, `[[`, 'mean')
