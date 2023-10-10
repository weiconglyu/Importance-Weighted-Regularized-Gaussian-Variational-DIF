library(torch)
library(matrixStats)

ics <- c('AIC', 'BIC', 'GIC')

results <- lapply(1:8, function(r) {
  j <- 20 * c(0.2, 0.6)[(r - 1) %% 2 + 1]
  reps <- torch_stack(lapply(readRDS(paste0('simulation_', r, '.rds')), function(rep) {
    torch_stack(lapply(ics, function(ic) {
      beta <- rep[[paste0('result.', ic)]]$beta[-1, ]
      gamma <- rep[[paste0('result.', ic)]]$gamma[-1, , ]
      dif <- (beta != 0) | (rowSums(gamma != 0, dims = 2) != 0)
      d <- dif[, 1:j]
      power <- c(mean(colAnys(d)), mean(d[1, ]), mean(d[2, ]))
      d <- dif[, -(1:j)]
      typeI <- c(mean(colAnys(d)), mean(d[1, ]), mean(d[2, ]))
      prop <- rbind(power, typeI)
      colnames(prop) <- c('Omnibus', 'Low', 'High')
      prop
    }))
  }))
  props <- lapply(torch_cat(torch_std_mean(reps, 1))$unbind(), function(prop) {
    prop <- as.array(prop)
    rownames(prop) <- c('Power', 'TypeI')
    colnames(prop) <- c('Omnibus', 'Low', 'High')
    prop
  })
  mean <- setNames(props[4:6], ics)
  sd <- setNames(props[1:3], ics)
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
