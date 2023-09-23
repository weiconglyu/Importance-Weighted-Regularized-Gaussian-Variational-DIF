library(matrixStats)

J.dif <- list(c(4, 5, 12, 13), c(4:9, 12:17))
ics <- c('aic', 'bic', 'gic')
summarize <- function(sim, cond) {
  J <- J.dif[[(cond - 1) %% 2  + 1]]
  dif <- torch_stack(lapply(1:20, function(rep) {
    load(paste0('GVEM_Sim', sim, '_Condition', cond, '_', rep, '.RData'))
    torch_stack(lapply(ics, function(ic) {
      result <- eval(parse(text = paste0('result.', ic)))
      #dif <- apply(result$gamma[-1, , ] != 0, 1:2, any) | (result$beta[-1, ] != 0)
      dif <- result$beta[-1, ] != 0
      d <- dif[, J]
      power <- c(mean(colAnys(d)), mean(d[1, ]), mean(d[2, ]))
      d <- dif[, -J]
      typeI <- c(mean(colAnys(d)), mean(d[1, ]), mean(d[2, ]))
      rbind(power, typeI)
    }))
  }))
  setNames(lapply(dif$mean(1)$unbind(), function(mean) {
    result <- as_array(mean)
    rownames(result) <- c('Power', 'TypeI')
    colnames(result) <- c('Total', 'LowDIF', 'HighDIF')
    result
  }), ics)
}

summarize(3, 8)
