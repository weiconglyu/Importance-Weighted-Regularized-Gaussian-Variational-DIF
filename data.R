library(abind)
library(mvnfast)
set.seed(202310)

R <- 20

Sigma <- matrix(c(1, 0.85, 0.85, 1), 2)
J <- 20
D <- cbind(rep(1:0, J / 2), rep(0:1, J / 2))
j <- J * c(0.2, 0.6)

dat.unif <- unlist(lapply(c(500, 1000), function(n) {
  X <- rep(1:3, each = n)
  lapply(j, function(j) {
    replicate(R, {
      a <- matrix(runif(J * 2, 1.5, 2.5), ncol = 2) * D
      a <- abind(a, a, a, along = 0)
      b <- rnorm(J)
      b <- unname(rbind(b, b, b))
      b[-1, 1:(j / 2)] <- b[-1, 1:(j / 2)] - c(0.5, 1)
      b[-1, (j / 2 + 1):j] <- b[-1, (j / 2 + 1):j] + c(0.5, 1)
      theta <- rmvn(n * 3, rep(0, 2), Sigma)
      Y <- t(sapply(1:(n * 3), function(n) {
        rbinom(J, 1, 1 / (1 + exp(-(a[X[n], , ] %*% theta[n, ] - b[X[n], ]))))
      }))
      list(Y = Y, D = D, X = X, params = list(a = a, b = b, theta = theta))
    }, simplify = F)
  })
}), F)

dat.nonunif <- unlist(lapply(c(500, 1000), function(n) {
  X <- rep(1:3, each = n)
  lapply(j, function(j) {
    replicate(M, {
      a <- matrix(runif(J * 2, 1.5, 2.5), ncol = 2) * D
      a <- abind(a, a, a, along = 0)
      a[-1, 1:(j / 2), ] <- a[-1, 1:(j / 2), ] - c(0.4, 0.6)
      a[-1, (j / 2 + 1):j, ] <- a[-1, (j / 2 + 1):j, ] + c(0.4, 0.6)
      a[-1, , ] <- a[-1, , ] * abind(D, D, along = 0)
      b <- rnorm(J)
      b <- unname(rbind(b, b, b))
      b[-1, 1:(j / 2)] <- b[-1, 1:(j / 2)] - c(0.25, 0.6)
      b[-1, (j / 2 + 1):j] <- b[-1, (j / 2 + 1):j] + c(0.25, 0.6)
      theta <- rmvn(n * 3, rep(0, 2), Sigma)
      Y <- t(sapply(1:(n * 3), function(n) {
        rbinom(J, 1, 1 / (1 + exp(-(a[X[n], , ] %*% theta[n, ] - b[X[n], ]))))
      }))
      list(Y = Y, D = D, X = X, params = list(a = a, b = b, theta = theta))
    }, simplify = F)
  })
}), F)

saveRDS(c(dat.unif, dat.nonunif), file = 'data.rds')
