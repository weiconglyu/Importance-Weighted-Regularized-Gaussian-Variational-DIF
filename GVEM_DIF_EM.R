library(torch)

# rm(list = ls())
# load('em.RData')

varlist <- function(...) {
  setNames(list(...), as.character(substitute(alist(...)))[-1])
}

prox <- function(x, lambda) {
  sign(x) * (abs(x) - lambda)$maximum(0)
}

EM <- function(Y, D, X, eta, Sigma, Mu, a, b, gamma, beta, lambda, iter = 100, eps = 1e-3) {
  # save(list = ls(), file = 'em.RData')
  # stop()
  
  init <- function() {
    with(parent.frame(), {
      N <- nrow(Y)
      J <- ncol(Y)
      K <- max(X)
      Y <- torch_tensor(Y)
      eta <- torch_tensor(eta)
      Sigma <- torch_tensor(Sigma)$permute(c(3, 1, 2))$contiguous()
      SIGMA <- Sigma[X]
      Mu <- torch_tensor(t(Mu))
      MU <- Mu[X]
      a.mask <- torch_tensor(t(D) != 0)
      a <- torch_tensor(a)$mul(a.mask)
      a.mask.diag <- a.mask$diag_embed()$type_as(a)
      b <- torch_tensor(b)$view(-1)
      gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(dim(gamma)[1], a.mask)))
      gamma <- torch_cat(list(torch_zeros_like(gamma.mask[1:1]), torch_tensor(gamma)$transpose(2, 3)))$mul(gamma.mask)
      beta.mask <- torch_cat(list(torch_zeros_like(beta[, 1, drop = F]), torch_ones_like(beta)), 2)$t()$bool()
      beta <- torch_cat(list(torch_zeros_like(beta.mask[1:1]), torch_tensor(t(beta))))$mul(beta.mask)
    })
  }
  
  parameters <- function() {
    varlist(SIGMA, MU, Sigma, Mu, a, b, gamma, beta)
  }
  
  update <- function(lambda, gamma.mask, beta.mask) {
    with(parent.frame(), {
      args <- sys.frame(-4)
      aG <- (a + gamma)$unsqueeze(4)
      aG.t <- aG$transpose(3, 4)
      AG <- aG[X]
      AG.t <- AG$transpose(3, 4)
      BB <- (b - beta)[X]
      
      Sigma.inv <- Sigma$inverse()
      SIGMA.inv <- Sigma.inv[X] + 2 * (eta$view(c(N, -1, 1, 1)) * aG$matmul(aG.t)[X])$sum(2)
      SIGMA <- SIGMA.inv$inverse()
      MU <- SIGMA$matmul(((Y - 0.5 + 2 * eta * BB)$view(c(N, -1, 1, 1)) * AG)$sum(2) + Sigma.inv$matmul(Mu$unsqueeze(3))[X])$squeeze(3)
      Mu <- torch_stack(tapply(1:N, X, function(n) {
        MU[n]$mean(1)
      }))
      mu <- MU - Mu[X]
      sigma <- SIGMA + mu$unsqueeze(3)$matmul(mu$unsqueeze(2))
      Sigma <- torch_stack(tapply(1:N, X, function(n) {
        sigma[n]$mean(1)
      }))
      
      MU$unsqueeze_(2)
      sigma.mu <- a.mask.diag$transpose(2, 3)$matmul((SIGMA + MU$transpose(2, 3)$matmul(MU))$unsqueeze(2))$matmul(a.mask.diag)
      mu <- MU$mul(a.mask)
      xi <- (BB$square() - 2 * BB * AG.t$matmul(mu$unsqueeze(4))$view(c(N, -1)) + AG.t$matmul(sigma.mu)$matmul(AG)$view(c(N, -1)))$sqrt()#$mul((AG.t$matmul(mu$unsqueeze(4))$view(c(N, -1)) - BB)$sign())
      eta <- torch_where(abs(xi) < 1e-3, 0.125, (1 / (1 + exp(-xi)) - 0.5) / (2 * xi))
      a <- (2 * eta$view(c(N, -1, 1, 1)) * sigma.mu)$sum(1)$pinverse()$matmul(
        ((Y$unsqueeze(3) - 0.5) * mu + 2 * eta$unsqueeze(3) * (BB$unsqueeze(3) * mu - sigma.mu$matmul(gamma[X]$unsqueeze(4))$squeeze(4)))$sum(1)$unsqueeze(3))$squeeze(3)$masked_fill(!a.mask, 0)
      b <- ((0.5 - Y) + 2 * eta * (beta[X] + mu$unsqueeze(3)$matmul(AG)$view(c(N, -1))))$sum(1) / (2 * eta$sum(1))
      gamma <- torch_stack(tapply(1:N, X, function(n) {
        prox(((Y[n]$unsqueeze(3) - 0.5) * mu[n] + 2 * eta[n]$unsqueeze(3) *
                          (BB[n]$unsqueeze(3) * mu[n] - sigma.mu[n]$matmul(a$unsqueeze(3))$squeeze(4)))$sum(1), args$lambda) /
          (2 * eta[n]$view(c(length(n), -1, 1, 1)) * sigma.mu[n])$sum(1)$diagonal(dim1 = 2, dim2 = 3)
      }))$masked_fill(!args$gamma.mask, 0)
      beta <- torch_stack(tapply(1:N, X, function(n) {
        prox(((Y[n] - 0.5) + 2 * eta[n] * (b$unsqueeze(1) - mu[n]$unsqueeze(3)$matmul(AG[n])$view(c(length(n), -1))))$sum(1), args$lambda) / (2 * eta[n]$sum(1))
      }))$masked_fill(!args$beta.mask, 0)
      MU$squeeze_(2)
      
      mu <- Mu[1]$clone()
      MU$sub_(mu)
      Mu$sub_(mu)
      b$sub_(a$matmul(mu))
      beta$add_(gamma$matmul(mu))
      rescale <- Sigma[1]$diag()$sqrt()
      a$mul_(rescale)
      gamma$mul_(rescale)
      rescale.inv <- (1 / rescale)$diag()
      sigma <- rescale.inv$t()$matmul(SIGMA)$matmul(rescale.inv)
      SIGMA <- (sigma + sigma$transpose(2, 3)) / 2
      sigma <- rescale.inv$t()$matmul(Sigma)$matmul(rescale.inv)
      Sigma <- (sigma + sigma$transpose(2, 3)) / 2
    })
  }
  
  init()
  params.old <- parameters()
  for (i in 1:iter) {
    update(lambda, gamma.mask, beta.mask)
    update(0, gamma != 0, beta != 0)

    params.new <- parameters()
    dis <- mapply(params.old, params.new, FUN = function(x, y) {
      as.array(max(abs(x - y)))
    })

    if (all(dis < eps))
      break
    else
      params.old <- params.new
  }
  c(n = i, varlist(xi, eta, SIGMA, MU, Sigma, Mu, a, b, gamma, beta))
}

# library(cubature)
# library(mvnfast)
# IC <- function(Y, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, cons) {
#   N <- nrow(Y)
#   Q <- sum(sapply(1:N, function(n) {
#     ag <- a + gamma[X[n], , ]
#     bb <- b - beta[X[n]]
#     mu <- Mu[X[n], ]
#     sigma <- Sigma[X[n], , ]
#     density <- function(theta) {
#       t(colProds(dbinom(Y[n, ], 1, 1 / (1 + exp(-(ag %*% theta - bb))))) * dmvn(t(theta), mu, sigma))
#     }
#     log(hcubature(density, rep(-2, 2), rep(2, 2), vectorInterface = T)$integral)
#   }))
#   l0 <- as.array(sum(gamma != 0) + sum(beta != 0))
#   c(ll = Q, BIC = -2 * Q + l0 * log(N), GIC = -2 * Q + cons * l0 * log(N) * log(log(N)))
# }


IC <- function(Y, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, cons) {
  N <- nrow(Y)
  K <- max(X)
  Y <- torch_tensor(Y)
  AG <- (a + gamma)[X]$unsqueeze(4)
  AG.t <- AG$transpose(3, 4)
  BB <- (b - beta)[X]
  MU <- MU$clone()$unsqueeze_(2)
  SIGMA.MU <- (SIGMA + MU$transpose(2, 3)$matmul(MU))$unsqueeze(2)

  MU$unsqueeze_(4)
  xi <- (BB$square() - 2 * BB * AG.t$matmul(MU)$view(c(N, -1)) + AG.t$matmul(SIGMA.MU)$matmul(AG)$view(c(N, -1)))$sqrt()$mul((AG.t$matmul(MU)$view(c(N, -1)) - BB)$sign())
  mu <- MU$squeeze(2) - Mu[X]$unsqueeze(3)
  Q <- as_array((nnf_logsigmoid(xi) + (0.5 - Y) * (BB - AG.t$matmul(MU)$view(c(N, -1))) - xi / 2)$sum()
                - (Sigma$logdet()[X]$sum() + (linalg_solve(Sigma[X], SIGMA + mu$matmul(mu$transpose(2, 3)))$diagonal(dim1 = 2, dim2 = 3))$sum()) / 2)
  l0 <- as_array(sum(gamma != 0) + sum(beta != 0))
  c(ll = Q, AIC = -2 * Q + l0 * 2, BIC = -2 * Q + l0 * log(N), GIC = -2 * Q + cons * l0 * log(N) * log(log(N)))
}


#tmp <- EM(Y, D, X, eta, Sigma, Mu, a, b, gamma, beta, lambda)
