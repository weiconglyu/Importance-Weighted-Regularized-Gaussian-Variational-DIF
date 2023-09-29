library(torch)

# rm(list = ls())
# load('em.RData')

varlist <- function(...) {
  setNames(list(...), as.character(substitute(alist(...)))[-1])
}

`%*%.torch_tensor` <- function(e1, e2) {
  torch_matmul(e1, e2)
}

`%<-%` <- function(e1, e2) {
  e1$set_data(e2)
  invisible(e1)
}

t.torch_tensor <- function(e) {
  torch_transpose(e, -1, -2)
}

sym <- function(e) {
  (e + t(e)) / 2
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
      Sigma <- sym(torch_tensor(aperm(Sigma, c(3, 1, 2))))
      SIGMA <- Sigma[X]
      Mu <- torch_tensor(t(Mu))
      MU <- Mu[X]
      a.mask <- torch_tensor(t(D) != 0)
      a <- torch_tensor(a) * a.mask
      a.mask.diag <- a.mask$diag_embed()$type_as(a)
      j.mask <- a.mask$diag_embed()$type_as(a)
      b <- torch_tensor(b)$view(-1)
      gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(dim(gamma)[1], a.mask)))
      gamma <- torch_cat(list(torch_zeros_like(gamma.mask[1:1]), t(torch_tensor(gamma)))) * gamma.mask
      beta.mask <- torch_cat(list(torch_zeros_like(beta[, 1, drop = F]), torch_ones_like(beta)), 2)$t()$bool()
      beta <- torch_cat(list(torch_zeros_like(beta.mask[1:1]), torch_tensor(t(beta)))) * beta.mask
      
      # Sigma <- torch_stack(replicate(Sigma$shape[1], torch_eye(Sigma$shape[2])))
      # Mu$mul_(0)
      # SIGMA <- torch_stack(replicate(SIGMA$shape[1], torch_eye(SIGMA$shape[2])))
      # MU$mul_(0)
      # a <- torch_ones_like(a) * a.mask
      # b$mul_(0)
      # gamma$mul_(0)
      # beta$mul_(0)
    })
  }
  
  parameters <- function() {
    varlist(SIGMA, MU, Sigma, Mu, a, b, gamma, beta)
  }
  
  update <- function(lambda, gamma.mask, beta.mask) {
    with(parent.frame(), {
      args <- sys.frame(-4)
      
      aG <- (a + gamma)$unsqueeze(4)
      aG.t <- t(aG)
      AG <- aG[X]
      AG.t <- t(AG)
      BB <- (b - beta)[X]
      sigma.mu <- (SIGMA + MU$unsqueeze(3) %*% MU$unsqueeze(2))$unsqueeze(2)
      xi <- sqrt(BB$square() - 2 * BB * (AG.t %*% MU$view(c(N, 1, -1, 1)))$view(c(N, -1)) + (AG.t %*% sigma.mu %*% AG)$view(c(N, -1)))
      eta <- torch_where(abs(xi) < 1e-3, 0.125, (1 / (1 + exp(-xi)) - 0.5) / (2 * xi))
      
      Sigma.inv <- Sigma$inverse()
      SIGMA.inv <- Sigma.inv[X] + 2 * (eta$view(c(N, -1, 1, 1)) * (aG %*% aG.t)[X])$sum(2)
      SIGMA <- sym(SIGMA.inv$inverse())
      MU <- (SIGMA %*% (((Y - 0.5 + 2 * eta * BB)$view(c(N, -1, 1, 1)) * AG)$sum(2) + (Sigma.inv %*% Mu$unsqueeze(3))[X]))$squeeze(3)
      Mu <- torch_stack(tapply(1:N, X, function(n) {
        MU[n]$mean(1)
      }))
      mu <- MU - Mu[X]
      sigma.mu <- SIGMA + mu$unsqueeze(3) %*% mu$unsqueeze(2)
      Sigma <- torch_stack(tapply(1:N, X, function(n) {
        sigma.mu[n]$mean(1)
      }))
      
      mu <- j.mask %*% MU$view(c(N, 1, -1, 1))
      sigma.mu <- j.mask %*% SIGMA$unsqueeze(2) %*% j.mask + mu %*% t(mu)
      a <- ((2 * eta$view(c(N, -1, 1, 1)) * sigma.mu)$sum(1)$pinverse() %*% ((Y - 0.5)$view(c(N, -1, 1, 1)) * mu + 2 * eta$view(c(N, -1, 1, 1)) * (BB$view(c(N, -1, 1, 1)) * mu - sigma.mu %*% gamma[X]$unsqueeze(4)))$sum(1))$squeeze(3) * a.mask
      b <- (0.5 - Y + 2 * eta * (beta[X] + (AG.t %*% mu)$view(c(N, -1))))$sum(1) / (2 * eta$sum(1))
      beta.sigma  <- torch_stack(tapply(1:N, X, function(n) {
        N <- length(n)
        torch_cat(list((prox(((Y[n] - 0.5) + 2 * eta[n] * (b - (AG.t[n] %*% mu[n])$view(c(N, -1))))$sum(1), args$lambda) / (2 * eta[n]$sum(1)))$unsqueeze(2),
                       prox(((Y[n] - 0.5)$unsqueeze(3) * mu[n]$squeeze(4) + 2 * eta[n]$unsqueeze(3) * (BB[n]$unsqueeze(3) * mu[n]$squeeze(4) - (sigma.mu[n] %*% a$unsqueeze(3))$squeeze(4)))$sum(1), args$lambda) / (2 * eta[n]$view(c(N, -1, 1, 1)) * sigma.mu[n])$sum(1)$diagonal(dim1 = 2, dim2 = 3)), 2)
      }))
      gamma <- beta.sigma[, , 2:beta.sigma$shape[3]]$masked_fill(!args$gamma.mask, 0)
      beta <- beta.sigma[, , 1] * args$beta.mask
      
      mu <- Mu[1]$clone()
      MU$sub_(mu)
      Mu$sub_(mu)
      b$sub_(a %*% mu)
      beta$add_(gamma %*% mu)
      rescale <- Sigma[1]$diag()$sqrt()
      a$mul_(rescale)
      gamma$mul_(rescale)
      rescale.inv <- (1 / rescale)$diag()
      SIGMA <- sym(rescale.inv %*% SIGMA %*% rescale.inv)
      Sigma <- sym(rescale.inv %*% Sigma %*% rescale.inv)
    })
  }
  
  init()
  params.old <- parameters()
  for (i in 1:iter) {
    update(lambda, gamma.mask, beta.mask)
    #update(0, gamma != 0, beta != 0)

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
  AG.t <- t(AG)
  BB <- (b - beta)[X]
  MU <- MU$clone()$unsqueeze(2)
  xi <- sqrt(BB$square() - 2 * BB * (AG.t %*% MU$view(c(N, 1, -1, 1)))$view(c(N, -1)) + (AG.t %*% (SIGMA + t(MU) %*% MU)$unsqueeze(2) %*% AG)$view(c(N, -1)))
  MU$unsqueeze_(4)
  mu <- MU$squeeze(2) - Mu[X]$unsqueeze(3)
  Q <- as.array((nnf_logsigmoid(xi) + (0.5 - Y) * (BB - (AG.t %*% MU)$view(c(N, -1))) - xi / 2)$sum() - (Sigma$logdet()[X]$sum() + (linalg_solve(Sigma[X], SIGMA + mu %*% t(mu))$diagonal(dim1 = 2, dim2 = 3))$sum()) / 2)
  l0 <- as.array(sum(gamma != 0) + sum(beta != 0))
  c(ll = Q, AIC = -2 * Q + l0 * 2, BIC = -2 * Q + l0 * log(N), GIC = -2 * Q + cons * l0 * log(N) * log(log(N)))
}


# tmp <- EM(Y, D, X, eta, Sigma, Mu, a, b, gamma, beta, lambda)
# list2env(tmp, environment())
# as_array(beta) |> round(3)
# as_array(a) |> t() |> round(3)
