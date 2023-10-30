library(torch)

`%*%.torch_tensor` <- function(e1, e2) {
  torch_matmul(e1, e2)
}

t.torch_tensor <- function(e) {
  torch_transpose(e, -1, -2)
}

diagonal <- function(e) {
  torch_diagonal(e, dim1 = -1, dim2 = -2)
}

sym <- function(e) {
  (e + t(e)) / 2
}

distance <- function(x, y) {
  mapply(function(x, y) {
    max(abs(x - y))
  }, x, y)
}

EM <- function(Y, D, X, iter = 1000, eps = 1e-3) {
  init <- function() {
    with(parent.frame(), {
      N <- nrow(Y)
      J <- ncol(Y)
      K <- ncol(D)
      G <- max(X)
      
      Y <- torch_tensor(Y)
      eta <- torch_tensor(matrix(0.125, N, J))
      Sigma <- torch_stack(replicate(G, diag(K), simplify = F))
      Mu <- torch_zeros(c(G, K))
      a.mask <- torch_tensor(D * 1.0)
      a.mask.diag <- a.mask$diag_embed()
      a <- a.mask$clone()
      b <- torch_zeros(J)
      gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(G - 1, a.mask)))
      gamma <- torch_zeros_like(gamma.mask)
      beta.mask <- torch_cat(list(torch_zeros(c(1, J)), torch_ones(c(G - 1, J))))
      beta <- torch_zeros_like(beta.mask)
    })
  }
  
  parameters <- function() {
    lapply(list(SIGMA = SIGMA, MU = MU, Sigma = Sigma, Mu = Mu, a = a, b = b, gamma = gamma, beta = beta), as.array)
  }

  update <- function() {
    with(parent.frame(), {
      aG <- (a + gamma)$unsqueeze(4)
      aG.t <- t(aG)
      AG <- aG[X]
      AG.t <- t(AG)
      BB <- (b - beta)[X]
      
      Sigma.inv <- Sigma$inverse()
      SIGMA.inv <- Sigma.inv[X] + 2 * (eta$view(c(N, -1, 1, 1)) * (aG %*% aG.t)[X])$sum(2)
      SIGMA <- sym(SIGMA.inv$inverse())
      MU <- (SIGMA %*% (((Y - 0.5 + 2 * eta * BB)$view(c(N, -1, 1, 1)) * AG)$sum(2) + (Sigma.inv %*% Mu$unsqueeze(3))[X]))$squeeze(3)
      Mu <- torch_stack(tapply(1:N, X, function(n) {
        MU[n]$mean(1)
      }))
      mu <- MU - Mu[X]
      sigma.mu <- SIGMA + mu$unsqueeze(3) %*% mu$unsqueeze(2)
      Sigma <- sym(torch_stack(tapply(1:N, X, function(n) {
        sigma.mu[n]$mean(1)
      })))
      
      mu <- a.mask.diag %*% MU$view(c(N, 1, -1, 1))
      sigma.mu <- a.mask.diag %*% SIGMA$unsqueeze(2) %*% a.mask.diag + mu %*% t(mu)
      xi <- sqrt(BB$square() - 2 * BB * (AG.t %*% mu)$view(c(N, -1)) + (AG.t %*% sigma.mu %*% AG)$view(c(N, -1)))
      eta <- torch_where(abs(xi) < 1e-3, 0.125, (1 / (1 + exp(-xi)) - 0.5) / (2 * xi))
      a <- ((2 * eta$view(c(N, -1, 1, 1)) * sigma.mu)$sum(1)$pinverse() %*% ((Y - 0.5)$view(c(N, -1, 1, 1)) * mu + 2 * eta$view(c(N, -1, 1, 1)) * (BB$view(c(N, -1, 1, 1)) * mu - sigma.mu %*% gamma[X]$unsqueeze(4)))$sum(1))$squeeze(3) * a.mask
      b <- (0.5 - Y + 2 * eta * (beta[X] + (AG.t %*% mu)$view(c(N, -1))))$sum(1) / (2 * eta$sum(1))
      gamma.beta <- torch_stack(tapply(1:N, X, function(n) {
        N <- length(n)
        torch_cat(list(((2 * eta[n]$view(c(N, -1, 1, 1)) * sigma.mu[n])$sum(1)$pinverse() %*% ((Y[n] - 0.5)$unsqueeze(3) * mu[n]$squeeze(4) + 2 * eta[n]$unsqueeze(3) * (BB[n]$unsqueeze(3) * mu[n]$squeeze(4) - (sigma.mu[n] %*% a$unsqueeze(3))$squeeze(4)))$sum(1)$unsqueeze(3))$squeeze(3),
                       (((Y[n] - 0.5) + 2 * eta[n] * (b - (AG.t[n] %*% mu[n])$view(c(N, -1))))$sum(1) / (2 * eta[n]$sum(1)))$unsqueeze(2)), 2)
      }))
      gamma$set_data(gamma.beta[, , 1:K] * gamma.mask)
      beta$set_data(gamma.beta[, , K + 1] * beta.mask)

      mu <- Mu[1]$clone()
      MU$sub_(mu)
      Mu$sub_(mu)
      b$sub_(a %*% mu)
      beta$add_(gamma %*% mu)
      sigma <- Sigma[1]$diag()$sqrt()
      a$mul_(sigma)
      gamma$mul_(sigma)
      sigma.inv <- (1 / sigma)$diag()
      SIGMA$set_data(sym(sigma.inv %*% SIGMA %*% sigma.inv))
      Sigma$set_data(sym(sigma.inv %*% Sigma %*% sigma.inv))
    })
  }

  init()
  params.old <- NULL
  for (i in 1:iter) {
    update()
    params <- parameters()
    if (!is.null(params.old) && all(distance(params, params.old) < eps))
      break
    params.old <- params
  }
  print(paste('EM:', i))
  params
}

IW <- function(Y, D, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, lambda, z, theta.logd, S, M, iter = 1000, eps = 1e-3, lr = 0.1) {
  init <- function() {
    with(parent.frame(), {
      N <- nrow(Y)
      J <- ncol(Y)
      K <- ncol(D)
      G <- max(X)
      
      Y <- torch_tensor(Y)$unsqueeze(2)$bool()
      SIGMA.L <- linalg_cholesky(SIGMA)$unsqueeze(2)$contiguous()
      MU <- torch_tensor(MU)$unsqueeze(2)
      Sigma.L <- linalg_cholesky(Sigma)
      Sigma1.L.mask <- torch_tensor(lower.tri(matrix(0, K, K)))$bool()
      Sigma1.L.v <- (Sigma.L[1] / torch_cat(list(torch_ones(c(K, 1)), sqrt(1 - Sigma.L[1]$square()$cumsum(2))[, 1:(K - 1)]), 2))$masked_select(Sigma1.L.mask)$arctanh()$requires_grad_(T)
      Sigma2.L.mask <- torch_stack(replicate(G - 1, lower.tri(matrix(0, K, K), T), F))$bool()
      Sigma2.L.v <- Sigma.L[2:G]$masked_select(Sigma2.L.mask)$requires_grad_(T)
      Mu.mask <- torch_cat(list(matrix(0, 1, K), matrix(1, G - 1, K)))$bool()
      Mu.v <- torch_tensor(Mu)$masked_select(Mu.mask)$requires_grad_(T)
      a.mask <- torch_tensor(D)$bool()
      a.v <- torch_tensor(a)$masked_select(a.mask)$requires_grad_(T)
      b <- torch_tensor(b)$requires_grad_(T)
      gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(G - 1, a.mask)))
      gamma.v <- torch_tensor(gamma)$masked_select(gamma.mask)$requires_grad_(T)
      beta.mask <- torch_cat(list(torch_zeros(c(1, J)), torch_ones(c(G - 1, J))))$bool()
      beta.v <- torch_tensor(beta)$masked_select(beta.mask)$requires_grad_(T)
      theta.logd <- torch_tensor(theta.logd)
      theta <- (SIGMA.L %*% torch_tensor(z)$unsqueeze(4))$squeeze(4) + MU
    })
  }
  
  assemble <- function() {
    with(parent.frame(), {
      z <- Sigma1.L.v$tanh()
      Sigma1.L <- torch_eye(K)$masked_scatter(Sigma1.L.mask, z) * torch_cat(c(torch_ones(c(K, 1)), (1 - torch_zeros(c(K, K))$masked_scatter(Sigma1.L.mask, z)[, 1:(K- 1)]$square())$cumprod(2)$sqrt()), 2)
      Sigma2.L <- torch_zeros(Sigma2.L.mask$shape)$masked_scatter(Sigma2.L.mask, Sigma2.L.v)
      Sigma.L <- torch_cat(c(Sigma1.L$unsqueeze(1), Sigma2.L))
      Mu <- torch_zeros(Mu.mask$shape)$masked_scatter(Mu.mask, Mu.v)
      a <- torch_zeros(a.mask$shape)$masked_scatter(a.mask, a.v)
      gamma <- torch_zeros(gamma.mask$shape)$masked_scatter(gamma.mask, gamma.v)
      beta <- torch_zeros(beta.mask$shape)$masked_scatter(beta.mask, beta.v)
    })
  }
  
  proximal <- function(x, state, params, lambda) {
    beta.t <- params$betas ^ state$step
    lr <- params$lr / (1 - beta.t[1]) / (sqrt(state$exp_avg_sq / (1 - beta.t[2])) + params$eps)
    x$set_data(sign(x) * (abs(x) - lr * lambda)$maximum(0))
  }
  
  parameters <- function() {
    with_no_grad({
      assemble()
      lapply(list(Sigma = Sigma.L %*% t(Sigma.L), Sigma.L = Sigma.L, Mu = Mu, a = a, b = b, gamma = gamma, beta = beta), as.array)
    })
  }
  
  update <- function(lambda, gamma.v.mask, beta.v.mask) {
    force(gamma.v.mask)
    force(beta.v.mask)
    
    with(parent.frame(), {
      args <- sys.frame(-4)
      
      opt <- optim_adam(list(gamma.v, beta.v, Sigma1.L.v, Sigma2.L.v, Mu.v, a.v, b), lr)
      params.old <- NULL
      for (i in 1:iter) {
        assemble()
        xi <- ((a + gamma)[X]$unsqueeze(2) * theta$unsqueeze(3))$sum(4) - (b - beta)[X]$unsqueeze(2)
        # with_no_grad(Sigma.L$add_(diagonal(Sigma.L)$sign()$diag_embed() * 1e-8))
        log.w <- torch_where(Y, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$sum(3) +
          torch_cat(lapply(1:G, function(g) {
            distr_multivariate_normal(Mu[g], scale_tril = Sigma.L[g])$log_prob(theta[X == g])
          })) - theta.logd
        Q <- (log.w$view(c(N, S, M))$logsumexp(3) - log(M))$mean(2)$sum()
        opt$zero_grad()
        (-Q)$backward()
        with_no_grad({
          gamma.v$grad$mul_(args$gamma.v.mask)
          beta.v$grad$mul_(args$beta.v.mask)
          opt$step()
          if (args$lambda != 0) {
            opt.state <- opt$state_dict()$state
            proximal(gamma.v, opt.state[[1]], opt$param_groups[[1]], args$lambda)
            proximal(beta.v, opt.state[[2]], opt$param_groups[[1]], args$lambda)
          }
        })

        params <- parameters()
        if (!is.null(params.old) && all(distance(params, params.old)[-1] < eps))
          break
        params.old <- params
      }
      print(paste('IW:', i))
      params
    })
  }
  
  init()
  # z <- array(rnorm(N * S * M * K), c(N, S * M, K))
  # theta.logd <- torch_tensor(rowSums(dnorm(z, log = T), dim = 2))
  # theta <- (SIGMA.L %*% z)$squeeze(4) + MU
  update(lambda, torch_ones_like(gamma.v)$bool(), torch_ones_like(beta.v)$bool())
  update(0, gamma.v != 0, beta.v != 0)
}

IC <- function(Y, X, SIGMA, MU, Sigma.L, Mu, a, b, gamma, beta, c, z, theta.logd, S, M) {
  N <- nrow(Y)
  G <- max(X)
  
  Y <- torch_tensor(Y)$unsqueeze(2)$bool()
  SIGMA.L <- linalg_cholesky(SIGMA)$unsqueeze(2)
  MU <- torch_tensor(MU)$unsqueeze(2)
  Sigma.L <- torch_tensor(Sigma.L)
  Mu <- torch_tensor(Mu)
  a <- torch_tensor(a)
  b <- torch_tensor(b)
  gamma <- torch_tensor(gamma)
  beta <- torch_tensor(beta)
  theta.logd <- torch_tensor(theta.logd)
  theta <- (SIGMA.L %*% torch_tensor(z)$unsqueeze(4))$squeeze(4) + MU
  
  xi <- ((a + gamma)[X]$unsqueeze(2) * theta$unsqueeze(3))$sum(4) - (b - beta)[X]$unsqueeze(2)
  log.w <- torch_where(Y, nnf_logsigmoid(xi), nnf_logsigmoid(-xi))$sum(3) +
    torch_cat(lapply(1:G, function(g) {
      distr_multivariate_normal(Mu[g], scale_tril = Sigma.L[g])$log_prob(theta[X == g])
    })) - theta.logd
  Q <- as.array((log.w$view(c(N, S, M))$logsumexp(3) - log(M))$mean(2)$sum())
  l0 <- as.array(sum(gamma != 0) + sum(beta != 0))
  c(ll = Q, AIC = -2 * Q + l0 * 2, BIC = -2 * Q + l0 * log(N), GIC = -2 * Q + c * l0 * log(N) * log(log(N)))
}

estimate <- function(Y, D, X, Lambda0, c, S, M, iter = 1000, eps = 1e-3, lr = 0.1) {
  N <- nrow(Y)
  K <- ncol(D)
  z <- array(rnorm(N * S * M * K), c(N, S * M, K))
  theta.logd <- rowSums(dnorm(z, log = T), dim = 2)

  print(system.time(em <- EM(Y, D, X, iter, eps)))
  lapply(Lambda0, function(lambda0) {
    lambda <- sqrt(N) * lambda0
    print(c(lambda0, lambda))
    list2env(em, environment())
    print(system.time(iw <- IW(Y, D, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, lambda, z, theta.logd, S, M, iter, eps, lr)))
    list2env(iw, environment())
    list2env(as.list(IC(Y, X, SIGMA, MU, Sigma.L, Mu, a, b, gamma, beta, c, z, theta.logd, S, M)), environment())
    list(lambda0 = lambda0, lambda = lambda, SIGMA = SIGMA, MU = MU, Sigma = Sigma, Mu = Mu, a = a, b = b, gamma = gamma, beta = beta, ll = ll, AIC = AIC, BIC = BIC, GIC = GIC)
  })
}
