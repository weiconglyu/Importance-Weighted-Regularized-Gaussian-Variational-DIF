library(torch)
library(mvnfast)

# rm(list = ls())
load('iw.RData')
# D <- A
# X <- G
# eps <- 1e-8

IW <- function(Y, D, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, lambda, S = 10, M = 10, lr = 0.05, iter = 1000, R = 1) {
  init <- function() {
    with(parent.frame(), {
      N <- nrow(Y)
      J <- ncol(Y)
      K <- max(X)
      o <- order(X)
      Y <- torch_tensor(Y[o, ])$unsqueeze(2)$bool()
      X <- X[o]
      SIGMA.L <- linalg_cholesky(torch_tensor(SIGMA[, , o])$permute(c(3, 1, 2)))$contiguous()
      MU <- torch_tensor(t(MU[, o]))
      Sigma.L <- linalg_cholesky(torch_tensor(Sigma)$permute(c(3, 1, 2)))$contiguous()$requires_grad_(T)
      Mu <- torch_tensor(t(Mu), requires_grad = T)
      a.mask <- torch_tensor(t(D) != 0)
      a <- torch_tensor(a)$masked_select(a.mask)$requires_grad_(T)
      b <- torch_tensor(b, requires_grad = T)
      gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(dim(gamma)[1], a.mask)))
      gamma <- torch_cat(list(torch_zeros_like(gamma.mask[1:1]), torch_tensor(gamma)$transpose(2, 3)))$masked_select(gamma.mask)$requires_grad_(T)
      beta.mask <- torch_cat(list(torch_zeros_like(beta[, 1, drop = F]), torch_ones_like(beta)), 2)$t()$bool()  
      beta <- torch_cat(list(torch_zeros_like(beta.mask[1:1]), torch_tensor(t(beta))))$masked_select(beta.mask)$requires_grad_(T)
    })
  }
  
  parameters <- function() {
    list(Sigma = as_array(Sigma.L$matmul(Sigma.L$transpose(2, 3))$permute(c(2, 3, 1))),
         Mu = as_array(Mu$t()), a = as_array(unzip(a, a.mask)), b = as_array(b),
         gamma = as_array(unzip(gamma, gamma.mask)), beta = as_array(unzip(beta, beta.mask)))
  }
  
  zip <- function(Z, z, z.mask) {
    z$set_data(Z$masked_select(z.mask))
    z$detach_()$requires_grad_(T)
  }
  
  unzip <- function(z, z.mask) {
    torch_zeros(z.mask$shape)$masked_scatter_(z.mask, z)
  }
  
  
  objective <- function(a, gamma, beta) {
    z <- ((a + gamma)[X]$unsqueeze(2) * theta$unsqueeze(3))$sum(4) - (b - beta)[X]$unsqueeze(2)
    log.w <- torch_where(Y, nnf_logsigmoid(z), nnf_logsigmoid(-z))$sum(3) +
      torch_cat(lapply(1:K, function(k) {
        distr_multivariate_normal(Mu[k], scale_tril = Sigma.L[k])$log_prob(theta[X == k])
      })) - theta.logd
    (log.w$view(c(N, S, M))$logsumexp(3) - log(M))$mean(2)$sum()
  }
  
  update.l1 <- function() {
    prox <- function(x, state, params) {
      psi <- state$max_exp_avg_sq$sqrt() / (sqrt(1 - params$betas[2] ^ state$step) + params$eps)
      psi.max <- psi$max()
      phi <- psi / psi.max
      eta <- lambda * params$lr / psi.max
      z <- x
      while (T) {
        y <- z - phi * (z - x)
        z.new <- y$sign() * (y$abs() - eta)$maximum(0)
        if (as_array((z.new - z)$norm(Inf) < eps))
          break
        else
          z <- z.new
      }
      x$set_data(z.new)
    }
    
    # prox <- function(x, state, params) {
    #   betat <- params$beta ^ state$step
    #   eta <- params$lr / (1 - betat[1]) / (sqrt(state$max_exp_avg_sq / (1 - betat[2])) + params$eps)
    #   x$set_data(x$sign() * (x$abs() - lambda * eta)$maximum(0))
    # }

    Q <- objective(unzip(a, a.mask), unzip(gamma, gamma.mask), unzip(beta, beta.mask))
    opt$zero_grad()
    (-Q)$backward()
    opt$step()
    opt.state <- opt$state_dict()$state[5:6]
    with_no_grad({
      prox(gamma, opt.state[[1]], opt$param_groups[[1]])
      prox(beta, opt.state[[2]], opt$param_groups[[1]])
    })
  }
  
  update.mask <- function() {
    with_no_grad({
      gamma.Mask <- gamma == 0
      beta.Mask <- beta == 0
    })
    Q <- objective(unzip(a, a.mask), unzip(gamma, gamma.mask), unzip(beta, beta.mask))
    opt$zero_grad()
    (-Q)$backward()
    with_no_grad({
      gamma$grad$masked_fill_(gamma.Mask, 0)
      beta$grad$masked_fill_(beta.Mask, 0)
    })
    opt$step()
    with_no_grad({
      A <- unzip(a, a.mask)
      gamma$masked_fill_(gamma.Mask, 0)
      Gamma <- unzip(gamma, gamma.mask)
      beta$masked_fill_(beta.Mask, 0)
      Beta <- unzip(beta, beta.mask)
      mu <- Mu[1]$clone()
      MU$sub_(mu)
      Mu$sub_(mu)
      b$sub_(A$matmul(mu))
      Beta$add_(Gamma$matmul(mu))
      rescale <- Sigma.L[1]$norm(2, 2)
      A$mul_(rescale)
      Gamma$mul_(rescale)
      rescale.inv <- (1 / rescale)$diag()
      Sigma.L$set_data(rescale.inv$matmul(Sigma.L))
      SIGMA.L$set_data(rescale.inv$matmul(SIGMA.L))
      zip(A, a, a.mask)
      zip(Gamma, gamma, gamma.mask)
      zip(Beta, beta, beta.mask)
    })
  }
  
  init()
  opt <- optim_adam(list(Sigma.L, Mu, a, b, gamma, beta), lr, amsgrad = T)
  params <- list()
  params.old <- parameters()
  for (i in 1:iter) {
    sigma.L <- aperm(as_array(SIGMA.L), c(2, 3, 1))
    mu <- t(as_array(MU))
    theta <- sapply(1:N, function(n) {
      rmvn(S * M, mu[, n], sigma.L[, , n], isChol = T)
    }, simplify = 'array')
    theta.logd <- torch_stack(lapply(1:N, function(n) {
      dmvn(theta[, , n], mu[, n], sigma.L[, , n], T, isChol = T)
    }))
    theta <- torch_tensor(theta)$permute(c(3, 1, 2))
    
    update.l1()
    update.mask()
    
    params[[(i - 1) %% R + 1]] <- parameters()
    if (i %% R == 0) {
      params.new <- lapply(Reduce(params, f = function(x, y) {
        mapply(`+`, x, y)
      }), function(x) {
        x / R
      })
      dis <- mapply(params.old, params.new, FUN = function(x, y) {
        max(abs(x - y))
      })
      print(round(t(params.new$a), 3))
      print(round(params.new$beta, 3))
      print(dis)
      params.old <- params.new
      if (all(dis < 1e-2))
        break
    }
  }
  params.old
}

#self <- model$parameters
#lambda <- 15
# lambda <- 10
# tmp <- IW(Y, A, G, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, lambda)
