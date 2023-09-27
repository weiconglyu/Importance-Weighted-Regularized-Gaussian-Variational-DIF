library(torch)
library(mvnfast)

# rm(list = ls())
# load('iw.RData')
# for (var in c('SIGMA', 'MU', 'Sigma', 'Mu', 'a', 'b', 'gamma', 'beta'))
#   assign(var, eval(parse(text = paste0('torch_tensor(', var, ')'))))

`%*%.torch_tensor` <- function(e1, e2) {
  torch_matmul(e1, e2)
}

`%<-%` <- function(e1, e2) {
  e1$set_data(e2)
}

t.torch_tensor <- function(e) {
  torch_transpose(e, -1, -2)
}

prox <- function(x, lambda) {
  sign(x) * (abs(x) - lambda)$maximum(0)
}


IW <- function(Y, D, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, lambda, S = 10, M = 10, lr = 0.05, iter = 100, eps = 0.05) {
  # for (var in c('SIGMA', 'MU', 'Sigma', 'Mu', 'a', 'b', 'gamma', 'beta'))
  #   assign(var, eval(parse(text = paste0('as.array(', var, ')'))))
  # save(list = ls(), file = 'iw.RData')
  # stop()
  
  init <- function() {
    with(parent.frame(), {
      N <- nrow(Y)
      J <- ncol(Y)
      K <- max(X)
      Y <- torch_tensor(Y)$unsqueeze(2)$bool()
      SIGMA.L <- linalg_cholesky(SIGMA)$contiguous()
      Sigma.L <- linalg_cholesky(Sigma)$contiguous()$requires_grad_(T)
      Mu <- Mu$clone()$requires_grad_(T)
      a.mask <- torch_tensor(t(D) != 0)
      a <- a$masked_select(a.mask)$requires_grad_(T)
      a.mask.diag <- a.mask$diag_embed()$type_as(a)
      b <- b$clone()$requires_grad_(T)
      gamma.mask <- torch_stack(c(torch_zeros_like(a.mask), replicate(gamma$shape[1] - 1, a.mask)))
      gamma <- gamma$masked_select(gamma.mask)$requires_grad_(T)
      beta.mask <- torch_cat(list(torch_zeros_like(beta[1:1]), torch_ones_like(beta[2:beta$shape[1]])))$bool()  
      beta <- beta$masked_select(beta.mask)$requires_grad_(T)
    })
  }
  
  parameters <- function() {
    list(Sigma = Sigma.L %*% t(Sigma.L), Sigma.L = Sigma.L, Mu = Mu,
         a = unzip(a, a.mask), b = b, gamma = unzip(gamma, gamma.mask), beta = unzip(beta, beta.mask))
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
    log.w$view(c(N, S, M))$logsumexp(3)$mean(2)$sum()
  }
  
  update.l1 <- function() {
    adaprox <- function(x, state, params, iter = 100, eps = 1e-3) {
      psi <- sqrt(state$max_exp_avg_sq / (1 - params$betas[2] ^ state$step)) + params$eps
      psi.max <- max(psi)
      Psi <- psi / psi.max
      r.lambda <- params$lr / psi.max * lambda
      z <- x
      for (i in 1:iter) {
        z.new <- prox(z - Psi * (z - x), r.lambda)
        if (as.array((z.new - z)$norm(Inf)) < eps)
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
      adaprox(gamma, opt.state[[1]], opt$param_groups[[1]])
      adaprox(beta, opt.state[[2]], opt$param_groups[[1]])
    })
  }
  
  update.mask <- function() {
    with_no_grad({
      gamma.zero <- gamma == 0
      beta.zero <- beta == 0
    })
    Q <- objective(unzip(a, a.mask), unzip(gamma, gamma.mask), unzip(beta, beta.mask))
    opt$zero_grad()
    (-Q)$backward()
    with_no_grad({
      gamma$grad$masked_fill_(gamma.zero, 0)
      beta$grad$masked_fill_(beta.zero, 0)
    })
    opt$step()
    with_no_grad({
      A <- unzip(a, a.mask)
      gamma$masked_fill_(gamma.zero, 0)
      Gamma <- unzip(gamma, gamma.mask)
      beta$masked_fill_(beta.zero, 0)
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
      SIGMA.L$set_data(rescale.inv %*% SIGMA.L)
      Sigma.L$set_data(rescale.inv %*% Sigma.L)
      zip(A, a, a.mask)
      zip(Gamma, gamma, gamma.mask)
      zip(Beta, beta, beta.mask)
      
      # aG <- (A + Gamma)$unsqueeze(4)
      # aG.t <- aG$transpose(3, 4)
      # AG <- aG[X]
      # AG.t <- AG$transpose(3, 4)
      # BB <- (b - Beta)[X]
      # mu <- MU$unsqueeze_(2)
      # sigma.mu <- a.mask.diag$transpose(2, 3)$matmul((SIGMA + mu$transpose(2, 3)$matmul(mu))$unsqueeze(2))$matmul(a.mask.diag)
      # xi <- (BB$square() - 2 * BB * AG.t$matmul(mu$unsqueeze(4))$view(c(N, -1)) + AG.t$matmul(sigma.mu)$matmul(AG)$view(c(N, -1)))$sqrt()$mul((AG.t$matmul(mu$unsqueeze(4))$view(c(N, -1)) - BB)$sign())
      # eta <- torch_where(abs(xi) < 1e-3, 0.125, (1 / (1 + exp(-xi)) - 0.5) / (2 * xi))
      # Sigma.inv <- Sigma$inverse()
      # SIGMA.inv <- Sigma.inv[X] + 2 * (eta$view(c(N, -1, 1, 1)) * aG$matmul(aG.t)[X])$sum(2)
      # SIGMA$set_data(SIGMA.inv$inverse())
      # MU$set_data(SIGMA$matmul(((Y$type_as(eta)$squeeze(2) - 0.5 + 2 * eta * BB)$view(c(N, -1, 1, 1)) * AG)$sum(2) + Sigma.inv$matmul(Mu$unsqueeze(3))[X])$squeeze(3))
    })
  }
  
  init()
  opt <- optim_adam(list(Sigma.L, Mu, a, b, gamma, beta), lr, amsgrad = T)
  params <- parameters()
  for (i in 1:iter) {
    sigma.L <- as.array(SIGMA.L$permute(c(2, 3, 1)))
    mu <- as.array(t(MU))
    theta <- sapply(1:N, function(n) {
      rmvn(S * M, mu[, n], sigma.L[, , n], isChol = T)
    }, simplify = 'array')
    theta.logd <- torch_stack(lapply(1:N, function(n) {
      dmvn(theta[, , n], mu[, n], sigma.L[, , n], T, isChol = T)
    }))
    theta <- torch_tensor(aperm(theta, c(3, 1, 2)))
    update.l1()
    update.mask()
    
    params.new <- parameters()
    ok <- mapply(params[-1], params.new[-1], FUN = function(x, y) {
      as.array((abs(x - y) <= abs(x) * eps)$all())
    })
    # print(round(as.array(t(params.new$a)), 3))
    # print(round(as.array(params.new$beta), 3))
    # print(ok)
    if (all(ok))
      break
    else
      params <- params.new
  }
  c(n = i, params.new)
}

#self <- model$parameters
#lambda <- 15
# lambda <- 10
# load('iw.RData')
# for (var in c('SIGMA', 'MU', 'Sigma', 'Mu', 'a', 'b', 'gamma', 'beta')) {
#   assign(var, eval(parse(text = paste0('torch_tensor(', var, ')'))))
# }
# eps <- 1e-3
# tmp <- IW(Y, D, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, lambda)

# load('iw.RData')
# SIGMA<-torch_tensor(SIGMA)
# Sigma<-torch_tensor(Sigma)
# MU<-torch_tensor(MU)
# Mu<-torch_tensor(Mu)
# a<-torch_tensor(a)
# b<-torch_tensor(b)
# gamma<-torch_tensor(gamma)
# beta<-torch_tensor(beta)
# IW(Y, D, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, lambda)

#tmp <- IW(Y, D, X, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, lambda)
