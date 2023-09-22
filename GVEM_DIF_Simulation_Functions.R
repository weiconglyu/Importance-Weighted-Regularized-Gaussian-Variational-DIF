source('GVEM_DIF_EM.R')
source('GVEM_DIF_IW.R')

DIF_init=function(resp,Group,indic,Unif){
  m=2 ##fixed, 2pl only
  N=nrow(resp)
  J=ncol(resp)
  domain=nrow(indic)
  y=length(unique(Group)) 
  y.allgroup=rbind(rep(0,y-1),diag(y-1)) 
  G=matrix(0,N,y-1)
  for (yy in 1:y){
    vec=which(Group==sort(unique(Group))[yy])
    for (i in 1:length(vec)){
      G[vec[i],]=y.allgroup[yy,]
    }
  }
  # defalt for no impact (when using mirt to estimate MLE, fix the mean and variance for all groups)
  COV <- matrix(TRUE,domain,domain); diag(COV)=FALSE
  model <- mirt.model(t(indic), COV=COV) ##
  if (Unif==T){
    md.noncons0 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('slopes'))
    starting5new=cbind(coef(md.noncons0,simplify=T)[[1]]$items[,1:(domain+m-1)])
    for (yy in 2:y){
      starting5new=cbind(starting5new,coef(md.noncons0,simplify=T)[[yy]]$items[,domain+1]-coef(md.noncons0,simplify=T)[[1]]$items[,domain+1])
      
    }
    gra00=as.matrix(starting5new[,1:domain])
    rownames(gra00) <- c()
    colnames(gra00) <- c()
    grd00=matrix(starting5new[,domain+1],J,1)
    grgamma00=array(0,dim=c((y-1),domain,J))
    grbeta00=as.matrix(starting5new[,(domain+1+1):(domain+1+1+y-1-1)])
    rownames(grbeta00) <- c()
    colnames(grbeta00) <- c()
    #Sigma0=array(double(domain*domain*y),dim = c(domain,domain,y))
    Sigma0=matrix(0,domain*y,domain)
    #Mu0 = matrix(0,domain,y)
    Mu0 = numeric(domain*y)
    for (yy in 1:y){
      Sigma0[((yy-1)*domain+1):(yy*domain),]=coef(md.noncons0,simplify=T)[[yy]]$cov
      #Mu0[,yy]=coef(md.noncons0,simplify=T)[[yy]]$means
    }
  } else {
    md.noncons0 <- multipleGroup(resp, model, group = Group,SE=TRUE)
    starting5new=cbind(coef(md.noncons0,simplify=T)[[1]]$items[,1:(domain+m-1)])
    for (yy in 2:y){
      starting5new=cbind(starting5new,coef(md.noncons0,simplify=T)[[yy]]$items[,1:domain]-coef(md.noncons0,simplify=T)[[1]]$items[,1:domain])
    }
    for (yy in 2:y){
      starting5new=cbind(starting5new,coef(md.noncons0,simplify=T)[[yy]]$items[,domain+1]-coef(md.noncons0,simplify=T)[[1]]$items[,domain+1])
      
    }
    gra00=as.matrix(starting5new[,1:domain])
    rownames(gra00) <- c()
    colnames(gra00) <- c()
    grd00=matrix(starting5new[,domain+1],J,1)
    grgamma00=array(0,dim=c((y-1),domain,J))
    for (yy in 2:y){
      grgamma00[(yy-1),,]=t(starting5new[,(domain+2+(yy-2)*domain):(2*domain+1+(yy-2)*domain)])
    }
    grbeta00=as.matrix(starting5new[,(2*domain+2+(y-2)*domain):(2*domain+1+(y-2)*domain+(y-1))])
    rownames(grbeta00) <- c()
    colnames(grbeta00) <- c()
    #Sigma0=array(double(domain*domain*y),dim = c(domain,domain,y))
    Sigma0=matrix(0,domain*y,domain)
    Mu0 = matrix(0,domain,y)
    for (yy in 1:y){
      Sigma0[((yy-1)*domain+1):(yy*domain),]=coef(md.noncons0,simplify=T)[[yy]]$cov
      #Mu0[,yy]=coef(md.noncons0,simplify=T)[[yy]]$means
    }
  }
  
  return(list(G=G,y=y,r=domain,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Sigma0=Sigma0,Mu0 =Mu0))
}

init.else<-function(u,indic,a0,b0){
  domain=nrow(indic)
  person=dim(u)[1]
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  #person*item
  xi=array(1,person)%*%t(b0)-theta%*%t(a0)
  eta0=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
  eta0[abs(xi) < 0.01] <- 0.125
  eps0=xi
  return(list(domain=domain,a0=a0,b0=b0,eta0=eta0,eps0=eps0))
}

Reg_VEMM_DIF <- function(resp,indic,lbd,eta,eps,Group,Unif,domain,y,G,N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list,cons)
{
  resp=as.matrix(resp)
  Gind=Group
  Gind=as.numeric(as.factor(Gind))
  y.allgroup=rbind(rep(0,y-1),diag(y-1)) 
  N.vec=as.vector(table(Group))
  person=nrow(resp)
  item=ncol(resp)
  # free slope and intercept
  new_a=gra00
  new_b=-grd00
  new_gamma=grgamma00
  new_beta=grbeta00
  eta=eta
  eps=eps
  Mu=matrix(0,domain,y)
  Sigma = array(0,c(domain,domain,y))
  for (yy in 1:y){
    Mu[,yy]=Mu.list[((yy-1)*domain+1):((yy-1)*domain+domain)]
    Sigma[,,yy]=Sig.list[((yy-1)*domain+1):((yy-1)*domain+domain),]
  }
  
  Gind <- as.integer(as.factor(Gind))
  o <- order(Gind)
  resp <- resp[o, ]
  Gind <- Gind[o]
  eta <- eta[o, ]

  lambda <- sqrt(person) * lbd * 5
  list2env(lapply(EM(resp, indic, Gind, eta, Sigma, Mu, new_a, new_b, new_gamma, new_beta, lambda), as.array), environment())
  list2env(as.list(IC(resp, Gind, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, cons)), environment())

  varlist <- function(...) {
    setNames(lapply(list(...), as.array), as.character(substitute(alist(...)))[-1])
  }
  varlist(lbd, lambda, SIGMA, MU, Sigma, Mu, a, b, gamma, beta, n, ll, BIC, GIC)
}
