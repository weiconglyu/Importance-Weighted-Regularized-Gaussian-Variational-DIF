rm(list = ls())
library(psych)
library(MASS)
library(mirt) 
mirtCluster(4)
library(testit)
library(Matrix)
library(gtools)
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
source('GVEM_DIF_Simulation_Functions.R')

lbd.vec <- seq(0.4, 0.9, by = 0.05)
# constant used in calculating GIC. 
cons <- 1

#######################
#    Simulation 1
#######################

###################################
# Unif DIF, test length 20
###################################

# Study 1 four conditions 
condition=0
Nvec=c(500,1000)
J=20 #test length
for (N_condition in 1:2){
  N1=N2=N3=Nvec[N_condition] # sample size each group
  N=N1+N2+N3
  for (DIF_condition in 1:2){
    # read in MBR data
    params=read.csv(paste0("../Para",DIF_condition,".csv"),row.names = 1)
    condition=condition+1
    resp.all=read.csv(paste0("../RESP",condition,".csv"),row.names = 1)
    
    reps=20 # number of replications each condition
    indic=matrix(0,2,20);indic[1,1]=1;indic[2,2]=1;indic[1,3:11]=1;indic[2,12:20]=1 #loading structure indicator
    Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3)) # group indicator
    Unif=T # uniform DIF 
    
    for (rep in 1:reps){
      # read in data for each replication
      u=resp.all[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
      u=as.matrix(u)
      resp=u
      # estimate starting values
      Gind=Group
      Gind=as.numeric(as.factor(Gind))
      y=length(unique(Group)) 
      y.allgroup=rbind(rep(0,y-1),diag(y-1)) 
      N.vec=as.vector(table(Group))
      person=nrow(resp)
      item=ncol(resp)
      eta=matrix(0,person,item)
      eps=matrix(0,person,item)
      init=DIF_init(resp=as.data.frame(resp),Group=Group,indic=indic,Unif=Unif)
      G=init$G
      y.allgroup=rbind(rep(0,y-1),diag(y-1)) 
      Groupind=matrix(0,person,y-1)
      for (yy in 1:y){
        vec=which(Group==sort(unique(Group))[yy])
        for (i in 1:length(vec)){
          Groupind[vec[i],]=y.allgroup[yy,]
        }
      }
      r=init$r
      m=2 #fixed 2pl
      domain=init$r
      gra00=init$gra00
      grd00=init$grd00
      grbeta00=init$grbeta00
      grgamma00=init$grgamma00
      grbetamle=init$grbeta00
      grgammamle=init$grgamma00
      Mu.list=init$Mu0
      Sig.list=init$Sigma0
      
      eta=matrix(0,person,item)
      eps=matrix(0,person,item)
      for (yy in 1:y){
        Ny=N.vec[yy]
        Ny0=sum(N.vec[0:(yy-1)])
        uy=resp[which(Group==(unique(Group))[yy]),]
        initelse=init.else(uy,indic,gra00,-grd00)
        eta[(Ny0+1):(Ny0+Ny),]=initelse$eta0
        eps[(Ny0+1):(Ny0+Ny),]=initelse$eps0
      }
      
      #lbd.vec=seq(10,100,10)
      results <- lapply(lbd.vec, function(lbd) {
        print(system.time(result <- Reg_VEMM_DIF(resp,indic,lbd,eta,eps,Group,Unif,domain,y,Groupind,N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list,cons)))
        result
      })
      
      aic <- sapply(results, `[[`, 'AIC')
      bic <- sapply(results, `[[`, 'BIC')
      gic <- sapply(results, `[[`, 'GIC')
      result.aic <- results[[which.min(aic)]]
      result.bic <- results[[which.min(bic)]]
      result.gic <- results[[which.min(gic)]]
      print(c(result.aic$lbd, result.bic$lbd, result.gic$lbd))
      print(round(result.aic$beta, 3))
      print(round(result.bic$beta, 3))
      print(round(result.gic$beta, 3))
      save(results, aic, bic, gic, result.aic, result.bic, result.gic, file = paste0('GVEM_Sim1_Condition', condition, '_', rep, '.RData'))
    }
  }
}



# Simulation 2
# Non-Unif DIF on slope and intercept
condition=4
Nvec=c(500,1000)
J=20 #test length
for (N_condition in 1:2){
  N1=N2=N3=Nvec[N_condition] # sample size each group
  N=N1+N2+N3
  for (DIF_condition in 3:4){
    # read in MBR data
    params=read.csv(paste0("../Para",DIF_condition,"new.csv"),row.names = 1)
    condition=condition+1
    resp.all=read.csv(paste0("../RESP",condition,"new.csv"),row.names = 1)
    
    reps=20 # number of replications each condition
    indic=matrix(0,2,20);indic[1,1]=1;indic[2,2]=1;indic[1,3:11]=1;indic[2,12:20]=1 #loading structure indicator
    Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3)) # group indicator
    Unif=F # uniform DIF 
    
    for (rep in 1:reps){
      # read in data for each replication
      u=resp.all[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
      u=as.matrix(u)
      resp=u
      # estimate starting values
      Gind=Group
      Gind=as.numeric(as.factor(Gind))
      y=length(unique(Group)) 
      y.allgroup=rbind(rep(0,y-1),diag(y-1)) 
      N.vec=as.vector(table(Group))
      person=nrow(resp)
      item=ncol(resp)
      eta=matrix(0,person,item)
      eps=matrix(0,person,item)
      init=DIF_init(resp=as.data.frame(resp),Group=Group,indic=indic,Unif=Unif)
      G=init$G
      y.allgroup=rbind(rep(0,y-1),diag(y-1)) 
      Groupind=matrix(0,person,y-1)
      for (yy in 1:y){
        vec=which(Group==sort(unique(Group))[yy])
        for (i in 1:length(vec)){
          Groupind[vec[i],]=y.allgroup[yy,]
        }
      }
      r=init$r
      m=2 #fixed 2pl
      domain=init$r
      gra00=init$gra00
      grd00=init$grd00
      grbeta00=init$grbeta00
      grgamma00=init$grgamma00
      grbetamle=init$grbeta00
      grgammamle=init$grgamma00
      Mu.list=init$Mu0
      Sig.list=init$Sigma0
      
      eta=matrix(0,person,item)
      eps=matrix(0,person,item)
      for (yy in 1:y){
        Ny=N.vec[yy]
        Ny0=sum(N.vec[0:(yy-1)])
        uy=resp[which(Group==(unique(Group))[yy]),]
        initelse=init.else(uy,indic,gra00,-grd00)
        eta[(Ny0+1):(Ny0+Ny),]=initelse$eta0
        eps[(Ny0+1):(Ny0+Ny),]=initelse$eps0
      }

      #lbd.vec=seq(20,48,4)
      results <- lapply(lbd.vec, function(lbd) {
        print(system.time(result <- Reg_VEMM_DIF(resp,indic,lbd,eta,eps,Group,Unif,domain,y,Groupind,N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list,cons)))
        result
      })
      
      aic <- sapply(results, `[[`, 'AIC')
      bic <- sapply(results, `[[`, 'BIC')
      gic <- sapply(results, `[[`, 'GIC')
      result.aic <- results[[which.min(aic)]]
      result.bic <- results[[which.min(bic)]]
      result.gic <- results[[which.min(gic)]]
      print(c(result.aic$lbd, result.bic$lbd, result.gic$lbd))
      print(round(result.aic$beta, 3))
      print(round(result.bic$beta, 3))
      print(round(result.gic$beta, 3))
      save(results, aic, bic, gic, result.aic, result.bic, result.gic, file = paste0('GVEM_Sim2_Condition', condition, '_', rep, '.RData'))
    }
  }
}


# Simulation 3
# Non-Unif DIF only on slope
condition=4
Nvec=c(500,1000)
J=20 #test length
for (N_condition in 1:2){
  N1=N2=N3=Nvec[N_condition] # sample size each group
  N=N1+N2+N3
  for (DIF_condition in 3:4){
    # read in MBR data
    params=read.csv(paste0("../Para",DIF_condition,".csv"),row.names = 1)
    condition=condition+1
    resp.all=read.csv(paste0("../RESP",condition,".csv"),row.names = 1)
    
    reps=20 # number of replications each constion
    indic=matrix(0,2,20);indic[1,1]=1;indic[2,2]=1;indic[1,3:11]=1;indic[2,12:20]=1 #loading structure indicator
    Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3)) # group indicator
    Unif=F # uniform DIF 
    
    for (rep in 1:reps){
      # read in data for each replication
      u=resp.all[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
      u=as.matrix(u)
      resp=u
      # estimate starting values
      Gind=Group
      Gind=as.numeric(as.factor(Gind))
      y=length(unique(Group)) 
      y.allgroup=rbind(rep(0,y-1),diag(y-1)) 
      N.vec=as.vector(table(Group))
      person=nrow(resp)
      item=ncol(resp)
      eta=matrix(0,person,item)
      eps=matrix(0,person,item)
      init=DIF_init(resp=as.data.frame(resp),Group=Group,indic=indic,Unif=Unif)
      G=init$G
      y.allgroup=rbind(rep(0,y-1),diag(y-1)) 
      Groupind=matrix(0,person,y-1)
      for (yy in 1:y){
        vec=which(Group==sort(unique(Group))[yy])
        for (i in 1:length(vec)){
          Groupind[vec[i],]=y.allgroup[yy,]
        }
      }
      r=init$r
      m=2 #fixed 2pl
      domain=init$r
      gra00=init$gra00
      grd00=init$grd00
      grbeta00=init$grbeta00
      grgamma00=init$grgamma00
      grbetamle=init$grbeta00
      grgammamle=init$grgamma00
      Mu.list=init$Mu0
      Sig.list=init$Sigma0
      
      eta=matrix(0,person,item)
      eps=matrix(0,person,item)
      for (yy in 1:y){
        Ny=N.vec[yy]
        Ny0=sum(N.vec[0:(yy-1)])
        uy=resp[which(Group==(unique(Group))[yy]),]
        initelse=init.else(uy,indic,gra00,-grd00)
        eta[(Ny0+1):(Ny0+Ny),]=initelse$eta0
        eps[(Ny0+1):(Ny0+Ny),]=initelse$eps0
      }

      #lbd.vec=seq(10,50,5)
      results <- lapply(lbd.vec, function(lbd) {
        print(system.time(result <- Reg_VEMM_DIF(resp,indic,lbd,eta,eps,Group,Unif,domain,y,Groupind,N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list,cons)))
        result
      })
      
      aic <- sapply(results, `[[`, 'AIC')
      bic <- sapply(results, `[[`, 'BIC')
      gic <- sapply(results, `[[`, 'GIC')
      result.aic <- results[[which.min(aic)]]
      result.bic <- results[[which.min(bic)]]
      result.gic <- results[[which.min(gic)]]
      print(c(result.aic$lbd, result.bic$lbd, result.gic$lbd))
      print(round(result.aic$beta, 3))
      print(round(result.bic$beta, 3))
      print(round(result.gic$beta, 3))
      save(results, aic, bic, gic, result.aic, result.bic, result.gic, file = paste0('GVEM_Sim3_Condition', condition, '_', rep, '.RData'))
    }
  }
}
