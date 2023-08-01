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

#######################
#    Simulation 1
#######################

###################################
# Unif DIF, test length 20
###################################

setwd("/Users/ruoyizhu/Documents/GitHub/GVEM_DIF/GVEM_sim_results")
# Study 1 four conditions 
condition=0
Nvec=c(500,1000)
J=20 #test length
for (N_condition in 1:2){
  N1=N2=N3=Nvec[N_condition] # sample size each group
  N=N1+N2+N3
  for (DIF_condition in 1:2){
    # read in MBR data
    params=read.csv(paste0("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Para",DIF_condition,".csv"),row.names = 1)
    condition=condition+1
    resp.all=read.csv(paste0("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/RESP",condition,".csv"),row.names = 1)
    
    reps=20 # number of replications each constion
    #store results for all replications to compute power 
    lbd.1=lbd.2=numeric(reps) #eta selected in each replication
    Betas.1=Betas.2=array(double(J*2*reps),dim = c(J,2,reps)) # beta (dif on intercept) estimated in each replication
    ADmat.1=ADmat.2=array(double(J*3*reps),dim = c(J,3,reps)) # item parameters cbind(a,d) estimated in each replication
    means.1=means.2=array(double(1*6*reps),dim=c(1,6,reps)) # mean vec est each replication
    covs.1=covs.2=array(double(6*2*reps),dim=c(6,2,reps)) # covariance est each replication
    
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
      # constant used in calculating GIC. 
      cons=5
      
      lbd.vec=seq(20,42,3)
      #store results for all tuning parameter values
      lls=bics=gics=rep(0,length(lbd.vec))
      ADmat=array(double(J*3*length(lbd.vec)),dim = c(J,3,length(lbd.vec)))
      Betas=array(double(J*2*length(lbd.vec)),dim = c(J,2,length(lbd.vec)))
      means=array(double(1*6*length(lbd.vec)),dim=c(1,6,length(lbd.vec)))
      covs=array(double(6*2*length(lbd.vec)),dim=c(6,2,length(lbd.vec)))
      for (k in 1:length(lbd.vec))
      {
        lbd=lbd.vec[k]
        ptm <- proc.time()
        sim=Reg_VEMM_DIF(resp=resp,indic=indic,lbd=lbd,eta=eta,eps=eps,Group=Group,Unif=Unif,domain=domain,y,G=Groupind,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list,cons=cons)
        print(proc.time() - ptm)
        bics[k]=sim$bic
        gics[k]=sim$gic
        lls[k]=sim$Likelihood
        ADmat[,,k]=sim$est
        Betas[,,k]=sim$Beta
        means[,,k]=sim$means
        covs[,,k]=sim$Covs
      }
      
      kk=which.min(bics)
      lbd.1[rep]=lbd.vec[kk]
      Betas.1[,,rep]=Betas[,,kk]
      ADmat.1[,,rep]=ADmat[,,kk]
      means.1[,,rep]=means[,,kk]
      covs.1[,,rep]=covs[,,kk]
      
      
      kk2=which.min(gics)
      lbd.2[rep]=lbd.vec[kk2]
      Betas.2[,,rep]=Betas[,,kk2]
      ADmat.2[,,rep]=ADmat[,,kk2]
      means.2[,,rep]=means[,,kk2]
      covs.2[,,rep]=covs[,,kk2]
      print(lbd.2[rep])
      write.csv(lbd.1[rep],file = paste0("GVEM_Sim1_Condition",condition,"_lbd1_",rep))
      write.csv(ADmat.1[,,rep],file = paste0("GVEM_Sim1_Condition",condition,"_ADmat1_",rep))
      write.csv(Betas.1[,,rep],file = paste0("GVEM_Sim1_Condition",condition,"_Beta1_",rep))
      write.csv(means.1[,,rep],file = paste0("GVEM_Sim1_Condition",condition,"_Mean1_",rep))
      write.csv(covs.1[,,rep],file = paste0("GVEM_Sim1_Condition",condition,"_Cov1_",rep))
      
      write.csv(lbd.2[rep],file = paste0("GVEM_Sim1_Condition",condition,"_lbd2_",rep))
      write.csv(ADmat.2[,,rep],file = paste0("GVEM_Sim1_Condition",condition,"_ADmat2_",rep))
      write.csv(Betas.2[,,rep],file = paste0("GVEM_Sim1_Condition",condition,"_Beta2_",rep))
      write.csv(means.2[,,rep],file = paste0("GVEM_Sim1_Condition",condition,"_Mean2_",rep))
      write.csv(covs.2[,,rep],file = paste0("GVEM_Sim1_Condition",condition,"_Cov2_",rep))
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
    params=read.csv(paste0("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Para",DIF_condition,"new.csv"),row.names = 1)
    condition=condition+1
    resp.all=read.csv(paste0("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/RESP",condition,"new.csv"),row.names = 1)
    
    reps=25 # number of replications each constion
    #store results for all replications to compute power 
    lbd.1=lbd.2=numeric(reps) #eta selected in each replication
    Gammas.1=Gammas.2=array(double(2*J*m*reps),dim = c(2,2,J,reps))
    Betas.1=Betas.2=array(double(J*2*reps),dim = c(J,2,reps)) # beta (dif on intercept) estimated in each replication
    ADmat.1=ADmat.2=array(double(J*3*reps),dim = c(J,3,reps)) # item parameters cbind(a,d) estimated in each replication
    means.1=means.2=array(double(1*6*reps),dim=c(1,6,reps)) # mean vec est each replication
    covs.1=covs.2=array(double(6*2*reps),dim=c(6,2,reps)) # covariance est each replication
    
    indic=matrix(0,2,20);indic[1,1]=1;indic[2,2]=1;indic[1,3:11]=1;indic[2,12:20]=1 #loading structure indicator
    Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3)) # group indicator
    Unif=F # uniform DIF 
    
    for (rep in 1:25){
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
      # constant used in calculating GIC. 
      cons=0.1
      
      lbd.vec=seq(24,60,4)
      #store results for all tuning parameter values
      lls=lls2=lls3=bics=gics=rep(0,length(lbd.vec))
      ADmat=array(double(J*3*length(lbd.vec)),dim = c(J,3,length(lbd.vec)))
      Gammas=array(double(2*J*m*length(lbd.vec)),dim = c(2,2,J,length(lbd.vec)))
      Betas=array(double(J*2*length(lbd.vec)),dim = c(J,2,length(lbd.vec)))
      means=array(double(1*6*length(lbd.vec)),dim=c(1,6,length(lbd.vec)))
      covs=array(double(6*2*length(lbd.vec)),dim=c(6,2,length(lbd.vec)))
      for (k in 1:length(lbd.vec))
      {
        lbd=lbd.vec[k]
        ptm <- proc.time()
        sim=Reg_VEMM_DIF(resp=resp,indic=indic,lbd=lbd,eta=eta,eps=eps,Group=Group,Unif=Unif,domain=domain,y,G=Groupind,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list,cons=cons)
        print(proc.time() - ptm)
        bics[k]=sim$bic
        gics[k]=sim$gic
        lls[k]=sim$Likelihood
        lls2[k]=sim$Likelihood2
        lls3[k]=sim$Likelihood3
        ADmat[,,k]=sim$est
        Gammas[,,,k]=sim$Gamma
        Betas[,,k]=sim$Beta
        means[,,k]=sim$means
        covs[,,k]=sim$Covs
      }
      
      #plot(lbd.vec,lls)
      #plot(lbd.vec,lls2)
      #plot(lbd.vec,lls3)
      
      (kk=which.min(bics))
      lbd.1[rep]=lbd.vec[kk]
      Gammas.1[,,,rep]=Gammas[,,,kk]
      (Betas.1[,,rep]=Betas[,,kk])
      ADmat.1[,,rep]=ADmat[,,kk]
      means.1[,,rep]=means[,,kk]
      covs.1[,,rep]=covs[,,kk]
      
      
      (kk2=which.min(gics))
      lbd.2[rep]=lbd.vec[kk2]
      Gammas.2[,,,rep]=Gammas[,,,kk2]
      (Betas.2[,,rep]=Betas[,,kk2])
      ADmat.2[,,rep]=ADmat[,,kk2]
      means.2[,,rep]=means[,,kk2]
      covs.2[,,rep]=covs[,,kk2]
      print(lbd.2[rep])
      write.csv(lbd.1[rep],file = paste0("GVEM_Sim2_Condition",condition,"_lbd1_",rep))
      write.csv(ADmat.1[,,rep],file = paste0("GVEM_Sim2_Condition",condition,"_ADmat1_",rep))
      write.csv(Gammas.1[,,,rep],file =  paste0("GVEM_Sim2_Condition",condition,"_Gamma1_",rep))
      write.csv(Betas.1[,,rep],file = paste0("GVEM_Sim2_Condition",condition,"_Beta1_",rep))
      write.csv(means.1[,,rep],file = paste0("GVEM_Sim2_Condition",condition,"_Mean1_",rep))
      write.csv(covs.1[,,rep],file = paste0("GVEM_Sim2_Condition",condition,"_Cov1_",rep))
      
      write.csv(lbd.2[rep],file = paste0("GVEM_Sim2_Condition",condition,"_lbd2_",rep))
      write.csv(ADmat.2[,,rep],file = paste0("GVEM_Sim2_Condition",condition,"_ADmat2_",rep))
      write.csv(Gammas.1[,,,rep],file =  paste0("GVEM_Sim2_Condition",condition,"_Gamma2_",rep))
      write.csv(Betas.2[,,rep],file = paste0("GVEM_Sim2_Condition",condition,"_Beta2_",rep))
      write.csv(means.2[,,rep],file = paste0("GVEM_Sim2_Condition",condition,"_Mean2_",rep))
      write.csv(covs.2[,,rep],file = paste0("GVEM_Sim2_Condition",condition,"_Cov2_",rep))
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
    params=read.csv(paste0("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/Para",DIF_condition,".csv"),row.names = 1)
    condition=condition+1
    resp.all=read.csv(paste0("/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData/RESP",condition,".csv"),row.names = 1)
    
    reps=20 # number of replications each constion
    #store results for all replications to compute power 
    lbd.1=lbd.2=numeric(reps) #eta selected in each replication
    Gammas.1=Gammas.2=array(double(2*J*m*reps),dim = c(2,2,J,reps))
    Betas.1=Betas.2=array(double(J*2*reps),dim = c(J,2,reps)) # beta (dif on intercept) estimated in each replication
    ADmat.1=ADmat.2=array(double(J*3*reps),dim = c(J,3,reps)) # item parameters cbind(a,d) estimated in each replication
    means.1=means.2=array(double(1*6*reps),dim=c(1,6,reps)) # mean vec est each replication
    covs.1=covs.2=array(double(6*2*reps),dim=c(6,2,reps)) # covariance est each replication
    
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
      # constant used in calculating GIC. 
      cons=3
      
      lbd.vec=seq(10,50,5)
      #store results for all tuning parameter values
      bics=gics=rep(0,length(lbd.vec))
      ADmat=array(double(J*3*length(lbd.vec)),dim = c(J,3,length(lbd.vec)))
      Gammas=array(double(2*J*m*length(lbd.vec)),dim = c(2,2,J,length(lbd.vec)))
      Betas=array(double(J*2*length(lbd.vec)),dim = c(J,2,length(lbd.vec)))
      means=array(double(1*6*length(lbd.vec)),dim=c(1,6,length(lbd.vec)))
      covs=array(double(6*2*length(lbd.vec)),dim=c(6,2,length(lbd.vec)))
      for (k in 1:length(lbd.vec))
      {
        lbd=lbd.vec[k]
        ptm <- proc.time()
        sim=Reg_VEMM_DIF(resp=resp,indic=indic,lbd=lbd,eta=eta,eps=eps,Group=Group,Unif=Unif,domain=domain,y,G=Groupind,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list,cons=cons)
        print(proc.time() - ptm)
        bics[k]=sim$bic
        gics[k]=sim$gic
        ADmat[,,k]=sim$est
        Gammas[,,,k]=sim$Gamma
        Betas[,,k]=sim$Beta
        means[,,k]=sim$means
        covs[,,k]=sim$Covs
      }
      
      kk=which.min(bics)
      lbd.1[rep]=lbd.vec[kk]
      Gammas.1[,,,rep]=Gammas[,,,kk]
      Betas.1[,,rep]=Betas[,,kk]
      ADmat.1[,,rep]=ADmat[,,kk]
      means.1[,,rep]=means[,,kk]
      covs.1[,,rep]=covs[,,kk]
      
      
      kk2=which.min(gics)
      lbd.2[rep]=lbd.vec[kk2]
      Gammas.2[,,,rep]=Gammas[,,,kk2]
      Betas.2[,,rep]=Betas[,,kk2]
      ADmat.2[,,rep]=ADmat[,,kk2]
      means.2[,,rep]=means[,,kk2]
      covs.2[,,rep]=covs[,,kk2]
      print(lbd.2[rep])
      write.csv(lbd.1[rep],file = paste0("GVEM_Sim3_Condition",condition,"_lbd1_",rep))
      write.csv(ADmat.1[,,rep],file = paste0("GVEM_Sim3_Condition",condition,"_ADmat1_",rep))
      write.csv(Gammas.1[,,,rep],file =  paste0("GVEM_Sim3_Condition",condition,"_Gamma1_",rep))
      write.csv(Betas.1[,,rep],file = paste0("GVEM_Sim3_Condition",condition,"_Beta1_",rep))
      write.csv(means.1[,,rep],file = paste0("GVEM_Sim3_Condition",condition,"_Mean1_",rep))
      write.csv(covs.1[,,rep],file = paste0("GVEM_Sim3_Condition",condition,"_Cov1_",rep))
      
      write.csv(lbd.2[rep],file = paste0("GVEM_Sim3_Condition",condition,"_lbd2_",rep))
      write.csv(ADmat.2[,,rep],file = paste0("GVEM_Sim3_Condition",condition,"_ADmat2_",rep))
      write.csv(Gammas.1[,,,rep],file =  paste0("GVEM_Sim3_Condition",condition,"_Gamma2_",rep))
      write.csv(Betas.2[,,rep],file = paste0("GVEM_Sim3_Condition",condition,"_Beta2_",rep))
      write.csv(means.2[,,rep],file = paste0("GVEM_Sim3_Condition",condition,"_Mean2_",rep))
      write.csv(covs.2[,,rep],file = paste0("GVEM_Sim3_Condition",condition,"_Cov2_",rep))
    }
    
  }
}
