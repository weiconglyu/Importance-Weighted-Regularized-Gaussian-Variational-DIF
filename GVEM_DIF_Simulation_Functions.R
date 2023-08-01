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
#####cfa_2p functions###
#e step
e_cfa2pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List ecfa2(const arma::mat&u,const arma::vec& Mu, const arma::mat& Sigma, const int& domain, const int& person, 
const int& item, const arma::mat& eta, const arma::mat& new_a,const arma::vec& new_b){
arma::mat MU=mat(domain,person,fill::zeros);
arma::cube SIGMA=zeros(domain,domain,person);
arma::mat Spart=mat(domain,domain,fill::zeros);
arma::mat xi=mat(person,item,fill::zeros);
for(int i=0; i<person; ++i){
arma::mat sigma_part=arma::mat(domain,domain,arma::fill::zeros);
//mu_part.zeros();
arma::vec mu_part= arma::zeros(domain);
for(int j=0; j<item; ++j){
if(NumericVector::is_na(eta(i,j))){
continue;}
sigma_part=sigma_part+eta(i,j)*trans(new_a.row(j))*new_a.row(j);
mu_part=mu_part+trans((2*eta(i,j)*new_b(j)+u(i,j)-0.5)*new_a.row(j));
}
arma::mat sigmahat=solve((solve(Sigma,eye(domain,domain))+2*sigma_part),eye(domain,domain));
arma::vec muhat=sigmahat*(mu_part+(solve(Sigma,eye(domain,domain))*Mu));
SIGMA.slice(i)=sigmahat;
MU.col(i)=muhat;
arma::mat apro=new_a*(sigmahat+muhat*trans(muhat))*trans(new_a);
xi.row(i)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+apro.diag()));
Spart=Spart+sigmahat+(muhat-Mu)*trans((muhat-Mu));
}
return List::create(Named("Spart") = Spart,Named("SIGMA") = SIGMA,
Named("MU") = MU,Named("xi") = xi);
}
'
sourceCpp(code=e_cfa2pl)

#loglikelihood
#loglikelihood
llloglikelihood<-'
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace arma;
const double log2pi = std::log(2.0*M_PI);

// [[Rcpp::export]]
double dmvnrm2(arma::rowvec x, arma::rowvec mean, arma::mat sigma, bool logd =false){
  int xdim = x.n_elem;
  double out;
  arma::mat upper = inv(trimatu(chol(sigma))).t();
  double uppersum = sum(log(upper.diag()));
  double constants = -(static_cast<double>(xdim)/2.0)*log2pi;
  arma::vec z = upper*(x - mean).t();
  out = constants - 0.5*sum(z%z) + uppersum;
  if (logd == false){
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double loglikelihood(const arma::mat&resp,const arma::mat&new_a,const arma::mat&new_b,
const arma::mat&new_beta,arma::vec Aallgroups,arma::mat X2,arma::mat G,
int domain,int person,int Gs,int item)
{
double ll1 = 0; 
for(int i=0; i<person; i++){
double ll2 = 0; 
for(int g=0; g<Gs; g++){
double ll3 = 0; 
for(int j=0; j<item; j++){
double Pigj=1/(1+exp(-(as_scalar(( new_a.row(j))*(X2.row(g).t()))-(as_scalar(new_b.row(j))-as_scalar(new_beta.row(j)*(G.row(i).t()))))));
ll3 = ll3 + (resp(i,j)*log(Pigj)+(1-resp(i,j))*log(1-Pigj));
}
double fi=Aallgroups(i*Gs+g);
ll2=ll2+exp(ll3+log(fi))*pow(0.2,domain);
}
ll1=ll1+log(ll2);
}
return(ll1);
}
'
sourceCpp(code=llloglikelihood)




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
  eta0[is.na(u)]=NA
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
  
  
  SIGMA = array(double(domain*domain*person),dim = c(domain,domain,person))
  MU = matrix(0,domain,person)
  
  xi = matrix(0,person,item)
  rankSigma=domain
  converge=0
  singularity = 0
  n=0
  # EM cycle
  while(converge == 0  && rankSigma == domain)
  {
    #Step 1: Update mu and sigma & Step 2.1: Update xi and eta(xi). Meanwhile, updata Sigma. -2*eta[i,j]*(new_beta[,j]%*%G[i,])
    par_Sigma=Sigma 
    par_Mu=Mu 
    yy=1
    uy=resp[which(Gind==yy),]
    Muy=Mu[,yy]
    Sigmay=as.matrix(Sigma[,,yy])
    persony=N.vec[yy]
    etay=eta[which(Gind==yy),]
    new_ay=new_a #new_a[,j]+G[i,]%*%new_gamma[,,j]
    new_by=as.vector(new_b)#new_b[j]-new_beta[,j]
    rs1<-ecfa2(as.matrix(uy),Muy,Sigmay,domain, persony, item, etay, new_ay, new_by)
    SIGMA[,,which(Gind==yy)]=rs1$SIGMA
    MU[,which(Gind==yy)]=rs1$MU
    xi[which(Gind==yy),]=rs1$xi
    Sigma[,,yy]=rs1$Spart/N.vec[1]
    for (yy in 2:y){
      Ny=N.vec[yy]
      uy=resp[which(Gind==yy),]
      Muy=Mu[,yy]
      Sigmay=Sigma[,,yy]
      persony=N.vec[yy]
      etay=eta[which(Gind==yy),]
      new_ay=new_a+  t(new_gamma[yy-1,,]) #new_a[,j]+G[i,]%*%new_gamma[,,j]
      new_by=as.vector(new_b-new_beta[,yy-1])#new_b[j]-new_beta[,j]
      rs1<-ecfa2(as.matrix(uy),Muy,Sigmay,domain, persony, item, etay, new_ay, new_by)
      SIGMA[,,which(Gind==yy)]=rs1$SIGMA
      MU[,which(Gind==yy)]=rs1$MU
      xi[which(Gind==yy),]=rs1$xi
      Sigma[,,yy]=rs1$Spart/Ny 
      Mu[,yy]=rowMeans(rs1$MU)
    }
    d_temp1=sqrt(diag(diag(Sigma[,,1])))
    for (yy in 1:y){
      Sigma[,,yy]=solve(d_temp1)%*%Sigma[,,yy]%*%solve(d_temp1)
    }
    
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(eta)< 0.01,0.125)
    
    # Step 2.2: Update b.
    par_b=new_b
    b_part = matrix(0,person,item)
    for (i in 1:person){
      new_gamma_i=NULL
      for (j in 1:item){
        new_gamma_ij=G[i,]%*%new_gamma[,,j]
        new_gamma_i=cbind(new_gamma_i,t(new_gamma_ij))
      }
      b_part[i,]=(new_a+t(new_gamma_i))%*%MU[,i]+new_beta%*%G[i,] 
    }
    new_b=colSums(0.5-resp+2*eta*b_part)/colSums(2*eta)
    
    # Step 2.3: Update a.
    par_a=new_a
    for (j in 1:item){
      a_nu=0;a_de=0
      Ind=indic[,j]
      iind=which(Ind==1)
      for (i in 1:person){
        sigma=SIGMA[,,i]
        sigma=sigma[iind,iind]
        mu=MU[iind,i]
        a_de=a_de+eta[i,j]*sigma+eta[i,j]*(mu*t(mu))
        a_nu=a_nu+(resp[i,j]-0.5+2*(new_b[j]-new_beta[j,]%*%G[i,])*eta[i,j])*mu 
        -eta[i,j]*2*G[i,]%*%new_gamma[,  iind,j]*(sigma+(mu*t(mu)))
      }
      new_a[j,iind]=(1/(a_de))*a_nu/2
      #rescale
      new_a[j,]=new_a[j,]%*%d_temp1
    }
    
    #L1-penalty on gamma
    par_gamma = new_gamma
    for (yy in 1:(y-1)){
      for (j in 1:item){
        if (sum(par_gamma[yy,,j])!=0){
          rr=which(par_gamma[yy,,j]!=0)
          delta = 0; deriv2 = 0
          for (i in 1:person){
            sigma=SIGMA[rr,rr,i]
            mu=MU[rr,i]
            sigmumu = sigma + mu%*%t(mu)
            delta = delta + (((as.numeric(resp[i,j])*G[i,yy]-0.5*G[i,yy])%*%mu)+(2*(new_b[j]-new_beta[j,]%*%G[i,])*eta[i,j]*mu*G[i,yy])
                             -(2*eta[i,j]*G[i,yy]*new_a[j,rr]%*%sigmumu))
            deriv2 = deriv2 + 2*eta[i,j]*sigmumu*G[i,yy]^2
          }
          if (abs(delta) > lbd  ){
            S = sign(delta)*(abs(delta) - lbd)     
            new_gamma[yy,rr,j] = qr.solve(deriv2,t(S))
            new_gamma[yy,,j]=new_gamma[yy,,j]%*%d_temp1
          } else {
            new_gamma[yy,rr,j] = 0
          }
        } else {
          new_gamma[yy,,j] = 0
        }
      }
    }
    
    #L1-penalty on beta
    par_beta = new_beta
    for (yy in 1:(y-1)){
      for (j in 1:item){
        if (par_beta[j,yy]!=0){
          delta = 0; deriv2 = 0
          for (i in 1:person){
            delta = delta + resp[i,j]*G[i,yy]-0.5*G[i,yy]-eta[i,j]*(2*G[i,yy]*(new_a[j,]+G[i,yy]%*%new_gamma[yy,,j])%*%MU[,i]-2*new_b[j]*G[i,yy])
            deriv2 = deriv2 + 2*eta[i,j]*(G[i,yy]^2)
          }
          if (abs(delta) > lbd  ){
            S = sign(delta)*(abs(delta) - lbd)     
            new_beta[j,yy] = S/deriv2
          } else{
            new_beta[j,yy] = 0
          }
        } else {
          new_beta[j,yy] = 0
        }
      }
    }
    sparsity1=new_gamma
    for (j in 1:item){
      for (rr in 1:domain){
        for (nn in 1:(y-1)){
          sparsity1[nn,rr,j]=ifelse(new_gamma[nn,rr,j]==0,0,1)
        }
      }
    }
    sparsity2=new_beta
    for (j in 1:item){
      for (rr in 1:(y-1)){
        sparsity2[j,rr]=ifelse(new_beta[j,rr]==0,0,1)
      }
    } 
    
    #Second M-step for on gamma
    #par_gamma = new_gamma
    #L1-penalty on gamma
    for (yy in 1:(y-1)){
      for (j in 1:item){
        if (sum(par_gamma[yy,,j])!=0){
          rr=which(par_gamma[yy,,j]!=0)
          delta = 0; deriv2 = 0
          for (i in 1:person){
            sigma=SIGMA[rr,rr,i]
            mu=MU[rr,i]
            sigmumu = sigma + mu%*%t(mu)
            delta = delta + (((as.numeric(resp[i,j])*G[i,yy]-0.5*G[i,yy])%*%mu)+(2*(new_b[j]-new_beta[j,]%*%G[i,])*eta[i,j]*mu*G[i,yy])
                             -(2*eta[i,j]*G[i,yy]*new_a[j,rr]%*%sigmumu))
            deriv2 = deriv2 + 2*eta[i,j]*sigmumu*G[i,yy]^2
          }
          if (abs(delta) > 0  ){
            S = sign(delta)*(abs(delta) )     
            new_gamma[yy,rr,j] = qr.solve(deriv2,t(S))
            new_gamma[yy,,j]=new_gamma[yy,,j]%*%d_temp1
          } else {
            new_gamma[yy,rr,j] = 0
          }
        } else {
          new_gamma[yy,,j] = 0
        }
      }
    }
    new_gamma=new_gamma*sparsity1
    
    
    
    #Second M-step for beta
    #par_beta = new_beta
    for (yy in 1:(y-1)){
      for (j in 1:item){
        if (par_beta[j,yy]!=0){
          delta = 0; deriv2 = 0
          for (i in 1:person){
            delta = delta + resp[i,j]*G[i,yy]-0.5*G[i,yy]-eta[i,j]*(2*G[i,yy]*(new_a[j,]+G[i,yy]%*%new_gamma[yy,,j])%*%MU[,i]-2*new_b[j]*G[i,yy])
            deriv2 = deriv2 + 2*eta[i,j]*(G[i,yy]^2)
          }
          if (abs(delta) > 0  ){
            S = sign(delta)*(abs(delta))     
            new_beta[j,yy] = S/deriv2
          } else{
            new_beta[j,yy] = 0
          }
        } else {
          new_beta[j,yy] = 0
        }
      }
    }
    new_beta=new_beta*sparsity2
    
    if (max(abs(new_a-par_a))< 0.001 & max(abs(new_b-par_b))< 0.001 & max(abs(new_gamma-par_gamma)) < 0.001 & max(abs(new_beta-par_beta)) < 0.001 &  max(abs(Mu-par_Mu)) < 0.001 &  max(abs(Sigma-par_Sigma)) < 0.001){
      converge=1
    } 
    n = n+1
    #print(n)
    #print(new_beta)
    #print(new_gamma)
    #rankSigma=rankMatrix(Sigma[,,1])[1]
    
  
  }
  
  
  #second EM cycles
  
  converge=0
  singularity = 0
  n=0
  # EM cycle
  while(converge == 0  && rankSigma == domain)
  {
    #Step 1: Update mu and sigma & Step 2.1: Update xi and eta(xi). Meanwhile, updata Sigma. -2*eta[i,j]*(new_beta[,j]%*%G[i,])
    par_Sigma=Sigma 
    par_Mu=Mu 
    yy=1
    uy=resp[which(Gind==yy),]
    Muy=Mu[,yy]
    Sigmay=as.matrix(Sigma[,,yy])
    persony=N.vec[yy]
    etay=eta[which(Gind==yy),]
    new_ay=new_a #new_a[,j]+G[i,]%*%new_gamma[,,j]
    new_by=as.vector(new_b)#new_b[j]-new_beta[,j]
    rs1<-ecfa2(as.matrix(uy),Muy,Sigmay,domain, persony, item, etay, new_ay, new_by)
    SIGMA[,,which(Gind==yy)]=rs1$SIGMA
    MU[,which(Gind==yy)]=rs1$MU
    xi[which(Gind==yy),]=rs1$xi
    Sigma[,,yy]=rs1$Spart/N.vec[1]
    for (yy in 2:y){
      Ny=N.vec[yy]
      uy=resp[which(Gind==yy),]
      Muy=Mu[,yy]
      Sigmay=Sigma[,,yy]
      persony=N.vec[yy]
      etay=eta[which(Gind==yy),]
      new_ay=new_a+  t(new_gamma[yy-1,,]) #new_a[,j]+G[i,]%*%new_gamma[,,j]
      new_by=as.vector(new_b-new_beta[,yy-1])#new_b[j]-new_beta[,j]
      rs1<-ecfa2(as.matrix(uy),Muy,Sigmay,domain, persony, item, etay, new_ay, new_by)
      SIGMA[,,which(Gind==yy)]=rs1$SIGMA
      MU[,which(Gind==yy)]=rs1$MU
      xi[which(Gind==yy),]=rs1$xi
      Sigma[,,yy]=rs1$Spart/Ny 
      Mu[,yy]=rowMeans(rs1$MU)
    }
    d_temp1=sqrt(diag(diag(Sigma[,,1])))
    for (yy in 1:y){
      Sigma[,,yy]=solve(d_temp1)%*%Sigma[,,yy]%*%solve(d_temp1)
    }
    
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(eta)< 0.01,0.125)
    
    # Step 2.2: Update b.
    par_b=new_b
    b_part = matrix(0,person,item)
    for (i in 1:person){
      new_gamma_i=NULL
      for (j in 1:item){
        new_gamma_ij=G[i,]%*%new_gamma[,,j]
        new_gamma_i=cbind(new_gamma_i,t(new_gamma_ij))
      }
      b_part[i,]=(new_a+t(new_gamma_i))%*%MU[,i]+new_beta%*%G[i,] 
    }
    new_b=colSums(0.5-resp+2*eta*b_part)/colSums(2*eta)
    
    # Step 2.3: Update a.
    par_a=new_a
    for (j in 1:item){
      a_nu=0;a_de=0
      Ind=indic[,j]
      iind=which(Ind==1)
      for (i in 1:person){
        sigma=SIGMA[,,i]
        sigma=sigma[iind,iind]
        mu=MU[iind,i]
        a_de=a_de+eta[i,j]*sigma+eta[i,j]*(mu*t(mu))
        a_nu=a_nu+(resp[i,j]-0.5+2*(new_b[j]-new_beta[j,]%*%G[i,])*eta[i,j])*mu 
        -eta[i,j]*2*G[i,]%*%new_gamma[,  iind,j]*(sigma+(mu*t(mu)))
      }
      new_a[j,iind]=(1/(a_de))*a_nu/2
      #rescale
      new_a[j,]=new_a[j,]%*%d_temp1
    }
    
    #L1-penalty on gamma
    par_gamma = new_gamma
    for (yy in 1:(y-1)){
      for (j in 1:item){
        if (sum(par_gamma[yy,,j])!=0){
          rr=which(par_gamma[yy,,j]!=0)
          delta = 0; deriv2 = 0
          for (i in 1:person){
            sigma=SIGMA[rr,rr,i]
            mu=MU[rr,i]
            sigmumu = sigma + mu%*%t(mu)
            delta = delta + (((as.numeric(resp[i,j])*G[i,yy]-0.5*G[i,yy])%*%mu)+(2*(new_b[j]-new_beta[j,]%*%G[i,])*eta[i,j]*mu*G[i,yy])
                             -(2*eta[i,j]*G[i,yy]*new_a[j,rr]%*%sigmumu))
            deriv2 = deriv2 + 2*eta[i,j]*sigmumu*G[i,yy]^2
          }
          if (abs(delta) > 0  ){
            S = sign(delta)*(abs(delta) - 0)     
            new_gamma[yy,rr,j] = qr.solve(deriv2,t(S))
            new_gamma[yy,,j]=new_gamma[yy,,j]%*%d_temp1
          } else {
            new_gamma[yy,rr,j] = 0
          }
        } else {
          new_gamma[yy,,j] = 0
        }
      }
    }
    
    #L1-penalty on beta
    par_beta = new_beta
    for (yy in 1:(y-1)){
      for (j in 1:item){
        if (par_beta[j,yy]!=0){
          delta = 0; deriv2 = 0
          for (i in 1:person){
            delta = delta + resp[i,j]*G[i,yy]-0.5*G[i,yy]-eta[i,j]*(2*G[i,yy]*(new_a[j,]+G[i,yy]%*%new_gamma[yy,,j])%*%MU[,i]-2*new_b[j]*G[i,yy])
            deriv2 = deriv2 + 2*eta[i,j]*(G[i,yy]^2)
          }
          if (abs(delta) > 0  ){
            S = sign(delta)*(abs(delta) - 0)     
            new_beta[j,yy] = S/deriv2
          } else{
            new_beta[j,yy] = 0
          }
        } else {
          new_beta[j,yy] = 0
        }
      }
    }
    sparsity1=new_gamma
    for (j in 1:item){
      for (rr in 1:domain){
        for (nn in 1:(y-1)){
          sparsity1[nn,rr,j]=ifelse(new_gamma[nn,rr,j]==0,0,1)
        }
      }
    }
    sparsity2=new_beta
    for (j in 1:item){
      for (rr in 1:(y-1)){
        sparsity2[j,rr]=ifelse(new_beta[j,rr]==0,0,1)
      }
    } 
    
    
    if (max(abs(new_a-par_a))< 0.001 & max(abs(new_b-par_b))< 0.001 & max(abs(new_gamma-par_gamma)) < 0.001 & max(abs(new_beta-par_beta)) < 0.001 &  max(abs(Mu-par_Mu)) < 0.001 &  max(abs(Sigma-par_Sigma)) < 0.001){
      converge=1
    } 
    n = n+1
  }
  
  Sigma.est=matrix(0,domain*y,domain)
  Mu.est = numeric(domain*y)
  for (yy in 1:y){
    Sigma.est[((yy-1)*domain+1):(yy*domain),]=Sigma[,,yy]
    Mu.est[((yy-1)*domain+1):(yy*domain)]=Mu[,yy]
  }
  
  
  
  
  
  # BIC
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
  Sig.est.slice=array(0,c(r,r,y))
  for (yy in 1:y){
    Sig.est.slice[,,yy]=Sigma.est[((yy-1)*r+1):((yy-1)*r+r),]
  }
  X1=seq(-3,3,by=0.2)
  G2=length(X1)^r
  gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
  idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
  X2 <- matrix(gh[idx,1],nrow(idx),r)
  Xijk=array(double(N*J*m),dim = c(N,J,m))
  for(i in 1:N){
    for(j in 1:J){
      for(k in 1:m){
        Xijk[i,j,k]=ifelse(resp[i,j]==k,1,0)
      }
    }
  }
  gra=new_a
  grd=as.matrix(-new_b)
  grbeta=new_beta
  grgamma=new_gamma
  y.allgroup=rbind(rep(0,y-1),diag(y-1)) #y1, y2, y3
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X2,y=y,G=G2,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=new_a, grd=as.matrix(-new_b), grbeta=new_beta, grgamma=new_gamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G2)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G2,N=N,m=m)
    sumoverk0=sumoverk(G=G2,rgky=rgk[(G2+1):(2*G2),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=numeric(r),X=X)
    for (yy in 2:y){
      sumoverk0=sumoverk0+sumoverk(G=G2,rgky=rgk[(yy*G2+1):((yy+1)*G2),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,(yy-1)],gamjy=grgamma[(yy-1),,j],X=X2)
    }
    temp=sumoverk0#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  
  
  # Likelihood 
  A.allgroups=numeric(nrow(X)*N) #A1, A2, A3
  Xij.all=matrix(0,N,(y-1))
  Xij.all[(N1+1):(N1+N2),1]=1
  Xij.all[(N1+N2+1):(N1+N2+N3),2]=1
  for (i in 1:N1){
    Mu2=Mu.est.mat[1,]
    Sig=Sig.est.slice[,,1]
    A.allgroups[((i-1)*nrow(X)+1):((i-1)*nrow(X)+nrow(X))]=dmvnorm(X, Mu2,Sig)
  }
  for (i in (N1+1):(N1+N2)){
    Mu2=Mu.est.mat[2,]
    Sig=Sig.est.slice[,,2]
    A.allgroups[((i-1)*nrow(X)+1):((i-1)*nrow(X)+nrow(X))]=dmvnorm(X, Mu2,Sig)
  }
  for (i in (N1+N2+1):(N1+N2+N3)){
    Mu2=Mu.est.mat[3,]
    Sig=Sig.est.slice[,,3]
    A.allgroups[((i-1)*nrow(X)+1):((i-1)*nrow(X)+nrow(X))]=dmvnorm(X, Mu2,Sig)
  }
  #loglikelihood(arma::cube Xijk, const arma::mat&new_a,const arma::mat&new_d,arma::cube grbeta, arma::vec Aallgroups,arma::vec ng,arma::mat X,arma::mat Xijall,int r,int N,int G,int J, int m, int y)
  Xij.all=matrix(0,N,y-1)
  Xij.all[(N/3+1):(N/3*2),1]=1
  Xij.all[(N/3*2+1):(N),2]=1
  grbeta2=array(0,c(20,2,1))
  grbeta2[,,1]=grbeta
  lh2=loglikelihood(as.matrix(resp2),new_a,as.matrix(new_b),new_beta,A.allgroups,X2,G,domain,person,G2,item)
    
  
  # GIC BIC
  ll=0
  for (j in 1:item){
    ll1=0
    for (i in 1:person){
      ll1=ll1+(log(exp(xi[i,j])/(1+exp(xi[i,j]))) + (0.5-resp[i,j])*new_b[j] - (0.5-resp[i,j])*(new_beta[j,]%*%G[i,]) 
               + (resp[i,j]-0.5)*((new_a[j,]+G[i,]%*%new_gamma[,,j])%*%MU[,i]) - 0.5*xi[i,j]
               - eta[i,j]*(new_b[j]^2 + (new_beta[j,]%*%G[i,])^2 - 2*new_b[j]*((new_a[j,]+G[i,]%*%new_gamma[,,j])%*%MU[,i]) 
                           + 2*(new_beta[j,]%*%G[i,])*((new_a[j,]+G[i,]%*%new_gamma[,,j])%*%MU[,i]) 
                           + (new_a[j,]+G[i,]%*%new_gamma[,,j])%*%(SIGMA[,,i]+MU[,i]%*%t(MU[,i]))%*%t(new_a[j,]+G[i,]%*%new_gamma[,,j]) 
                           - xi[i,j]^2))
    }
    ll=ll+ll1
  }
  ll2=0
  for (yy in 1:y){
    ll2=ll2+N.vec[yy]/2*log(det(solve(Sigma[,,yy])))
  }
  ll3=0
  for (yy in 1:y){
    Ny=N.vec[yy]
    Ny0=sum(N.vec[0:(yy-1)])
    for (i in 1:Ny){
      ll3=ll3+0.5*tr(solve(Sigma[,,yy])%*%(SIGMA[,,(Ny0+i)]+(MU[,(Ny0+i)]-Mu[,yy])%*%t(MU[,(Ny0+i)]-Mu[,yy]))) 
    }
  }
  ll=ll+ll2-ll3
  
  l0norm1=l0norm2=0
  for(rr in 1:domain){
    for(j in 1:item){
      for(yy in 1:(y-1)){
        l0norm1=l0norm1+(new_gamma[yy,rr,j]!=0)
      }
    }
  }
  for(i in 1:item){
    for(j in 1:(y-1)){
      l0norm2=l0norm2+(new_beta[i,j]!=0)
    }
  }
  l0norm=l0norm1+l0norm2
  GIC=-2*lh2+cons*log(log(person))*log(person)*l0norm
  BIC=-2*lh2+log(person)*l0norm
  #GIC=-2*sum(lh)+cons*log(log(person))*log(person)*l0norm
  #BIC=-2*sum(lh)+log(person)*l0norm
  #BIC=-2*ll+log(person)*l0norm
  
  #Mu.gp1=Mu.est[1:2];Mu.gp2=Mu.est[3:4];Mu.gp3=Mu.est[5:6]
  #Sig.gp1=Sig.est[1:2,];Sig.gp2=Sig.est[3:4,];Sig.gp3=Sig.est[5:6,]
  #return(list(est=cbind(new_a,new_b),Gamma=new_gamma,Beta=new_beta,n=n,bic=BIC,gic=GIC, means=Mu.est,Covs=Sigma.est))
  return(list(est=cbind(new_a,new_b),Gamma=new_gamma,Beta=new_beta,n=n,bic=BIC,gic=GIC,Likelihood=ll,Likelihood2=sum(lh),Likelihood3=lh2, means=Mu.est,Covs=Sigma.est))
}

l0norm=numeric(length(lbd.vec))
for (k in 1:length(lbd.vec)){
  l0norm[k]=sum(Betas[,,k]!=0)
}
which.min(-2*lls)
which.min(-2*lls2)
which.min(-2*lls3)
which.min(-2*lls3+log(person)*l0norm)
which.min(-2*lls3+0.1*log(log(person))*log(person)*l0norm)
cons.seq=seq(0.1,0.3,0.001)
for (cons in cons.seq){
  print(which.min((-2*lls3+cons*log(log(person))*log(person)*l0norm)))
}

which.min(-2*lls2+2*l0norm)


which.min((-2*lls2+cons*log(log(person))*log(person)*l0norm)[-1])
which.min(-2*lls3+cons*log(log(person))*log(person)*l0norm)
which.min(-2*lls+log(person)*l0norm)
which.min((-2*lls2+log(person)*l0norm)[-1])
which.min(-2*lls3+log(person)*l0norm)
lbd.vec[which.min(-2*lls+log(person)*l0norm)]
Betas[,,which.min(-2*lls+log(person)*l0norm)]
Betas[,,which.min(-2*lls2+log(person)*l0norm)]
lbd.vec[which.min(-2*lls+cons*log(log(person))*log(person)*l0norm)]
Betas[,,which.min(-2*lls2+cons*log(log(person))*log(person)*l0norm)]
cons=0.5

