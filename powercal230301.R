# Table 2 Lasso GVEM
power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM1_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/3
sd(power1[,2])/3
sd(power1[,3])/3
sd(TpI1[,1])/3
sd(TpI1[,2])/3
sd(TpI1[,3])/3

# Table 2 Lasso GVEMM
power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:10){
  lbdall[rep]=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM1_EMM_lbd_",rep))
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM1_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim1_Condition1_Beta1_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
}
colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)

# Table 2 Lasso GVEMM2
power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:10){
  lbdall[rep]=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM1_EMM_lbd_",rep))
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM1_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim1_Condition1_Beta2_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)


# Table 2 Lasso GVEM
power1=TpI1=matrix(0,10,3)
lbdall=numeric(rep)
for (rep in 1:10){
  lbdall[rep]=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM2_lbd_",rep))
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM2_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/3
sd(power1[,2])/3
sd(power1[,3])/3
sd(TpI1[,1])/3
sd(TpI1[,2])/3
sd(TpI1[,3])/3

# Table 2 Lasso GVEMM

power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:10){
  lbdall[rep]=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM2_EMM_lbd_",rep))
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM2_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim1_Condition2_Beta1_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)

# Table 2 Lasso GVEMM2
power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:10){
  lbdall[rep]=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM2_EMM_lbd_",rep))
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM2_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}
for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim1_Condition2_Beta2_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
}
colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)



# Table 2 Lasso GVEM
power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM3_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}


colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/3
sd(power1[,2])/3
sd(power1[,3])/3
sd(TpI1[,1])/3
sd(TpI1[,2])/3
sd(TpI1[,3])/3

# Table 2 Lasso GVEMM
power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:10){
  lbdall[rep]=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM3_EMM_lbd_",rep))
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM3_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}
for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim1_Condition3_Beta1_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)

# Table 2 Lasso GVEMM2
power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:10){
  lbdall[rep]=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM3_EMM_lbd_",rep))
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM3_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}
for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim1_Condition3_Beta2_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)



# Table 2 Lasso GVEM
power1=TpI1=matrix(0,5,3)
for (rep in 1:5){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM4_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/3
sd(power1[,2])/3
sd(power1[,3])/3
sd(TpI1[,1])/3
sd(TpI1[,2])/3
sd(TpI1[,3])/3

# Table 2 Lasso GVEMM

power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM4_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}
for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim1_Condition4_Beta1_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)

# Table 2 Lasso GVEMM2
power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM4_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim1_Condition4_Beta2_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)








# Table 4 Lasso GVEM
power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM5_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}



# Table 2 Lasso GVEMM
power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM5_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}

for (rep in 1:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEM_DIF/GVEM_sim_results/GVEM_Sim2_Condition6_Beta1_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
}


colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)

# Table 2 Lasso GVEMM2
power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM5_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim2_Condition5_Beta2_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
}


colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)


# Table 2 Lasso GVEM
power1=TpI1=matrix(0,25,3)
for (rep in 1:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEM_DIF/GVEM_sim_results/GVEM_Sim2_Condition8_Beta1_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/3
sd(power1[,2])/3
sd(power1[,3])/3
sd(TpI1[,1])/3
sd(TpI1[,2])/3
sd(TpI1[,3])/3

# Table 2 Lasso GVEMM

power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM6_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim2_Condition6_Beta1_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
}


colMeans(power1[11:25,])
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)

# Table 2 Lasso GVEMM2
power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM6_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim2_Condition6_Beta2_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
}


colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)




# Table 2 Lasso GVEM
power1=TpI1=matrix(0,10,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM7_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/3
sd(power1[,2])/3
sd(power1[,3])/3
sd(TpI1[,1])/3
sd(TpI1[,2])/3
sd(TpI1[,3])/3

# Table 2 Lasso GVEMM
power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM7_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim2_Condition7_Beta1_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
}


colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)


# Table 2 Lasso GVEMM2
power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM7_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim2_Condition7_Beta2_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,12,13),])!=0)!=0)/4
  power1[rep,2]=sum((beta[c(4,5,12,13),1])!=0)/4
  power1[rep,3]=sum((beta[c(4,5,12,13),2])!=0)/4
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,12,13),])!=0)!=0)/16
  TpI1[rep,2]=sum((beta[-c(4,5,12,13),1])!=0)/16
  TpI1[rep,3]=sum((beta[-c(4,5,12,13),2])!=0)/16
}


colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)



# Table 2 Lasso GVEM
power1=TpI1=matrix(0,5,3)
for (rep in 1:5){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM8_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/3
sd(power1[,2])/3
sd(power1[,3])/3
sd(TpI1[,1])/3
sd(TpI1[,2])/3
sd(TpI1[,3])/3

# Table 2 Lasso GVEMM

power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM8_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim2_Condition8_Beta1_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
}


colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)


# Table 2 Lasso GVEMM2
power1=TpI1=matrix(0,25,3)
for (rep in 1:10){
  beta=read.csv(paste("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/newGIC_GVEM8_EMM_Beta_",rep))
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
  
}

for (rep in 11:25){
  beta=read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEMRegDIF_SimData/GVEM DIF additional 15 rep/GVEM_Sim2_Condition8_Beta2_",rep))[,-1]
  power1[rep,1]=sum(rowSums((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/12
  power1[rep,2]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/12
  power1[rep,3]=sum((beta[c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/12
  TpI1[rep,1]=sum(rowSums((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),])!=0)!=0)/8
  TpI1[rep,2]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),1])!=0)/8
  TpI1[rep,3]=sum((beta[-c(4,5,6,7,8,9,12,13,14,15,16,17),2])!=0)/8
}


colMeans(power1)
colMeans(TpI1)
sd(power1[,1])/sqrt(24)
sd(power1[,2])/sqrt(24)
sd(power1[,3])/sqrt(24)
sd(TpI1[,1])/sqrt(24)
sd(TpI1[,2])/sqrt(24)
sd(TpI1[,3])/sqrt(24)


