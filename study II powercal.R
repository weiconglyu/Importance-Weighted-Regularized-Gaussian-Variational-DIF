power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:25){
  beta=(read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEM_DIF/GVEM_sim_results/GVEM_Sim2_Condition5_Beta2_",rep)))[,-1]
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


power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:25){
  beta=(read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEM_DIF/GVEM_sim_results/GVEM_Sim2_Condition6_Beta2_",rep)))[,-1]
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



power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:25){
  beta=(read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEM_DIF/GVEM_sim_results/GVEM_Sim2_Condition7_Beta2_",rep)))[,-1]
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


power1=TpI1=matrix(0,25,3)
lbdall=numeric(rep)
for (rep in 1:25){
  beta=(read.csv(paste0("/Users/ruoyizhu/Documents/Github/GVEM_DIF/GVEM_sim_results/GVEM_Sim2_Condition8_Beta2_",rep)))[,-1]
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
