#################################################################
#################################################################
## CWVS-E1 Data Example: 
## Time-varying effects of third-trimester GWG on 37-week EFW
#################################################################
#################################################################

library(mvtnorm)
library(msm)
library(stats)
library(stringr)
library(missMDA)

########################################
# Install R package from Github
library(devtools)
devtools::install_github('mstville/CWe1')
library(CWe1)

gest.ages=27:37
m=length(gest.ages)

data=read.csv("EFW37_Weight_Dataset.csv",header=T)  ## Set working directory appropriately to access data

head(data)
summary(data)
data$Education=as.factor(data$Education);summary(data$Education)
data$Race=as.factor(data$Race);summary(data$Race)
data$InfSex=as.factor(data$InfSex);summary(data$InfSex)
data$Parity=as.factor(data$Parity);summary(data$Parity)
data$GDM=as.factor(data$GDM);summary(data$GDM)
data$HTN=as.factor(data$HTN);summary(data$HTN)


pats=unique(data$ID)
n=length(pats);n

################################################
## Create GWG Exposure matrix
w.id=which(str_detect(names(data),"Weight")==TRUE)
W.dat=data[,w.id]

Z=matrix(-99,nrow=n,ncol=m)
for (t in 1:m){
  id=which(names(W.dat)==paste0("Weight",gest.ages[t]))
  Wt=W.dat[,c(id-1,id)]
  Zt=apply(Wt,1,diff)
  Z[,t]=Zt
}
colnames(Z)=paste0("GWG",gest.ages)
head(Z)
data=cbind.data.frame(data,as.data.frame(Z));head(data)

################################################
## Create covariate matrix
X=model.matrix(~Race+Age+Education+InfSex+Cotinine+GDM+HTN+Parity+PreBMI,data=data)
colnames(X)[c(1,15)]=c("Intercept","Parity2")
p=ncol(X);p

################################################
## Outcome: 37-week EFW
Y=data$EFW  


#-------------------------------------------------------------------------------------
###############################################################################################
## MODEL ESTIMATION
set.seed(221)

G=5000 ## MCMC samples

# initial values
a_e=b_e=0.01
sig2_e=1

sig2_b=10000
beta=rep(0,p)

gamma=rep(1,m)

delta1=rep(0,m)
delta2=rep(0,m)

a_phi=b_phi=1
phi1=phi2=-log(0.05)/(m-1)
temp_corr_info1=temporal_corr_fun(m,phi1)
temp_corr_info2=temporal_corr_fun(m,phi2)

sig2_a=1
A11=A22=1
A21=0

alpha=rep(-99,m)
for (t in 1:m){
  alpha[t]=gamma[t]*(A11*delta1[t])
}

eta=A21*delta1+A22*delta2



sig2.e.mcmc=rep(-99,G)
beta.mcmc=matrix(-99,nrow=G,ncol=p)
gamma.mcmc=matrix(-99,nrow=G,ncol=m)
delta1.mcmc=matrix(-99,nrow=G,ncol=m)
delta2.mcmc=matrix(-99,nrow=G,ncol=m)
phi1.mcmc=rep(-99,G)
phi2.mcmc=rep(-99,G)
A11.mcmc=rep(-99,G)
A22.mcmc=rep(-99,G)
A21.mcmc=rep(-99,G)
theta.mcmc=matrix(-99,nrow=G,ncol=m)
alpha.mcmc=matrix(-99,nrow=G,ncol=m)
eta.mcmc=matrix(-99,nrow=G,ncol=m)
lag.mcmc=rep(-99,G)

acctot_phi1_trans = acctot_phi2_trans = acctot_A11_trans = acctot_A22_trans = 0

#####################################
### GIBBS SAMPLER STARTS HERE
for (g in 1:G){
  ###########################
  # Update residual variance
  sig2_e=sigma2_epsilon_update(a_sigma2_epsilon = a_e, b_sigma2_epsilon = b_e,
                               y=Y, x=X, z=Z, beta_old = beta, gamma_old = gamma,
                               A11_old = A11, delta1_old = delta1)
  sig2.e.mcmc[g]=sig2_e
  
  w=rep((1/sig2_e),n)
  
  ###########################
  # Update beta
  beta=beta_update(x=X, z=Z, off_set=rep(0,n), sigma2_beta=sig2_b, w=w, gamma_l=Y, gamma_old=gamma, A11_old=A11, delta1_old=delta1)
  beta=as.vector(beta)
  beta.mcmc[g,]=beta
  
  ###########################
  # Update gamma
  gamma=gamma_update(x=X, z=Z, off_set=rep(0,n), w=w,gamma_l=Y,beta=beta,gamma_old=gamma,A11_old=A11,delta1_old=delta1,A21_old=A21,A22_old=A22,delta2_old = delta2)
  gamma=as.vector(gamma)
  # gamma=gam_update_arma(X,Z,w,Y,beta,gamma,A11,delta1,A21,A22,delta2,E1=TRUE)
  gamma.mcmc[g,]=gamma
  if (length(which(gamma==1))>0){
    tstar=max(which(gamma==1))
  } else {tstar=0}
  lag.mcmc[g]=m-tstar
  
  ###########################
  # Sample latent variables
  # from truncated normal
  gamma_star=gamma_star_update(gamma=gamma, delta1_old=delta1, A21_old=A21, 
                               A22_old=A22, delta2_old=delta2)
  gamma_star=as.vector(gamma_star)
  
  ###########################
  # Update delta1
  delta1=delta1_update(x=X, z=Z, off_set=rep(0,n), w=w, gamma_l=Y, beta=beta, 
                       gamma=gamma, gamma_star=gamma_star, A11_old=A11, A21_old=A21, A22_old=A22,
                       delta2_old=delta2, corr_inv1=temp_corr_info1[[1]])
  delta1=as.vector(delta1)
  delta1.mcmc[g,]=delta1
  
  ###########################
  # Update delta2
  delta2=delta2_update(gamma_star, delta1, A21_old=A21, A22_old=A22, corr_inv2=temp_corr_info2[[1]])
  delta2=as.vector(delta2)
  delta2.mcmc[g,]=delta2
  
  ###########################
  # Update phi1
  phi1_up=phi_update(phi_old=phi1, delta=delta1, temporal_corr_info=temp_corr_info1, 
                     alpha_phi=a_phi, beta_phi=b_phi, 
                     metrop_var_phi_trans=1.00, acctot_phi_trans=acctot_phi1_trans)
  temp_corr_info1=phi1_up$temporal_corr_info
  acctot_phi1_trans=phi1_up$acctot_phi_trans
  phi1=phi1_up$phi
  phi1.mcmc[g]=phi1
  
  ###########################
  # Update phi2
  phi2_up=phi_update(phi_old=phi2, delta=delta2, temporal_corr_info=temp_corr_info2, 
                     alpha_phi=a_phi, beta_phi=b_phi, 
                     metrop_var_phi_trans=1.00, acctot_phi_trans=acctot_phi2_trans)
  temp_corr_info2=phi2_up$temporal_corr_info
  acctot_phi2_trans=phi2_up$acctot_phi_trans
  phi2=phi2_up$phi
  phi2.mcmc[g]=phi2
  
  ###########################
  # Update A11
  A11_up=A11_update(A11_old=A11, x=X, z=Z, off_set=rep(0,n), w=w, gamma_l=Y, 
                    beta=beta, gamma=gamma, delta1=delta1, sigma2_A=sig2_a,
                    metrop_var_A11_trans=0.03, acctot_A11_trans=acctot_A11_trans)
  A11=A11_up$A11
  A11.mcmc[g]=A11
  acctot_A11_trans=A11_up$acctot_A11_trans
  
  ###########################
  # Update A22
  A22_up=A22_update(A22_old=A22, gamma_star=gamma_star, delta1=delta1, 
                    A21_old=A21, delta2=delta2, sigma2_A=sig2_a, 
                    metrop_var_A22_trans=0.30, acctot_A22_trans=acctot_A22_trans)
  A22=A22_up$A22
  A22.mcmc[g]=A22
  acctot_A22_trans=A22_up$acctot_A22_trans
  
  ###########################
  # Update A21
  A21=A21_update(gamma_star=gamma_star, delta1=delta1, A22=A22, delta2=delta2, sigma2_A=sig2_a)
  A21.mcmc[g]=A21
  
  ###########################
  # Recover theta
  theta=unlist(lapply(1:m,function(t){A11*delta1[t]}))
  theta.mcmc[g,]=theta
  
  ###########################
  # Recover alpha
  for (t in 1:m){alpha[t]=gamma[t]*(A11*delta1[t])}
  alpha.mcmc[g,]=alpha
  
  ###########################
  # Recover eta
  eta=A21*delta1+A22*delta2
  eta.mcmc[g,]=eta
  
  print(g)
}
### GIBBS SAMPLER ENDS HERE
#####################################

burn=2500:G

gamma=apply(gamma.mcmc[burn,],2,mean);gamma

alpha=apply(alpha.mcmc[burn,],2,mean);alpha
alpha.ci=apply(alpha.mcmc[burn,],2,quantile,probs=c(.025,.975));alpha.ci
