PRS = function(n,g_causal,g,b){
  sim=100
  corr = rep(0,sim)
  for(s in 1:sim){
    corr[s]=PRS_sim(n,g_causal,g,b)
    print(paste("Simulation ",s,"correlation: ",corr[s]))
  }
  return(mean(corr))
}

PRS_sim <- function(n,g_causal,g, b){
  #n= training sample size, g_causual= number of casual SNPs, g= total number of SNPS, 
  #create training set X (columns: Y, Causual SNP genotypes, non-causual SNP genotypes)
  dat = data.frame(matrix(0,nrow=n,ncol=g+1))
  
  #number of non-causal SNPs
  g_non = g - g_causal
  
  #simulate MAF for casual SNPs and non-casual SNPs
  p = c(runif(g_causal,0.1,0.2),runif(g_non,0.4,0.5))
  
  ##simulate genotype dataset for SNPs in training set X
  for(i in 1:n){
    dat[i,2:(g+1)] = rbinom(g,2,p)
  }
  
  ##Calculate Y for training set  
  beta=matrix(b) 
  dat[,1]=as.matrix(dat[,2:(g_causal+1)])%*%beta 
  
  ##Fit linear regression to find weights of PRS score
  beta.pred=matrix(NA,nrow=g, ncol=1) #create vector of predicted betas
  
  for(i in 2:(g+1)){
    beta.pred[i-1,] = coefficients(glm(dat[,1] ~ dat[,i], family=gaussian))[2]
  }
  
  #####Simulate test set Z ####
  datz = matrix(0,nrow=ceiling(n/2),ncol=g+1)
  for(i in 1:ceiling(n/2)){
    datz[i,2:(g+1)] = rbinom(g,2,p)
  }
  
  ##calculate PRS of the Z 
  riskScore = datz[,2:(g+1) ]%*%beta.pred
  
  #calculate true Y of Z 
  datz[,1]=as.matrix(datz[,2:(g_causal+1)])%*%beta 
  
  
  ######Correlation between true Y of test set and PRS score ####
  Corr=cor(datz[,1], riskScore)
  return(Corr)
}