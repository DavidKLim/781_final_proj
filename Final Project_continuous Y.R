#Bios 781 Final Project 
##### Outline of task 
#1. simulate g_causual x_i and simulate g-g_causal x_i (using the different MAF as before) for n individuals
#2. Calcuate y by setting g_casual betas 
#3. Fit g_causal glm models (Y=beta1+error... Y=beta_gcausal +error)
#4. Simulate another set of n individuals ( same way as #1) -- test set
#5. Calculate polygenic risk score for test set (use beta estimates of #3)
#6. Calculate actual disease status using set betas (set in #2)
#7. If this method is good, then the PRS and Disease status of test set will be correlated 
#I think we can just use sensitivity and FDR. I will think if there is another correlated stat that is better 
#change g_casual, g, and n 

########################Step 1: Creating Function###################################
set.seed(9)
# inputs: n, g_causal, g, prevalence

PRS <- function(n,g_causal,g){
  #n= training sample size, g_causual= number of casual SNPs, g= total number of SNPS, 
  #create training set X (columns: Y, Causual SNP genotypes, non-causual SNP genotypes)
  dat = data.frame(matrix(0,nrow=n,ncol=g+1))
  
  #number of non-causal SNPs
  g_non = g - g_causal
  
  #simulate MAF for casual SNPs and non-casual SNPs
  p = c(runif(g_causal,0.1,0.2),runif(g_non,0.4,0.5))
  
  ##simulate genotype dataset for SNPs in test set X
  for(i in 1:n){
    for(j in 1:g){
      dat[i,j+1] = rbinom(1,2,p[j])
    }
  }
  
  ##Calculate Y for test set  
  beta=matrix(c(.2, .4)) #could put beta into parameter list (dimensions change with g_causal)
  for(i in 1:n){
    dat[i,1]=as.matrix(dat[i,2:(g_causal+1)])%*%beta 
  }
  
  ##Fit linear regression to find weights of PRS score
  beta_pred=data.frame(matrix(0,nrow=g_causal,ncol=1))
  
  for(i in 2:g_causal+1){
  red.data=dat[, c(1:i)]
  beta.pred[i-1,] = glm(X1 ~ ., family=binomial("logit"), data=red.data) #oops, distribution not quite right. 
  }
  
  ##calculate PRS (select number of casual SNPs used to calculate )
  riskScore <- riskScore(weights=fit, data=dat,
                         cGenPreds=select_gen+1, Type="weighted")   # add 1 to account for intercept
  
  ##Determine predicted disease status 
  #calculate PRS 0.95 threshold 
  threshold = quantile(riskScore,1-prevalence)
  pred_disease = which(riskScore >= threshold) #predicted deiseased
  true_disease = which(disease==1) #actual diseased
  
  ### PRS performance statistics ###
  # Sensitivity
  sens = sum(pred_disease %in% true_disease)/length(true_disease)
  #False Positive 
  falsepos = sum(!(pred_disease %in% true_disease))/(n-length(true_disease))
  
  
  results = list(sens=sens,
                 falsepos=falsepos)
}