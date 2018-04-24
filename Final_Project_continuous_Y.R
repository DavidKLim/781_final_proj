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

PRS <- function(n,g_causal,g, b){
  #n= training sample size, g_causual= number of casual SNPs, g= total number of SNPS, 
  #create training set X (columns: Y, Causual SNP genotypes, non-causual SNP genotypes)
  dat = data.frame(matrix(0,nrow=n,ncol=g+1))
  
  #number of non-causal SNPs
  g_non = g - g_causal
  
  #simulate MAF for casual SNPs and non-casual SNPs
  p = c(runif(g_causal,0.1,0.2),runif(g_non,0.4,0.5))
  
  ##simulate genotype dataset for SNPs in training set X
  for(i in 1:n){
    for(j in 1:g){
      dat[i,j+1] = rbinom(1,2,p[j])
    }
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
    for(j in 1:g){
      datz[i,j+1] = rbinom(1,2,p[j])
    }
  }
  
  ##calculate PRS of the Z 
  riskScore = datz[,2:(g+1) ]%*%beta.pred
  
  #calculate true Y of Z 
  datz[,1]=as.matrix(datz[,2:(g_causal+1)])%*%beta 
  
  
  ######Correlation between true Y of test set and PRS score ####
  Corr=cor(datz[,1], riskScore)
  return(Corr)
}

######################## Step 6-7: Simulations ###################################
### Adjustments: 
#1. Change in Training sample size (n): 100, 500, 1000, 2500, 5000
#2. Change in # of SNPs simulated (g): 100, 1000, 100000, 50000
#3. Change in # of causal SNPs selected (g_causal): 0.05*g, .10*g, .25*g, .50*g, 0.75*g

#simulate betas per g_causual cateogry before simulation

## Simulation Table: 1 x 2
n = c(100,500,1000,2500,5000)
g = c(100,1000,10000,50000)
gcaus = c(0.05,0.1,0.25,0.5,0.75)

results = matrix(0,nrow=length(n)*length(g)*length(gcaus),ncol=4)
results[,1]=rep(n,times=length(g)*length(gcaus))
results[,2]=rep(rep(g,times=length(gcaus)),each=length(n))
results[,3]=rep(gcaus,each=length(n)*length(g))

sim=100

library(parallel)
no_cores=12

# # Non-parallelized
# for(i in 1:nrow(results)){
#   gcaus = results[i,3]*results[i,2]
#   corr = rep(0,sim)
#   b1 = runif(gcaus,0,1)
#   for(s in 1:sim){
#     X = PRS(n=results[i,1],g=results[i,2],g_causal=gcaus,b=b1)
#     corr[s] = X
#     print(paste("conditions:",results[i,],"Simulation s: corr =",corr[s]))   # Tracks progress. Can comment out
#   }
#   results[i,4] = mean(corr)
# }

par_sim_run = function(i){
  gcaus = results[i,3]*results[i,2]
  corr = rep(0,sim)
  b1 = runif(gcaus,0,1)
  for(s in 1:sim){
    start=Sys.time()
    X = PRS(n=results[i,1],g=results[i,2],g_causal=gcaus,b=b1)
    end=Sys.time()
    corr[s] = X
    print(paste("time elapsed:",start-end))
    #print(paste("conditions:",results[i,],"Simulation s: corr =",corr[s]))   # Tracks progress. Can comment out
  }
  results[i,4] = mean(corr)
  print(results[i,4])
}

cl<-makeCluster(no_cores,outfile="781_out.txt")
clusterExport(cl=cl,varlist=c(ls(),"PRS","par_sim_run"))

results[,4] = parSapply(cl, 1:nrow(results), par_sim_run)

stopCluster(cl)

save(results,file="781_res.out")


## Simulation Table: 1 x 3

## Simulation Table: 1 x 4

## Simulation Table 2 x 3 

## Simulation Table 2 x 4

## Simulation Table 3 x 4 ?? 