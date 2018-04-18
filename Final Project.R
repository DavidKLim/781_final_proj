#Bios 781 Final Project 
##### Outline of task 
##1: generate alleles to determine genotype (0,1,2) at SNPs (s) and disease status for n subjects
#(causal SNPs [s1] and regualr [s2] should be generate differnetly)
##2: run regression (linear or logit) on simulated data 
##3: select top c1 beta parameters (i.e. lowest p-value) to be number of causal SNPs (c1 < s)
##4: calcualte PRS using c1 beta parameters on each subject (n)
##5: calcualte prediction accuracy of PRS (need a threshold.. determine through 95th quantile of n subjects? )
##6: repeat r times, get average prediction accuracy for set up s, n, c1. 

##7: will need to change s (number of SNPS), n (number of subjects), and c1 (number of casual SNPs) 
#to see how it affects prediction accuracy of PRS. Not sure best way to represent this... 3 different 2x2 tables? 
#i.e. sxn, sxc1, nxc1
library(PredictABEL)
library(hdlm)

########### Example Code ################
# specify dataset with outcome and predictor variables
data(ExampleData)

# specify column numbers of genetic predictors
cGenPred <- c(11:16)

# fit a logistic regression model
# all steps needed to construct a logistic regression model are written in a function
# called 'ExampleModels', which is described on page 4-5

riskmodel <- ExampleModels()$riskModel2         # glm fit object
# compute unweighted risk scores
unweighted_riskScore <- riskScore(weights=riskmodel, data=ExampleData,
                       cGenPreds=cGenPred, Type="unweighted")

weighted_riskScore <- riskScore(weights=riskmodel, data=ExampleData,
                                  cGenPreds=cGenPred, Type="weighted")


########################Step 1-5: Creating Function###################################


set.seed(9)
# inputs: n, g_causal, g, prevalence

PRS <- function(n,g_causal,g,prevalence,k){
  if(g_causal > g){
    warning("Number of causal SNPs greater than number of SNPs. Setting them equal.")
    g_causal = g
  }
  #n= training sample size, g_causual= number of casual SNPs, g= total number of SNPS, prevalence=disease prev
  #create dataset (columns: disease status, Causual SNP genotypes, non-causual SNP genotypes)
  dat = data.frame(matrix(0,nrow=n,ncol=g+1))
  
  #number of non-causal SNPs
  g_non = g - g_causal
  
  #simulate the disease status using prevalence parameter
  disease = rbinom(n,1,prevalence)
  dat[,1] = disease
  
  #simulate MAF for casual SNPs and non-casual SNPs
  p = c(runif(g_causal,0.1,0.2),runif(g_non,0.4,0.5))
  
  ##simulate genotype dataset for SNPs
  #WLOG assume MAF associated iwth disease in causual SNPs
  for(i in 1:n){
    for(j in 1:g_causal){
      #for Casual SNPs
      if(disease[i]==1){
        dat[i,j+1] = rbinom(1,2,p[j])     
      } else {
        dat[i,j+1] = rbinom(1,2,1-p[j])     
      }
    }
      #for non-causal SNPs
    for(j in (g_causal+1):g){
      dat[i,j+1] = rbinom(1,2,p[j])
    }
  }
  
  #calculate weights for polygenic score
  fit = glm(X1 ~ ., family=binomial("logit"), data=dat)
  #fit = hdglm(X1 ~ ., family='binomial', data=dat)
  #fit = HDGLM_test(dat[,1], dat[,2:ncol(dat)], model = "logistic")
  
  #fit = fitLogRegModel(data=dat,cOutcome=1,
  #                     cNonGenPreds=c(0),cNonGenPredsCat=c(0),
  #                     cGenPreds=c(2:ncol(dat)),cGenPredsCat=c(0))     # same thing as glm() above.
                                                      # can include categorical predictors????
  
  pvals = coef(summary(fit))[-1,4]      # stores pvalues (omits intercept)
  #select_gen = which(pvals < 0.05)   # Select based on pvalue = 0.05 (is this threshold ok?? maybe bonferroni?)
  select_gen = order(pvals)[1:k]
  
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

######################## Step 6-7: Simulations ###################################
### Adjustments: 
#1. Change in Training sample size (n): 10, 100, 1000, 10000
#2. Change in # of SNPs simulated (g): 2, 10, 100, 1000 
#3. Change in # of causal SNPs selected in PRS calculation (k): 0.25*SNPs, 0.5*SNPs, 0.75*SNPs, total
#4. Change in # of causal SNPs simulated (g_causal): 0.25*SNPs, 0.5*SNPs, 0.75*SNPs
#5(?). Change in pvalue threshold to select risk SNPs?

# we could just do one of 3 or 4. not sure what the correct interpretation is...

## Simulation Table: 1 x 2

## Simulation Table: 1 x 3

## Simulation Table: 1 x 4

## Simulation Table 2 x 3 

## Simulation Table 2 x 4

## Simulation Table 3 x 4 ?? 

## Can also produce bar plots of the above table comparisons 


sim = 10
sens = rep(0,sim)
falsepos = rep(0,sim)
for(s in 1:sim){
  X = PRS(n=1000,g_causal=3,g=200,prevalence=0.2,k=3)
  sens[s] = X$sens
  falsepos[s] = X$falsepos
  print(paste("Simulation s: sens =",sens[s],"falsepos =",falsepos[s]))   # Tracks progress. Can comment out
}                   # increasing n prevents glm nonconvergence error
                    # also, i'm pretty sure n has to be > g..

mean(sens)
mean(falsepos)



