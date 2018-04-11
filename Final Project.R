#Bios 781 Final Project 
##### Outline of task 
##1: generate alleles to determine genotype (0,1,2) at SNPs (s) and disease status for n subjects
#(casual SNPs [s1] and regualr [s2] should be generate differnetly)
##2: run regression (linear or logit) on simulated data 
##3: select top c1 beta parameters (i.e. lowest p-value) to be number of causal SNPs (c1 < s)
##4: calcualte PRS using c1 beta parameters on each subject (n)
##5: calcualte prediction accuracy of PRS (need a threshold.. determine through 95th quantile of n subjects? )
##6: repeat r times, get average prediction accuracy for set up s, n, c1. 

##7: will need to change s (number of SNPS), n (number of subjects), and c1 (number of casual SNPs) 
#to see how it affects prediction accuracy of PRS. Not sure best way to represent this... 3 different 2x2 tables? 
#i.e. sxn, sxc1, nxc1


###########Example Code ################
library(PredictABEL)
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


#################################################


n = 100
n_SNP = 2
n_cSNP = 1

disease = rbinom(n,1,0.2)

SNP1 = rep(0,n)     # causal
SNP2 = rep(0,n)     # non causal

for(i in 1:n){
  if(disease[i]==1){
    SNP1[i] = rbinom(1,2,0.65)     # 0.8 for higher frequency of geno = 2?
  } else {
    SNP1[i] = rbinom(1,1,0.5)     # is this p = 0.5? H-W-E?
  }
  
  SNP2[i] = rbinom(1,2,0.5)      # is this p = 0.5? HWE?
}

dat = data.frame(cbind(disease,SNP1,SNP2))

fit = glm(disease~SNP1+SNP2, family=binomial("logit"), data=dat)

unweighted_riskScore <- riskScore(weights=fit, data=dat,
                                  cGenPreds=c(2:3), Type="unweighted")
weighted_riskScore <- riskScore(weights=fit, data=dat,
                                  cGenPreds=c(2:3), Type="weighted")




