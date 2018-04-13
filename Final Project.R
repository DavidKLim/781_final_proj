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


set.seed(9)
# inputs: n, g_causal, g, prevalence

PRS <- function(n,g_causal,g,prevalence){
  g_non = g - g_causal
  disease = rbinom(n,1,prevalence)
  
  p = c(runif(g_causal,0.1,0.15),runif(g_non,0.3,0.5))
  
  dat = data.frame(matrix(0,nrow=n,ncol=g+1))
  dat[,1] = disease
  
  for(i in 1:n){
    for(j in 1:g_causal){
      if(disease[i]==1){
        dat[i,j+1] = rbinom(1,2,p[j])     # 0.8 for higher frequency of geno = 2?
      } else {
        dat[i,j+1] = rbinom(1,2,1-p[j])     # is this p = 0.5? H-W-E?
      }
    }
    
    for(j in (g_causal+1):g){
      dat[i,j+1] = rbinom(1,2,p[j])
    }
  }
  
  fit = glm(X1 ~ ., family=binomial("logit"), data=dat)
  
  riskScore <- riskScore(weights=fit, data=dat,
                                    cGenPreds=c(2:4), Type="weighted")   # make way to SELECT cGenPreds from fit (significant betas)
  
  threshold = quantile(riskScore,1-prevalence)
  pred_disease = which(riskScore >= threshold)
  true_disease = which(disease==1)
  
  sens = sum(pred_disease %in% true_disease)/length(true_disease)
  falsepos = sum(!(pred_disease %in% true_disease))/(n-length(true_disease))

  results = list(sens=sens,
                 falsepos=falsepos)
}
sim = 100
for(s in 1:sim){
  X = PRS(n=1000,g_causal=3,g=10,prevalence=0.2)
}

