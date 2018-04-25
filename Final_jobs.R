
PRS_run = function(n = c(50,100,250,500,1000),g = c(100,500,1000,2500,5000,10000,25000,50000,100000),gcaus = c(0.05,0.1,0.25,0.5,0.75,1)){
  setwd("/netscr/deelim/out")
  for(i in 1:length(n)){for(j in 1:length(g)){for(k in 1:length(gcaus)){{
    cmd = rep(0, 5)
    cmd[1] = "unlink('.RData') \n source('PRS.R') \n"
    cmd[2] = sprintf("gcausal = ceiling(%d * %f) \n",g[j],gcaus[k])
    cmd[3] = "b1=runif(gcausal,0,1) \n"
    cmd[4] = sprintf("X = PRS(n=%d,g_causal=gcausal,g=%d,b=b1) \n",n[i],g[j])      # BRCA run here #
    out = sprintf("/netscr/deelim/final_%d_%d_%f",n[i],g[j],gcaus[k])
    out2 = sprintf("%s.out", out)
    cmd[5] = sprintf("save(X, file = '%s')", out2)
    cmdf = paste(cmd, collapse = "")
    write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)
    run = sprintf("bsub -M 8 -q week -o /netscr/deelim/dump R CMD BATCH %s", out)
    system(run)
  }}}}
}

PRS_collect = function(n = c(50,100,250,500,1000),g = c(100,500,1000,2500,5000,10000,25000,50000,100000),gcaus = c(0.05,0.1,0.25,0.5,0.75,1)){
  setwd("/netscr/deelim/out")
  tab = matrix(0,nrow=length(n)*length(g)*length(gcaus),ncol=4)
  index=0
  for(i in 1:length(n)){for(j in 1:length(g)){for(k in 1:length(gcaus)){{
    index=index+1
    out = sprintf("/netscr/deelim/final_%d_%d_%f",n[i],g[j],gcaus[k])
    out2 = sprintf("%s.out",out)
    print(out2)
    if(!file.exists(out2)) next
    print(out)
    load(out2)
    tab[index,4] = X
    tab[index,1] = n[i]
    tab[index,2] = g[j]
    tab[index,3] = gcaus[k]
    print(tab[index,])
  }}}}
  
  save(tab,file="781_final.out")
}
