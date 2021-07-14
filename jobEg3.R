rm(list = ls())


##### 

## Start with smaller values of n, p, and B to test the code, 
## for example, n=200, p=300, B=20

###

n=700          # sample size, 1000
p=2000          # number of predictors, 2000
tol=1e-4;      # tolerance
nloop = B = 500


source("0_HDCQR fns.R")

library(survival)
library(quantreg)
library(rqPen)
library(abind)
library(MASS);
library(glmnet);
library(doParallel);
library(foreach);

##############################
######## Simulation Example 3
##############################

jobname = paste0('_eg3_n=',n,'_p=',p)
wdir = getwd()

Sys.setenv(TZ="America/Los_Angeles")

print(c(n,p,nloop))
ncore = detectCores()
print(ncore)

kk = ncore-2


rho=0.3
tU=0.8;
tauL=0.1;   
source("1_modelSpec.R")
J;
ntau
TrueBeta[tau.idx,s1]
s1
plot(J, TrueBeta[,s1[1]], type = "l", ylim = c(-5,4))
lines(J, TrueBeta[,s1[2]], lty = 2)


nsim = 500
ns = 0
ORC <- SUMALL <- NULL
DAT <- list()

cl <- makeCluster(kk)    
registerDoParallel(cl)


print(getDoParWorkers())
while(ns < nsim){
  
  source("3_Eg3.R")
  
  if (ns %% 10 ==0) {
    print(paste(ns, "runs completed!"))
    save.image(paste0(wdir,"/Job",jobname,".RData"))
    
  }
}

stopCluster(cl)
stopImplicitCluster()
closeAllConnections()

