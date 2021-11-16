source("JASAfunctions.R")
##### the following variables need to be provided

# logx   #### the log of observed time
# delta  #### censoring indicator 0 denotes censoring
# zmat   #### the n*p covariance matrix, the intercept is not included
# nloop  #### number of re-samples

n = nrow(zmat)  #### sample size
p = ncol(zmat)  #### dimensionality

tauL = 0.1;  #### the lower limit of quantile regions
tU = 0.8;    #### the upper limit of quantile regions
J = seq(0,1,length.out=min(n/5,100));
taus = J[which((J>tauL)*(J<=tU)!=0)] #### quantile grid for estimation

J = (1:8)*0.1;  #### the quantile levels where the inference will be carried out


beta0=lam1=NA   #### initial estimates and tuning parameters if using Zheng et al 2018

fit1<- try(foreach(i = 1:nloop,.packages=c("MASS","glmnet","quantreg","rqPen","survival","foreach","doParallel")) 
           %dopar% onesplit(logx,delta,zmat,J,beta0,lam1,taus = taus))

sum1 <- sumfn(fit1)