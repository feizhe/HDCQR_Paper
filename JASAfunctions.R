# onesplit returns (1) set -- the selected set of variables (2) fitted coefficient for each variable j
#(3) the selected samples used in this split

# logx   #### the log of observed time
# delta  #### censoring indicator 0 denotes censoring
# zmat   #### the n*p covariance matrix, the intercept is not included
# J      #### the quantile levels where the inference will be carried out
# beta0  #### initial estimates at tauL if using Zheng et al 2018
# lam1   #### initial tuning parameters
# q      #### select q many variables out if using zheng et al 2018
# taus   #### quantile grid for estimation using Zheng et al 2018
# if using sis in the first step, the 20 variables with largested coefficients will be selected



onesplit<-function(logx,delta,zmat,J,beta0,lam1,q=sqrt(sum(delta)),thr=0.1,taus = 1:4/5){
  
  sis <- function(logx, delta, zmat, taus = (1:8)/10, maxonly = TRUE){
    
    crqj <- Vectorize(function(j,logx, delta, zmat, taus){
      tt=crq(Surv(logx,delta,type = "right")~zmat[,j],taus=taus,method="PengHuang") 
      coef(tt, taus = taus)[2,]
    }, vectorize.args = c("j"))
    

    fitcoef <- foreach(j = (1:p), .combine=cbind, .packages=c("MASS","quantreg","rqPen","survival")) %dopar%
      crqj(j,logx, delta, zmat, taus)
    coef1 <- apply(abs(fitcoef), 2, quantile, probs=0.95, na.rm=T)
    
    rank1 <- rank(-coef1)
    return(rank1)
  }
  
  rqj <- function (j, logx, delta, zmat, samp1, set, taus) {
    set1<-c(j,set)
    rq2 <- crq( Surv(logx[-samp1],delta[-samp1]) ~ zmat[-samp1,set1],  method = "PengHuang")
    as.numeric(coef(rq2,taus = taus)[2,] )   
  }
  
  zmat1<-cbind(rep(1,n),zmat)
  
  samp1<-sort(sample(n,n/2))
  
  sampidx<-rep(0,n)
  sampidx[-samp1] <- 1
  
  
  out <- sis(logx[samp1], delta[samp1], zmat[samp1,], taus = J)
  set = which(out <= 20)            #### indices of selected variables
  
  coefall<-  matrix(NA,nrow = 1+p,ncol = length(taus))
  
  rq2 <- crq( Surv(logx[-samp1],delta[-samp1]) ~ zmat[-samp1,set],  method = "PengHuang")
  coefall[c(1,1+set),] <- coef(rq2,taus = taus)
  
  tmp1 <- foreach(j = (1:p)[-set], .combine=rbind, .packages=c("MASS","glmnet","quantreg","rqPen","survival")) %dopar%
    rqj(j, logx, delta, zmat, samp1, set, taus)
  coefall[-c(1,1+set),] <- tmp1
  
  return(returnlist=list("set" = set, "bhat" = coefall, "sample" = sampidx ))
}




# variance to calculate variance following the formula from the paper
Vsub<-function(beta,Ycount,n,n1,B, unbiased = TRUE){
  v1 = (sum(cov(beta,t(Ycount))^2))*(n-1)*n/(n-n1)^2
  ifelse(unbiased, v1 - n*n1/B/(n-n1)*var(beta), v1)
}


# function to read values in the output of onesplit
extr<-function(fit1,k,pos){
  temp1<-fit1[[k]][[pos]]
}

# summarize the results from nloop splits, return (1) the average betas over splits, (2)the standard deviations of betahat
  #### (3) selected frequency
# fit1, the results from nloop onesplit 
sumfn<-function(fit1){
  B = length(fit1)                                            #### B splits
  
  BETA1<-sapply(1:B,extr,fit1=fit1,pos=2,simplify = "array")  #### obtained Beta's  (p+1)*length(taus)*B
  SET1<-lapply(1:B,extr,fit1=fit1,pos=1)                      #### obtained selected sets of variables
  Ycount<-sapply(1:B,extr,fit1=fit1,pos=3)                    #### selected sets of samples in nloops samples
  
  avgsel<- length(unlist(SET1))/B                             #### average length of selected variables
  
  bhat = apply((BETA1), c(1,2),mean,na.rm=T)                  #### fitted quantile coefficients   (p+1)*length(taus)
  
  vbetas1<-apply(BETA1, c(1,2), Vsub,Ycount=(Ycount),n=n,n1=n/2,B=B)  #### covariance matrix estimation
  sds1<- suppressWarnings(sqrt(vbetas1))                      #### sqrt of vbetas1
  
  #### selection relative frequency of each variable
  temptab<-table(unlist(SET1))
  sel_freq<-rep(0,p)
  sel_freq[as.numeric(names(temptab))]<-temptab/B
  
  return(cbind(bhat,sds1,c(0,sel_freq)) )                     #### (p+1)*(;ength(taus)*2+1)
}
