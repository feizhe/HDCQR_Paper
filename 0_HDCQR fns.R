
####  y: the observed variable( log(min(T,C)));
####  x: the covariates;
####	delta: the censor indicator;
####  JJ: the tau grid;
####  beta0: the initial estimator at tau0;
####  lambda: the tuning parameter

quantproc=function(y,x,delta,JJ,beta0,lambda,incr=0,	tol=1e-4){
  
  pp=dim(x)[2];                                                        #### the number of covariates
  tmpbeta=matrix(0,length(JJ),pp);                                     #### the coefficient matrix
  tmpbeta[1,]=beta0;
  
  sto_weights=matrix(0,length(JJ),dim(x)[1]);                          #### stochastic weights
  rproc=matrix(0,length(JJ),dim(x)[1]);                                #### the indictor (y>=x%*%betahat);
  sto_weights[1,]=JJ[1]*2;                                             #### the initial weight 2tau0
  rproc[1,]=1*(y>=x%*%beta0);                                          #### the initial indicator y>=x%*%betahat0;
  
  augy1=y[which(delta==1)];                                            #### the observed event time   
  augx1=x[which(delta==1),];                                           #### the corresponding covariates
  
  augx2=-apply(augx1,2,sum);                                           #### the 2nd part in the objective function                        
  
  augy=c(augy1, 1e+4, 1e+4);
  
  for(s in 2:length(JJ)){
    tuning=lambda+(s-1)*incr
    Hm = H(JJ[s])-H(JJ[s-1]);                                        #### H(tau[s])-H(tau[s-1]) 
    sto_weights[s,]=sto_weights[s-1,]+2*Hm*rproc[s-1,];            #### update the stochastic weight for tau[s]
    augx3=sto_weights[s,]%*%x;                                     #### the 3rd part in the objective function
    augx=rbind(augx1,augx2,augx3);
    tmpb=coef(rq(augy~0+augx,method="lasso",tau=JJ[s],lambda=rep(tuning,pp)));  #### quantile fit at tau[s];
    tmpbeta[s,]=tmpb*(abs(tmpb)>=tol);                             #### hard threshholding;
    rproc[s,] = 1*(y>x%*%tmpbeta[s,]);                               #### update the indicator y>=x%*%betahat at tau[s];	
  }
  return(tmpbeta);
}

####  y: the observed variable( log(min(T,C)));
####  x: the covariates;
####	delta: the censored indicator;
####  JJ: the tau grid;
####  beta: the estimates at JJ grid


#### function to calculate the loss at quantile tau
rho=function(err,tau){
  (err>0)*tau*err-(1-tau)*(err<0)*err;
}

calesterr=function(yP,xP,deltaP,JJ,beta){
  
  fit=xP%*%t(beta); 
  testerr=rep(0,length(JJ));
  
  sto_weights=matrix(0,length(JJ),dim(xP)[1]);                          #### stochastic weights
  rproc=matrix(0,length(JJ),dim(xP)[1]);                                #### the indictor (y>=x%*%betahat);
  sto_weights[1,]=JJ[1];                                             #### the initial weight 2tau0
  rproc[1,]=1*(yP>=fit[,1]); 
  
  countP=(1-rproc[1,])*deltaP;
  rhs=rep(JJ[1],length(countP));
  testerr[1]= mean( abs(countP-rhs));
  
  for(s in 2:length(JJ)){
    Hm=H(JJ[s])-H(JJ[s-1]);                                        #### H(tau[s])-H(tau[s-1]) ]']
    
    sto_weights[s,]=sto_weights[s-1,]+Hm*rproc[s-1,];  
    rproc[s,]=1*(yP>=fit[,s]);           #### update the stochastic weight for tau[s]
    countP=(1-rproc[s,])*deltaP;
    testerr[s]= mean(abs(countP-sto_weights[s,]));
  }
  return (testerr);
}

####  y: the observed variable( log(min(T,C)));
####  x: the covariates;
####	delta: the censored indicator;
####  JJ: the tau grid;
####  beta0: the initial estimator at tau0;
####  len: the length of the finer tuning parameter sequences
####  K: the K fold cross validation

rqUTlassoH_CV=function(y,x, delta, JJ, beta0,lambda, tau=0.5, len=sqrt(nrow(x)), K=5,tol=1e-4){
  # the tolerance
  ok<-complete.cases(x,y);
  x<-x[ok,]                            # get rid of na's
  y<-y[ok]                             # since regsubsets can't handle na's
  pp=ncol(x);                         
  nn=nrow(x);
  x<-as.matrix(x);                     # in case x is not a matrix
  
  ##sto_weights=matrix(0,length(JJ),nn);
  ##rproc=matrix(0,length(JJ),nn);
  ##sto_weights[1,]=JJ[1]*2;
  ##rproc[1,]=1*(y>=x%*%beta0);
  
  size <- (nn) %/% K;
  rdm <- runif(nn);
  ranked <- rank(rdm);
  block <- (ranked-1) %/% size+1;
  block <- as.factor(block);
  
  lam=seq(lambda/2, 5*lambda, length.out=len);
  
  
  cv=rep(Inf,len);
  earlyStop = 0
  for (i in 1: len){
    ka=lam[i];
    testerr=rep(0,K);             		
    for (k in 1:K) {
      
      yP=y[block==k];         #### the prediction response variables in k fold 
      xP=x[block==k,];	      #### the prediction covariates matrix in k fold 
      deltaP=delta[block==k]; #### the prediction delta in k fold 
      
      yT=y[block!=k];         #### the test response variable in k fold 
      xT=x[block!=k,];        #### the test covariates matrix in k fold 
      deltaT=delta[block!=k]; #### the test delta in k fold 
      
      tmpbeta=quantproc(yT,xT,deltaT,JJ,beta0,ka);  #### the estimates at JJ in k fold 
      
      
      testerr[k]=sum(calesterr(yP,xP,deltaP,JJ,tmpbeta));
      
    }
    
    cv[i]=sum(testerr);
    #print(c(ka, cv[i]))
    if (i>1 ) earlyStop = ifelse(cv[i] > cv[i-1], earlyStop+1, 0)
    if (earlyStop >2) break
  }
  step.cv=which.min(cv);
  tmpbeta=quantproc(y,x,delta,JJ,beta0,lam[step.cv]); 
  obj=list(beta=tmpbeta, lamb=lam[step.cv]);
  return(obj);
}



rqOAH=function(y,x,tau=0.5, len=2*sqrt(nrow(x)), K=5, tol=1e-4, traceprint=FALSE){
  #require(quantreg);
  ok<-complete.cases(x,y);
  x<-x[ok,]                            # get rid of na's
  y<-y[ok]                             # since regsubsets can't handle na's
  pp=ncol(x);                         
  nn=nrow(x);
  x<-as.matrix(x);                     # in case x is not a matrix
  
  ######## K-fold CV 
  size <- nn %/% K;
  rdm <- runif(nn);
  ranked <- rank(rdm);
  block <- (ranked-1) %/% size+1;
  block <- as.factor(block);
  
  
  #### 	1st level search (a bold search)
  comp=Inf;                             # Cross-Validation store
  count=0;					  # 	
  ptr1=0;                               # the smallest CV
  ptr2=0;          				  # the second smallest CV
  Ada=TRUE;
  #ka=1/nn/2;
  ka=ifelse(p<=nn, p/nn, log(p/nn));
  init_lam=NULL; 
  init_cv=NULL;
  while(Ada){
    #ka=ka*2;
    ka=ka*log(p);
    count=count+1;
    init_lam=c(init_lam,ka);
    tmpb=coef(rq(y~0+x,method="lasso",tau=tau,lambda=rep(ka,pp)));
    if(sum(abs(tmpb[-1]))<=tol) {Ada=FALSE;}
    if (traceprint) print(c(ka, sum(abs(tmpb[-1]))))
    
    testerr=rep(0,K);
    for (k in 1:K) {
      xT=x[block!=k,];
      xP=x[block==k,];
      yT=y[block!=k];
      yP=y[block==k];
      
      fitb=coef(rq(yT~0+xT,method="lasso",tau=tau,lambda=rep(ka,pp)));
      fitb=fitb*(abs(fitb)>=tol);
      err = yP-xP%*%fitb;
      testerr[k]=sum((err>0)*tau*err - (1-tau)*(err<0)*err); #### why (1-tau)?
    }
    cv_testerr=sum(testerr);
    init_cv=c(init_cv,cv_testerr);
    if (cv_testerr < comp) {
      comp=cv_testerr; 
      ptr2=ptr1;
      ptr1=count;
    }
  }
  
  #### 2nd search ( a finer one)
  ## construct a CV grid
  gL=min(ptr1,ptr2);             gU=max(ptr1,ptr2);
  if(gL==1) {gL=2;}              
  gridL=max(init_lam[1], init_lam[gL-1]);
  if(gU==count) {gU=count-1;}    
  gridU=min(init_lam[count], init_lam[gU+1]);
  
  lam=seq(gridL, gridU,length.out=len);
  cv1=rep(Inf,len);
  earlyStop = 0
  for (i in 1: len){
    ka=lam[i];
    testerr=rep(0,K);
    for (k in 1:K) {
      xT=x[block!=k,];
      xP=x[block==k,];
      yT=y[block!=k];
      yP=y[block==k];
      
      fitb=coef(rq(yT~0+xT,method="lasso",tau=tau,lambda=rep(ka,pp)));
      err= yP-xP%*%fitb;
      testerr[k]=sum((err>0)*tau*err-(1-tau)*(err<0)*err);
    }
    cv1[i]=sum(testerr);
    if (traceprint) print(c(ka, cv1[i]))
    if (i>1 ) earlyStop = ifelse(cv1[i] > cv1[i-1], earlyStop+1, 0)
    if (earlyStop >2) break
  }
  step.cv1=which.min(cv1);  
  
  ######## init beta by CV
  betatilde=coef(rq(y~0+x,method="lasso",tau=tau,lambda=rep(lam[step.cv1],pp)));
  betatilde=betatilde*(abs(betatilde)>=tol);
  obj=list(beta=betatilde,lambda=lam[step.cv1]);
  return(obj);
}



adquantproc=function(y,x,delta,JJ,beta0,lambda,w,	tol=1e-6){
  
  pp=dim(x)[2];                                                        #### the number of covariates
  tmpbeta=matrix(0,length(JJ),pp);                                     #### the coefficient matrix
  tmpbeta[1,]=beta0;
  
  sto_weights=matrix(0,length(JJ),dim(x)[1]);                          #### stochastic weights
  rproc=matrix(0,length(JJ),dim(x)[1]);                                #### the indictor (y>=x%*%betahat);
  sto_weights[1,]=JJ[1]*2;                                             #### the initial weight 2tau0
  rproc[1,]=1*(y>=x%*%beta0);                                          #### the initial indicator y>=x%*%betahat0;
  
  augy1=y[which(delta==1)];                                            #### the observed event time   
  augx1=x[which(delta==1),];                                           #### the corresponding covariates
  
  augx2=-apply(augx1,2,sum);                                      #### the 2nd part in the objective function                       
  augy=c(augy1, 1e+4, 1e+4);
  
  ka=lambda;
  for(s in 2:length(JJ)){
    Hm=H(JJ[s])-H(JJ[s-1]);                                        #### H(tau[s])-H(tau[s-1]) 
    sto_weights[s,]=sto_weights[s-1,]+2*Hm*rproc[s-1,];            #### update the stochastic weight for tau[s]
    augx3=sto_weights[s,]%*%x;                                     #### the 3rd part in the objective function
    augx=rbind(augx1,augx2,augx3);
    
    tka=ka/length(augy);
    tw=w[s,];
    tw=(abs(tw)<=tka)+(2.4*tka-abs(tw))*((2.4*tka-abs(tw))>0)/((2.4-1)*tka)*(abs(tw)>tka);  #SCAD weights  with a=2.4	
    
    tmpb=coef(rq(augy~0+augx,method="lasso",tau=0.5,lambda=ka*tw));  #### quantile fit at tau[s];
    tmpbeta[s,]=tmpb*(abs(tmpb)>=tol);                             #### hard threshholding;
    rproc[s,]=1*(y>x%*%tmpbeta[s,]);                               #### update the indicator y>=x%*%betahat at tau[s];	
  }
  return(tmpbeta);
}

rqACV=function(y,x, delta, JJ, pwei, lambda, len=sqrt(nrow(x)),K=1,tol=1e-4){
  pp=ncol(x);                         
  nn=nrow(x);
  
  #### pointwise weights or uniform weights	
  if(length(pwei)==pp){
    w=t(replicate(length(JJ), pwei))
  } else if(length(pwei)==pp*length(JJ)){
    w=pwei;
  } else {
    return ("invalida weights, the length of weights should be equal to # of columns of matrix or 
			   # of quantiles times # of columns of matrix");
  }
  
  w[which(w==0)]=1e-20;
  
  tmpbeta=matrix(0,length(JJ),pp);
  
  ##augy1=y[which(delta==1)];     
  ##augx1=x[which(delta==1),];  
  
  ##augx2=-apply(augx1,2,sum);
  ##augy=c(augy1, 1e+4, 1e+4);
  
  size <- (nn) %/% K;
  rdm <- runif(nn);
  ranked <- rank(rdm);
  block <- (ranked-1) %/% size+1;
  block <- as.factor(block);
  
  
  
  lam=seq(lambda/2, 8*lambda, length.out=len);
  cv=rep(Inf,len);
  earlyStop = 0
  for (i in 1: len){
    ka=lam[i];
    testerr=rep(0,K);   
    sel=rep(0,dim(x)[2]);         		
    for (k in 1:K) {
      
      yP=y[block==k];         #### the prediction response variables in k fold 
      xP=x[block==k,];	      #### the prediction covariates matrix in k fold 
      deltaP=delta[block==k]; #### the prediction delta in k fold 
      
      yT=y[block==k];         #### the test response variable in k fold 
      xT=x[block==k,];        #### the test covariates matrix in k fold 
      deltaT=delta[block==k]; #### the test delta in k fold 
      
      tka=ka/length(yT);
      tw=w[1,];
      tw=(abs(tw)<=tka)+(2.4*tka-abs(tw))*((2.4*tka-abs(tw))>0)/((2.4-1)*tka)*(abs(tw)>tka);  #SCAD weights  with a=2.4	
      tmpb=coef(rq(yT~0+xT,method="lasso",tau=JJ[1],lambda=ka*tw));
      tmpb=tmpb*(abs(tmpb)>=tol);
      
      tmpbeta=adquantproc(yT,xT,deltaT,JJ,tmpb,ka,w,tol=tol);  #### the estimates at JJ in k fold 
      sel[which(apply(abs(tmpbeta),2,max)>tol)]=1;
      testerr[k]=sum(calesterr(yP,xP,deltaP,JJ,tmpbeta));
      
    }
    
    cv[i]=sum(testerr);
    if(length(which(sel==1))>nn/4) {cv[i]=1e+8;}
    #print(c(ka, cv[i]))
    if (i>1 ) earlyStop = ifelse(cv[i] > cv[i-1], earlyStop+1, 0)
    if (earlyStop >2) break
  }
  step.cv=which.min(cv);
  
  ka=lam[step.cv];
  tka=ka/length(y);
  tw=w[1,];
  tw=(abs(tw)<=tka)+(2.4*tka-abs(tw))*((2.4*tka-abs(tw))>0)/((2.4-1)*tka)*(abs(tw)>tka);  #SCAD weights  with a=2.4	
  tmpb=coef(rq(y~0+x,method="lasso",tau=JJ[1],lambda=ka*tw));
  
  tmpb=tmpb*(abs(tmpb)>=tol);
  tmpbeta=adquantproc(y,x,delta,JJ,tmpb,lam[step.cv],w); 
  tmpbeta=tmpbeta*(abs(tmpbeta)>tol);
  obj=list(beta=tmpbeta, lamb=lam[step.cv]);
  
  return(obj);
}

#####################################
#
# Fused-HDCQR
#
##########################

progcombine<-function(){
  count<-0
  function(...) {
    count<<-count+length(list(...))
    setTxtProgressBar(pb,count)
    utils::flush.console()
    c(...)
  }
}
onesplit<-function(logx,delta,zmat,J,beta0,lam1,q=sqrt(sum(delta)),thr=0.1,taus = 1:4/5){
  
  sis <- function(logx, delta, zmat, taus = (1:8)/10, maxonly = TRUE){
    
    crqj <- Vectorize(function(j,logx, delta, zmat, taus){
      tt=crq(Surv(logx,delta,type = "right")~zmat[,j],taus=taus,method="PengHuang")  ####oracle
      #maxcoef <- max(abs(coef(tt, taus = taus)[2,]))
      coef(tt, taus = taus)[2,]
      
    }, vectorize.args = c("j"))
    
    #fitcoef = crqj(1:p, taus)
    fitcoef <- foreach(j = (1:p), .combine=cbind, .packages=c("MASS","quantreg","rqPen","survival")) %dopar%
      crqj(j,logx, delta, zmat, taus)
    #coef1 <- apply(abs(fitcoef), 2, max, na.rm=T)
    coef1 <- apply(abs(fitcoef), 2, quantile, probs=0.95, na.rm=T)
    
    rank1 <- rank(-coef1)
    return(rank1)
    #rank1[s1]
    #list("ranks_at_taus" = SISrank, "max_rank" = rank1)
  }
  
  rqj <- function (j, logx, delta, zmat, samp1, set, taus) 
  {
    set1<-c(j,set)
    
    rq2 <- crq( Surv(logx[-samp1],delta[-samp1]) ~ zmat[-samp1,set1],  method = "PengHuang")
    
    as.numeric(coef(rq2,taus = taus)[2,] )   ##rq2$sol[3,temp1]
    
  }
  zmat1<-cbind(rep(1,n),zmat)
  
  samp1<-sort(sample(n,n/2))
  
  sampidx<-rep(0,n)
  sampidx[-samp1] <- 1
  
  #onebeta=quantproc(logx[samp1],zmat1[samp1,],delta[samp1],J,beta0,lam1)  #### the estimates at JJ in k fold 
  #maxbeta<-apply(onebeta,2,function(a){max(abs(a))})
  #coefp<-maxbeta[-1]
  #q1 <- sort(abs(coefp),decreasing = T)[q]
  #set = as.numeric(which(abs(coefp) > ifelse( q1 > thr, q1, thr)))
  
  out <- sis(logx[samp1], delta[samp1], zmat[samp1,], taus = J)
  #out[s1]
  set = which(out <= 20)
  
  coefall<-  matrix(NA,nrow = 1+p,ncol = length(taus))
  
  rq2 <- crq( Surv(logx[-samp1],delta[-samp1]) ~ zmat[-samp1,set],  method = "PengHuang")
  coefall[c(1,1+set),] <- coef(rq2,taus = taus)
  
  tmp1 <- foreach(j = (1:p)[-set], .combine=rbind, .packages=c("MASS","glmnet","quantreg","rqPen","survival")) %dopar%
    rqj(j, logx, delta, zmat, samp1, set, taus)
  ##stopImplicitCluster()
  coefall[-c(1,1+set),] <- tmp1
  
  return(returnlist=list("set" = set, "bhat" = coefall, "sample" = sampidx ))
}


Vsub<-function(beta,Ycount,n,n1,B, unbiased = TRUE){
  v1 = (sum(cov(beta,t(Ycount))^2))*(n-1)*n/(n-n1)^2
  ifelse(unbiased, v1 - n*n1/B/(n-n1)*var(beta), v1)
}



extr<-function(fit1,k,pos){
  temp1<-fit1[[k]][[pos]]
}

sumfn<-function(fit1){
  B = length(fit1)
  
  BETA1<-sapply(1:B,extr,fit1=fit1,pos=2,simplify = "array")
  SET1<-lapply(1:B,extr,fit1=fit1,pos=1)
  Ycount<-sapply(1:B,extr,fit1=fit1,pos=3)
  
  avgsel<- length(unlist(SET1))/B
  
  bhat = apply((BETA1), c(1,2),mean,na.rm=T)
  
  vbetas1<-apply(BETA1, c(1,2), Vsub,Ycount=(Ycount),n=n,n1=n/2,B=B)
  sds1<- suppressWarnings(sqrt(vbetas1))
  
  temptab<-table(unlist(SET1))
  sel_freq<-rep(0,p)
  sel_freq[as.numeric(names(temptab))]<-temptab/B
  
  return(cbind(bhat,sds1,c(0,sel_freq)) )
  #return(rbind(c(int,var.int,var1.int,0),cbind(bhat,sds,sds1,sel_freq)))
}


rqj <- function (j, logx, delta, zmat, samp1, set, taus) 
{
  set1<-c(j,set)
  
  rq2 <- crq( Surv(logx[-samp1],delta[-samp1]) ~ zmat[-samp1,set1],  method = "PengHuang")
  
  as.numeric(coef(rq2,taus = taus)[2,] )   ##rq2$sol[3,temp1]
  
}


##############################
#
# MMB p-values
#
##############################

getp <- function(tmp1,set){
  pmin(as.numeric(tmp1$coefficients[-1,6])*length(set), 1)
  
}

mb_onesplit <- function(logx,delta, zmat, q=sqrt(sum(delta)), thr=0.1){
  zmat1<-cbind(rep(1,n),zmat)
  
  samp1<-sort(sample(n,n/2))
  sampidx<-rep(0,n)
  sampidx[samp1]<-1
  
  onebeta=quantproc(logx[samp1],zmat1[samp1,],delta[samp1],J,beta0,lam1)  #### the estimates at JJ in k fold 
  maxbeta<-apply(onebeta,2,function(a){max(abs(a))})
  
  coefp<-maxbeta[-1]
  q1 <- sort(abs(coefp),decreasing = T)[q]
  set = as.numeric(which(abs(coefp) > ifelse( q1 > thr, q1, thr)))
  
  rq2 <- crq( Surv(logx[-samp1],delta[-samp1]) ~ zmat[-samp1,set],  method = "PengHuang", taus = taus)
  tmp1 <- summary(rq2, taus = taus, se="boot")
  ##tmp1[[1]]$coefficients
  
  pmat <- matrix(1, nrow = p, ncol = length(taus))
  pmat[set,]<- sapply(tmp1, getp, set=set)
  
  return(pmat)
}

qj <- Vectorize(function(pvec, gamma = 0.05){
  min(1, quantile(pvec/gamma, probs=gamma))
}, vectorize.args = "gamma")

Pj <- function(pvec, gamma_seq = .05*(1:20)){
  min(1, (1-log(min(gamma_seq))))*min( qj(pvec, gamma_seq) )
}

#############################
#
# Zheng et al 2018 estimates
#
############################

zq_est <- function(logx,delta, zmat, q=sqrt(sum(delta)), thr=0.1){
  zmat1<-cbind(rep(1,n),zmat)
  
  rqbet=rqOAH(logx,zmat1,tau=tauL, K=5); 
  lam1 = rqbet$lambda
  beta0 = rqbet$beta
  print(lam1)
  #tmpbeta = quantproc(logx,zmat1, delta,J,beta0,lam1); 
  #UW2=apply(abs(tmpbeta[-1,]),2,max);
  tmp=rqUTlassoH_CV(logx,zmat1, delta, J, beta0, lam1);
  betatilde_UT=tmp$beta;
  lam2 <- tmp$lamb
  print(lam2)
  UW=apply(abs(betatilde_UT[-1,]),2,max);
  tmp=rqACV(logx,zmat1, delta, J,UW,lam1);
  #betahat_AUT=tmp$beta;
  
  obj=list(beta=tmp$beta, lamb=tmp$lamb);
  return(obj);
}

#------------------------------------------------------------------------#
### Compute Quantile Adaptive sure Independence Screening (QaSIS)
## for uncensored dataset
QaSIS<-function(y,x,tau)
{
  p=dim(x)[2]
  n=length(y)
  fit<-numeric(p)
  y<-y-quantile(y,tau) #centered y
  x<-scale(x)
  for(j in 1:p){
    x0<-x[,j]
    knots=quantile(x0,c(1/3,2/3))
    a=bs(x0, knots=knots,degree=1)
    b=rq(y~a,tau=tau)
    fit[j] <-sum((b$fitted)^2)/n #avg of f^2 for each j
  }
  return(fit)
}


QaSIS.surv = function(x,time,delta,tau)
{
  N = length(time)
  p = ncol(x)
  survy=survfit(Surv(time,delta)~1)
  medy=survy$time[min(which(survy$surv<tau))]
  
  survc=survfit(Surv(time,1-delta)~1)
  wt=rep(0,N)
  for(i in 1:N)
  {
    k=max(which(round(survc$time,4)<=round(time[i],4)))
    wt[i]=delta[i]/survc$surv[k]
  }
  
  fit=numeric(p)
  for(k in 1:p)
  {
    pix=bs(x[,k],df=3)
    betac=rq(time~pix, tau=tau, weight=wt)
    fit[k]=mean((predict(betac)-medy)^2)
  }
  return(fit)
}


##### functions for QA-SIS and SIS

maxwhich <- Vectorize(function(a, vec){
  max(which(vec<= a))
}, vectorize.args = c("a"))

minwhich <- Vectorize(function(a, vec){
  idx = ifelse(min(which(vec< a)) == Inf, length(vec), min(which(vec< a)) )
}, vectorize.args = c("a"))

fitj <- Vectorize(function(j,taus,wt, medy, sp = FALSE){
  if(sp) pix = bs(zmat[,j],df=3)
  else pix =  zmat[,j]
  betac=try(rq(logx~pix, tau=taus, weight=wt))
  #if (length(betac) > 1)  
  rowMeans((t(predict(betac)) - medy)^2 )
}, vectorize.args = c("j"))

qasis <- function(logx, delta, zmat, taus = (1:8)/10, sp=TRUE){
  
  survy=survfit(Surv(logx,delta)~1)
  medy=survy$time[minwhich(taus, survy$surv)]
  
  survc=survfit(Surv(logx,1-delta)~1)
  kvec <- maxwhich(logx, survc$time)
  wt = delta/survc$surv[kvec]
  
  fit= fitj(1:p, taus, wt, medy, sp=sp)
  #QA_rank <- apply(fit, 1, order, decreasing = TRUE)
  colnames(fit) <- 1:p
  rownames(fit) <- round(taus,3)
  return(fit)
}

########################
#
# Fan et al 2009 SIS
#
###############################

crqj <- Vectorize(function(j,taus){
  tt=crq(Surv(logx,delta,type = "right")~zmat[,j],taus=taus,method="PengHuang")  ####oracle
  #maxcoef <- max(abs(coef(tt, taus = taus)[2,]))
  coef(tt, taus = taus)[2,]
  
}, vectorize.args = c("j"))

sis <- function(logx, delta, zmat, taus = (1:8)/10, maxonly = TRUE){
  
  fitcoef = crqj(1:p, taus)
  #coef1 <- apply(abs(fitcoef), 2, max, na.rm=T)
  coef1 <- apply(abs(fitcoef), 2, quantile, probs=0.95, na.rm=T)
  
  #SISrank <- apply(abs(fitcoef), 1, order, decreasing = TRUE)
  #SISrank <- apply(-abs(fitcoef), 1, rank)
  #colnames(SISrank) <- round(taus,3)
  #rownames(SISrank) <- 1:p
  #SISrank[s1,]
  
  rank1 <- rank(-coef1)
  return(rank1)
  #rank1[s1]
  #list("ranks_at_taus" = SISrank, "max_rank" = rank1)
}


##############################
#
# Define phi functions
#
##############################


a=48;
b=3;
devfunc=function(x) {
  1/dnorm(qnorm(x))-a*(2*x-1);
}
cross.x=uniroot(devfunc, c(0.5, 0.9))$root;
dy=qnorm(cross.x)-a*(cross.x-1/2)^2

c.1=1/dnorm(qnorm(0.9));
d.1=qnorm(0.9)-dy-0.9*c.1;
c.2=1/dnorm(qnorm(0.1));
d.2=qnorm(0.1)+dy-0.1*c.2;


qnapprox.1=function(x) {
  1*(b*((a*(x-1/2)^2)*as.numeric(x<cross.x&x>0.5)+(qnorm(pmin(x, 0.9))-dy)*as.numeric(x>=cross.x&x<0.9)+(c.1*x+d.1)*as.numeric(x>=0.9)));
}

qnapprox.2=function(x) {
  1*(b*((-a*(x-1/2)^2)*as.numeric(x>(1-cross.x)&x<0.5)+(qnorm(pmax(x, 0.1))+dy)*as.numeric(x<=(1-cross.x)&x>0.1)+(c.2*x+d.2)*as.numeric(x<=0.1)));
}

qeffect.1=function(x) {
  qnapprox.1(x+0.05);
}

qeffect.2=function(x) {
  qnapprox.2(x-0.05);
}

H<-function(x){
  return (-log(1-x));
}