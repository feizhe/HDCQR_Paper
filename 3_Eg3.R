# get the job id from the environment
#jobID=Sys.getenv("JOB_ID", unset = as.numeric(Sys.time()))


TimeTracker = NULL
tstart = Sys.time()
TimeTracker = rbind(TimeTracker, c("start",as.character(tstart)))


source("2_genDat.R")
#################### Time Tracker ########
tnow = Sys.time()
TimeTracker = rbind(TimeTracker, c("Data",as.character(tnow)))
#######################

#t1=proc.time()
#rqbet=rqOAH(logx,zmat1,tau=tauL, K=5); 
#print(proc.time() - t1)####  initial estimate at tauL


#################### Time Tracker ########
tnow = Sys.time()
TimeTracker = rbind(TimeTracker, c("initialize",as.character(tnow)))
#######################

tt=crq(Surv(logx,delta,type = "right")~zmat[,s1],taus=taus,method="PengHuang")  # oracle estimator
orc.coef <- coef(tt,taus = taus)
rownames(orc.coef)[-1] <- s1


plot(taus, orc.coef[1,], type = "l",col="red", ylim = c(-4,4), main = "Oracle")
for(k in 2:6){
  lines(taus, orc.coef[k,])
  
}


#t1=proc.time()
#rqbet=rqOAH(logx,zmat1,tau=tauL, K=5); 
#print(proc.time() - t1)####  initial estimate at tauL

#lam1 = rqbet$lambda;
#beta0 = rqbet$beta
#print(paste("lam = ",lam1))

beta0=lam1=NA
################ ~1100 seconds to run each dataset on a macbook pro
J = (1:8)*0.1

print("Start B splits")

t1<-proc.time()
#pb<-txtProgressBar(0,nloop,style=3)
fit1<- try(foreach(i = 1:nloop,.packages=c("MASS","glmnet","quantreg","rqPen","survival","foreach","doParallel")) 
           %dopar% onesplit(logx,delta,zmat,J,beta0,lam1,taus = taus))
#close(pb)
print(proc.time() - t1)
print(Sys.time())

if (length(fit1) != nloop) next
## 78466.398/60/60 = 21.8 hrs to run B=500, n=700, p=2000, with 22 cores
## 9638.401/60/60 = 2.67 hrs to run B=300, n=700, p=2000, with 22 cores, J = (1:7)*0.1
#################### Time Tracker ########
tnow = Sys.time()
TimeTracker = rbind(TimeTracker, c("parallel",as.character(tnow)))
#######################

sum1 <- sumfn(fit1)
print(colSums(is.na(sum1)))

####################
#
# Output results
#
######################

ns = ns + 1
jobID = ns
print(paste("ns =",ns))

print( sum1[,ncol(sum1)][1+s1])

if(FALSE){
  png(paste0("../output/eg3_n=700_p=1000/coefs_",jobID, ".png"),height = 600, width = 600)
  par(mfrow=c(2,2))
  {
    j = s1[1]
    ylim = range(c(sum1[1+j,1:length(taus)],TrueBeta[,j], orc.coef[2,] ))
    plot(taus, sum1[1+j,1:length(taus)], xlab="tau",ylab="beta", type = "l", ylim=ylim)
    lines(taus, TrueBeta[,j],lty = 2)
    lines(taus, orc.coef[2,],lty = 3)
    
    j = s1[2]
    ylim = range(c(sum1[1+j,1:length(taus)],TrueBeta[,j], orc.coef[3,] ))
    plot(taus, sum1[1+j,1:length(taus)], xlab="tau",ylab="beta", type = "l", ylim=ylim)
    lines(taus, TrueBeta[,j],lty = 2)
    lines(taus, orc.coef[3,],lty = 3)
    #lines(taus, sum1[1,1:length(taus)], col="red")
    legend("bottomright",c("Fused","Oracle"),lty=c(1,3))
    
    
    plot(taus, orc.coef[1,], type = "l",col="red", ylim = c(-4,4), main = "Oracle")
    for(k in 2:6){
      lines(taus, orc.coef[k,])
      
    }
    
    plot(taus, sum1[1,1:length(taus)], type = "l", col="red",ylim = c(-4,4), main = "Fused")
    for(k in 1+s1){
      lines(taus, sum1[k,1:length(taus)])
      
    }
    
  }
  dev.off()
  
}
dat = list(y = logx, zmat = zmat, delta = delta, lam = lam1)

#saveRDS(dat, paste0("output_n=500_p=1000/dat_",jobID,".rds"))
#saveRDS(orc.coef, paste0("output_n=500_p=1000/oracle_",jobID,".rds"))
#saveRDS(sum1, paste0("output_n=500_p=1000/sum1_",jobID,".rds"))

tmp1 <- ifelse(ns == 1, 0,1)

ORC <- abind(ORC, orc.coef, along = tmp1)
SUMALL <- abind(SUMALL, sum1, along = tmp1)
DAT[[ns]] <- dat
print(dim(SUMALL))

