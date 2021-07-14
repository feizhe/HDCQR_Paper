
zmat <- mvrnorm(n,rep(0,p),Sigma = Sigma)
#zmat <- matrix(runif(n*p), nrow = n)
zmat[,1]<-abs(zmat[,1]) +0.5
zmat[,s1[2]]<-abs(zmat[,s1[2]]) +0.5
zmat1<-cbind(rep(1,n),zmat)

e.u<-runif(n);
logt=zmat%*%alpha+zmat[,1]*g1*qeffect.1(e.u)+zmat[,s1[2]]*g2*qeffect.2(e.u);

logc = rnorm(n,0,4) + rnorm(n,-4,1) + rnorm(n,10,0.5)
delta <- as.numeric(logt <= logc)
print(round(table(delta)/n,2))
logx <- pmin(logt,logc)

dat = list("y" = logx, "zmat"=zmat, "delta" = delta)