
###############AR1 correlation structure

HH <- abs(outer(1:p, 1:p, "-")) 
Sigma <- rho^HH
rm(HH)

Sigma[1:5,1:5]

#### Grid and other issues
J=seq(0,1,length.out=n/5);
J = J[which((J>tauL)*(J<=tU)!=0)];

length(J)

#### coeff
alpha=rep(0,p);
alpha[20]=1.5;
alpha[40]=2;
alpha[60]=2.5;

gamma.1=rep(0,p);
gamma.2=rep(0,p);
gamma.1[1]= g1 = 1;
gamma.2[10]= g2 = .7;


TrueBeta=matrix(rep(alpha, length(J)), length(J), p, byrow=T)+
  matrix(rep(gamma.1, length(J)), length(J), p, byrow=T)*qeffect.1(J)+
  matrix(rep(gamma.2, length(J)), length(J), p, byrow=T)*qeffect.2(J);

dim(TrueBeta)
beta=rep(0,p);
tmp=union(which(gamma.1!=0), which(gamma.2!=0));
beta[sort(union(tmp,which(alpha!=0)))]=1;
s1 = which(beta!=0)
#head(TrueBeta[,s1])
#tail(TrueBeta[,s1])

plot(J, TrueBeta[,s1[1]], type = "l", ylim = c(-5,4))
lines(J, TrueBeta[,s1[2]], lty = 2)

taus = J

ntau = ceiling(length(taus)/10)
tau.idx <- 10*(1:ntau)-9
