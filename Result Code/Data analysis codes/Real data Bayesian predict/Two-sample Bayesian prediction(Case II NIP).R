#Bayesian two-sample prediction for the real data in Case II under NIP. 
completedata <- c(0.1788,0.2892,0.3300,0.4152,0.4212,0.4560,0.4880,0.5184,0.5196,0.5412,0.5556,0.6780,0.6864,0.6864,0.6888,0.8412,0.9312,0.9864,1.0512,1.0584,1.2792,1.2804,1.7340)

datareal <- c(0.1788,0.2892,0.33,0.4152,0.4212,0.456,0.5184,0.5412,0.6888,0.8412,0.9312,0.9864,1.0512,1.2792)
R1 <- rep(0,14)
R1[c(4,8,14)] <- c(3,4,2)
#case2
n=length(completedata)
m=length(datareal)
c=11
Tim=1.2  
Datah<-hybriddata(Tim=1.2,m=14,c=11,datareal,R1)

a <-0
b <- 0
X=Datah$X
JX =Datah$JX
RXJ=Datah$RXJ
TimX=Datah$TimX
library(ars)
Nsim11 <- 1000
B <- 100
premiu <- array(dim=Nsim11)
prelamda <- array(dim=Nsim11)
if(Datah$type==2){
  lastR <- array(dim=Datah$JX+1)
  lastR[1:length(Datah$X)] <- R1[1:length(Datah$X)]
  lastR[Datah$JX+1] <- n-sum(R1[1:Datah$JX])-Datah$JX
  Xcuo <- array(dim=length(Datah$X)+1)
  Xcuo <- Datah$X
  X=Xcuo
  X[length(Datah$X)+1] <- Datah$TimX
}else{
  lastR <- Datah$lastR
}


library(ars)
l <-5
N2 <-10
S <-c(0,1,2,2,0)

Nsim11 <- 1000
B <- 100
premiu <- array(dim=Nsim11)
prelamda <- array(dim=Nsim11)
ysln <- array(dim=l)
for (s in 1:l) {
  cnsmatrix <- array(dim=s)
  for (i in 1:s) {
    if(i==1){
      cnsmatrix[1]<- N2
    }else{
      cnsmatrix[i] <- N2-sum(S[1:(i-1)])-i+1
    }
  }
  Cns <- prod(cnsmatrix)
  
  qiuhe <- array(0,dim=Nsim11)
  for (h in 1:Nsim11) {
    premiu[h]<- ars(1,logPaimiu,dlogPaimiu,x=0.1,m=1,lb=TRUE,xlb=0.01,ub=TRUE,xub=X[1])
    Gmiu<- sum((R1[1:JX]+1)*(X[1:JX]-premiu[h])^2)+RXJ*(TimX-premiu[h])^2
    prelamda[h] <- rgamma(1,shape = (JX+a),rate = (b+Gmiu))
    for (q2 in 0:(s-1)) {
      if(q2==(s-1)){
        Mqs <- N2
      }else{
        Mqs <- N2-sum(S[1:(s-q2-1)])-s+q2+1
      }
      if(q2!=0){
        part1 <- array(dim=q2)
        for (u in 1:q2) {
          part1[u] <- sum(S[(s-q2):(s-q2+u-1)]+1)
        }
      }else{
        part1 <- 1
      }
      if(q2!=(s-1)){
        part2 <- array(dim=s-q2-1)
        for(u in 1:(s-q2-1)){
          part2[u] <- sum(S[u:(s-q2-1)]+1)
        }
      }else{
        part2 <- 1
      }
      cqs1 <- (-1)^q2*(1/(prod(part1)*prod(part2)))
      
      fvalue <- function(x){
        return(exp(-Mqs*prelamda[h]*(x-premiu[h])^2))
      }
      intvalue <- integrate(fvalue,0,Inf)$value
      qiuhe[h] <- qiuhe[h]+cqs1*(1/Mqs)*intvalue*Cns
    }
  }
  ysln[s] <- mean(qiuhe[B:Nsim11])  #0.3232428 0.4546138 0.5896868 0.7446336 1.1330942
}
twosampleprevalue <- ysln

a1 <- c(0.975,0.025)
function4 <- function(t,a1,s){
  cnsmatrix <- array(dim=s)
  for (i in 1:s) {
    if(i==1){
      cnsmatrix[1]<- N2
    }else{
      cnsmatrix[i] <- N2-sum(S[1:(i-1)])-i+1
    }
  }
  Cns <- prod(cnsmatrix)
  
  qiuhe <- array(0,dim=Nsim11)
  for(h in 1:Nsim11){
    premiu[h]<- ars(1,logPaimiu,dlogPaimiu,x=0.1,m=1,lb=TRUE,xlb=0.01,ub=TRUE,xub=X[1])
    Gmiu<- sum((R1[1:JX]+1)*(X[1:JX]-premiu[h])^2)+RXJ*(TimX-premiu[h])^2
    prelamda[h] <- rgamma(1,shape = (JX+a),rate = (b+Gmiu))
    for (q2 in 0:(s-1)) {
      if(q2==(s-1)){
        Mqs <- N2
      }else{
        Mqs <- N2-sum(S[1:(s-q2-1)])-s+q2+1
      }
      if(q2!=0){
        part1 <- array(dim=q2)
        for (u in 1:q2) {
          part1[u] <- sum(S[(s-q2):(s-q2+u-1)]+1)
        }
      }else{
        part1 <- 1
      }
      if(q2!=(s-1)){
        part2 <- array(dim=s-q2-1)
        for(u in 1:(s-q2-1)){
          part2[u] <- sum(S[u:(s-q2-1)]+1)
        }
      }else{
        part2 <- 1
      }
      cqs1 <- (-1)^q2*(1/(prod(part1)*prod(part2)))
      
      qiuhe[h] <- qiuhe[h]+Cns*cqs1*(1/Mqs)*exp(-prelamda[h]*Mqs*(t-premiu[h])^2)
    }
  }
  return(mean(qiuhe)-a1)
}

twosampleET <- matrix(nrow=5,ncol=2)
for (a1 in c(0.975,0.025)) {
  for (s in 1:5) {
    if(a1==0.975){
      tryCatch({twosampleET[s,1] <- uniroot(function4,c(0,1000),a1=a1,s=s)$root },error=function(e){
        we=1
      })
    }else{
      twosampleET[s,2] <- uniroot(function4,c(0,1000),a1=a1,s=s)$root
    }
  }
}

xtest <- seq(0.08,0.1,0.00001)
ytest <- array(dim=length(xtest))
for (i1 in 1:length(xtest)) {
  ytest[i1] <- function4(xtest[i1],0.975,s=1)
}
plot(xtest,ytest,type = 'l')
abline(0, 0, untf = FALSE,col='red')
xtest[which.min(abs(ytest))]
twosampleET[1,1]<- xtest[which.min(abs(ytest))]

ff <- list(twosampleprevalue=twosampleprevalue,twosampleET=twosampleET)
library(aqp)
sink("E:\\progress result1119\\twosampleprecase2NIP.txt")
ff
sink()
