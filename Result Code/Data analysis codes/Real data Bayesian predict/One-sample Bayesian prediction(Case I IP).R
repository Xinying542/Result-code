#Bayesian one-sample prediction for the real data in Case I under IP. .
datareal <- c(0.1788,0.2892,0.33,0.4152,0.4212,0.456,0.5184,0.5412,0.6888,0.8412,0.9312,0.9864,1.0512,1.2792)
R1 <- rep(0,14)
R1[c(4,8,14)] <- c(3,4,2)

#case1
n=length(completedata)
m=length(datareal)
c=12
Tim=0.9  
Datah<-hybriddata(Tim=0.9,m=14,c=12,datareal,R1)

a <-13.6368
b <- 7.3856
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
onesampleprevalue <- vector('list',3)
for (j in c(4,8,12)){
  yirxj <- array(dim = lastR[j])
  for(i in 1:lastR[j]){
    qiuhe <- array(0,dim=Nsim11)
    cat(qiuhe,'\n')
    for (h in 1:Nsim11) {
      premiu[h]<- ars(1,logPaimiu,dlogPaimiu,x=0.1,m=1,lb=TRUE,xlb=0.01,ub=TRUE,xub=X[1])
      Gmiu<- sum((R1[1:JX]+1)*(X[1:JX]-premiu[h])^2)+RXJ*(TimX-premiu[h])^2
      prelamda[h] <- rgamma(1,shape = (JX+a),rate = (b+Gmiu))
      for (rou in 0:(i-1)) {
        omigarou <- rou+lastR[j]-i+1
        fvalue <- function(x){
          return(exp(-omigarou*prelamda[h]*(x-premiu[h])^2))
        }
        intvalue <- integrate(fvalue,X[j],Inf)
        qiuhe[h] <- qiuhe[h]+(factorial(lastR[j])/(factorial(i-1)*factorial(lastR[j]-i)))*((-1)^rou)*choose(i-1,rou)*exp(omigarou*prelamda[h]*(X[j]-premiu[h])^2)*(1/omigarou)*(intvalue$value+X[j]*exp(-prelamda[h]*omigarou*(X[j]-premiu[h])^2))
      }
    }
    yirxj[i] <- mean(qiuhe[B:Nsim11])
  }
  onesampleprevalue[[j]] <- yirxj
}
onesampleprevalue[c(4,8,12)]

preMatrixL <- matrix(ncol=4,nrow=12)
a1=0.975
for (j in c(4,8,12)) {
  if(j==8){
    for (i in 1:lastR[j]) {
      function3 <- function(x,a1){
        qiuhe <- array(0,dim=Nsim11)
        for(h in 1:Nsim11){
          premiu[h]<- ars(1,logPaimiu,dlogPaimiu,x=0.1,m=1,lb=TRUE,xlb=0.01,ub=TRUE,xub=X[1])
          Gmiu<- sum((R1[1:JX]+1)*(X[1:JX]-premiu[h])^2)+RXJ*(TimX-premiu[h])^2
          prelamda[h] <- rgamma(1,shape = (JX+a),rate = (b+Gmiu))
          for (rou in 0:(i-1)) {
            omigarou <- rou+lastR[j]-i+1 #(exp(prelamda[h]*(X[j]-premiu[h])^2))*
            qiuhe[h] <- qiuhe[h]+(factorial(lastR[j])/(factorial(i-1)*factorial(lastR[j]-i)))*((-1)^rou)*choose(i-1,rou)*(1/omigarou)*exp(-omigarou*prelamda[h]*((x-premiu[h])^2-(X[j]-premiu[h])^2))
          }
        }
        return(mean(qiuhe,na.rm = T)-a1)
      }
      tryCatch({preMatrixL[j,i]<-uniroot(function3,c(X[j],3),a1=a1)$root},error=function(e){
        we=1
      })
    }
  }else{
    for (i in 1:lastR[j]) {
      function3 <- function(x,a1){
        qiuhe <- array(0,dim=Nsim11)
        for(h in 1:Nsim11){
          premiu[h]<- ars(1,logPaimiu,dlogPaimiu,x=0.1,m=1,lb=TRUE,xlb=0.01,ub=TRUE,xub=X[1])
          Gmiu<- sum((R1[1:JX]+1)*(X[1:JX]-premiu[h])^2)+RXJ*(TimX-premiu[h])^2
          prelamda[h] <- rgamma(1,shape = (JX+a),rate = (b+Gmiu))
          for (rou in 0:(i-1)) {
            omigarou <- rou+lastR[j]-i+1 #(exp(prelamda[h]*(X[j]-premiu[h])^2))*
            qiuhe[h] <- qiuhe[h]+(factorial(lastR[j])/(factorial(i-1)*factorial(lastR[j]-i)))*((-1)^rou)*choose(i-1,rou)*(1/omigarou)*exp(-omigarou*prelamda[h]*((x-premiu[h])^2-(X[j]-premiu[h])^2))
          }
        }
        return(mean(qiuhe[B:Nsim11],na.rm = T)-a1)
      }
      tryCatch({preMatrixL[j,i]<-uniroot(function3,c(X[j],3),a1=a1)$root},error=function(e){
        we=1
      })
    }
  }
}
preMatrixL[c(4,8,12),]

preMatrixU <- matrix(ncol=4,nrow=14)
a1=0.025
for (j in c(4,8,12)) {
  if(j==8){
    for (i in 1:lastR[j]) {
      function3 <- function(x,a1){
        qiuhe <- array(0,dim=Nsim11)
        for(h in 1:Nsim11){
          premiu[h]<- ars(1,logPaimiu,dlogPaimiu,x=0.1,m=1,lb=TRUE,xlb=0.01,ub=TRUE,xub=X[1])
          Gmiu<- sum((R1[1:JX]+1)*(X[1:JX]-premiu[h])^2)+RXJ*(TimX-premiu[h])^2
          prelamda[h] <- rgamma(1,shape = (JX+a),rate = (b+Gmiu))
          for (rou in 0:(i-1)) {
            omigarou <- rou+lastR[j]-i+1 #(exp(prelamda[h]*(X[j]-premiu[h])^2))*
            qiuhe[h] <- qiuhe[h]+(factorial(lastR[j])/(factorial(i-1)*factorial(lastR[j]-i)))*((-1)^rou)*choose(i-1,rou)*(1/omigarou)*exp(-omigarou*prelamda[h]*((x-premiu[h])^2-(X[j]-premiu[h])^2))
          }
        }
        return(mean(qiuhe[B:Nsim11],na.rm = T)-a1)
      }
      tryCatch({preMatrixU[j,i]<-uniroot(function3,c(X[j],3),a1=a1)$root},error=function(e){
        we=1
      })
    }
  }else{
    for (i in 1:lastR[j]) {
      function3 <- function(x,a1){
        qiuhe <- array(0,dim=Nsim11)
        for(h in 1:Nsim11){
          premiu[h]<- ars(1,logPaimiu,dlogPaimiu,x=0.1,m=1,lb=TRUE,xlb=0.01,ub=TRUE,xub=X[1])
          Gmiu<- sum((R1[1:JX]+1)*(X[1:JX]-premiu[h])^2)+RXJ*(TimX-premiu[h])^2
          prelamda[h] <- rgamma(1,shape = (JX+a),rate = (b+Gmiu))
          for (rou in 0:(i-1)) {
            omigarou <- rou+lastR[j]-i+1 #(exp(prelamda[h]*(X[j]-premiu[h])^2))*
            qiuhe[h] <- qiuhe[h]+(factorial(lastR[j])/(factorial(i-1)*factorial(lastR[j]-i)))*((-1)^rou)*choose(i-1,rou)*(1/omigarou)*exp(-omigarou*prelamda[h]*((x-premiu[h])^2-(X[j]-premiu[h])^2))
          }
        }
        return(mean(qiuhe[B:Nsim11],na.rm = T)-a1)
      }
      tryCatch({preMatrixU[j,i]<-uniroot(function3,c(X[j],3),a1=a1)$root},error=function(e){
        we=1
      })
    }
  }
}
preMatrixU[c(4,8,12),] 

ff <- list(point=onesampleprevalue[c(4,8,12)],ETL=preMatrixL[c(4,8,12),],ETU=preMatrixU[c(4,8,12),])
library(aqp)
sink("E:\\progress result1119\\onesampleprecase1IP.txt")
ff
sink()
