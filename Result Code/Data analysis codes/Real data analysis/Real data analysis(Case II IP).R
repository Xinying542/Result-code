#The results of the real data analysis in Case II under IP.

completedata <- c(0.1788,0.2892,0.3300,0.4152,0.4212,0.4560,0.4880,0.5184,0.5196,0.5412,0.5556,0.6780,0.6864,0.6864,0.6888,0.8412,0.9312,0.9864,1.0512,1.0584,1.2792,1.2804,1.7340)
set.seed(19)
m <- 14
Yt <- completedata
Yr <- matrix(nrow = length(Yt),ncol=2)
Yr[,1]<- completedata
Yr[,2] <- 1:(length(Yt))
n <- length(Yt)
R1 <- rep(0,14)
R1[4] <- 3
R1[8] <- 4
R1[14] <- 2
for (i in 1:m) {
  cat(i,'\n')
  if(R1[i]!=0){
    k=n-sum(R1[1:(i-1)])
    deletesample <- sample(x = c((i+1):k),R1[i])
    cat('删除的是',Yr[deletesample,],'\n')
    Yr <- Yr[-deletesample,]
    Yr[,2] <- 1:(length(Yr[,1]))
    cat('Y的长度为',length(Yr[,1]),'\n')
    cat('Y为：',Yr[,1],'\n','\n')
  }
}

datareal <- c(0.1788,0.2892,0.33,0.4152,0.4212,0.456,0.5184,0.5412,0.6888,0.8412,0.9312,0.9864,1.0512,1.2792)
R1 <- rep(0,14)
R1[c(4,8,14)] <- c(3,4,2)


#case2
n=length(completedata)
m=length(datareal)
c=11
Tim=1.2  
Datah<-hybriddata(Tim=1.2,m=14,c=11,datareal,R1)
X=Datah$X
JX =Datah$JX
RXJ=Datah$RXJ
TimX=Datah$TimX

NRvalue<-NewtonsX(funs,Datah,theta=c(1.2,0.08),R1)

EMvalue<-EM12(R1,0.08,1.2) 

#####Bayes#####
library(ars)
library(HDInterval)

a <-13.6368
b <- 7.3856

Nsim6 <- 1000
X=Datah$X
JX=Datah$JX
RXJ=Datah$RXJ
TimX=Datah$TimX
gibbsmiu <- array(dim=Nsim6)
gibbslamda <- array(dim=Nsim6)
for (j in 1:Nsim6) {
  #生成miu值
  gibbsmiu[j]<- ars(1,logPaimiu,dlogPaimiu,x=0.1,m=1,lb=TRUE,xlb=0.01,ub=TRUE,xub=X[1])
  Gmiu<- sum((R1[1:JX]+1)*(X-gibbsmiu[j])^2)+RXJ*(TimX-gibbsmiu[j])^2
  gibbslamda[j] <- rgamma(1,shape = (JX+a),rate = (b+Gmiu))
  if(gibbsmiu[j]==0){
    gibbsmiu[j] <-  NA
    gibbslamda[j] <- NA
  }
}

selfbayesmiu <- mean(gibbsmiu[100:1000],na.rm = T)
selfbayeslamda <- mean(gibbslamda[100:1000],na.rm=T)

bayeslamdahpd <-hdi(gibbslamda[complete.cases(gibbslamda)],credMass = 0.95)
bayesmiuhpd<-hdi(gibbsmiu[complete.cases(gibbsmiu)],credMass = 0.95)
bayeslamdaCI <-CI(sort(gibbslamda[complete.cases(gibbslamda)]),tau=0.05,B=length(gibbslamda[complete.cases(gibbslamda)]))
bayesmiuCI <-CI(sort(gibbsmiu[complete.cases(gibbsmiu)]),tau=0.05,B=length(gibbsmiu[complete.cases(gibbsmiu)]))

theta0 = c(1.2,0.08) 
MLE <- NewtonsX(funs,Datah,theta0,R1)
miuhat1<-(MLE$root)[2]
lamdahat1 <- (MLE$root)[1]

FisherI <- matrix(data=c(JX/(lamdahat1^2),
                         -(2*sum(X-miuhat1)+2*sum(R1[1:JX]*(X-miuhat1))+2*RXJ*(TimX-miuhat1)),
                         -(2*sum(X-miuhat1)+2*sum(R1[1:JX]*(X-miuhat1))+2*RXJ*(TimX-miuhat1)),
                         -(-sum(1/((X-miuhat1)^2))-2*lamdahat1*JX-2*lamdahat1*sum(R1[1:JX])-2*RXJ*lamdahat1)),
                  2,2)

Taumatrix <- solve(FisherI)  

Lthreematrix <- matrix(data=c(-2*sum(1/(X-miuhat1)^3),-2*JX-2*sum(R1[1:JX])-2*RXJ,0,2*JX/lamdahat1^3),byrow=T,
                       2,2)   


Kmatrix <- c((a-1)/lamdahat1-b,0)

omiga1laglamda1 <- matrix(data=c(1,0),1,2)
omiga2laglamda1 <- matrix(data=c(0,0,0,0),2,2,byrow = T)
value1<- lindlyvalue(omiga1laglamda1,omiga2laglamda1,Taumatrix,Lthreematrix,Kmatrix)
selflindlylamda<-lamdahat1+value1

omiga1lagmiu1 <- matrix(data=c(0,1),1,2)
omiga2lagmiu1 <- matrix(data=c(0,0,0,0),2,2,byrow = T)
value2<- lindlyvalue(omiga1lagmiu1,omiga2lagmiu1,Taumatrix,Lthreematrix,Kmatrix)
selflindlymiu<-miuhat1+value2

###tk#####
thetatk0 <- c(1.2,0.08) 
parlianta<- optim(thetatk0,method='L-BFGS-B',liantafunction,lower=c(0.01,0.01),
                  hessian=T,X=X,R=R1)$par
lamdalianta <- parlianta[1]
miulianta <- parlianta[2]
mida<- liantafunction(parlianta,X,R1)
Omigamatrix <- matrix(data=c((1/n)*(JX+a+1)*(1/lamdalianta^2),
                             -(1/n)*(-sum(1/(X-miulianta)^2)-lamdalianta*(2*sum(R1[1:JX]+1)+2*RXJ)),
                             -(1/n)*(-sum(1/(X-miulianta)^2)-lamdalianta*(2*sum(R1[1:JX]+1)+2*RXJ)),
                             (1/n)*(sum((R1[1:JX]+1)*(X-miulianta)^2)+RXJ*(TimX-miulianta)^2)),
                      2,2)
OmigaValue <- det(solve(Omigamatrix))
liantaValue <- -mida

selfparliantaXlamda<- optim(thetatk0,method='L-BFGS-B',selfliantaXlamdafunction,lower=c(0.01,0.01),
                            hessian=T,X=X,R=R1)$par
selflamdaliantaXlamda <- selfparliantaXlamda[1]
selfmiuliantaXlamda <- selfparliantaXlamda[2]
selflamdamida<- selfliantaXlamdafunction(selfparliantaXlamda,X,R1)
selfOmigagXlamdamatrix <- matrix(data=c(1/(n*selflamdaliantaXlamda^2)+(1/n)*(JX+a+1)*(1/selflamdaliantaXlamda^2),
                                        -(1/n)*(-sum(1/(X-selfmiuliantaXlamda)^2)-selflamdaliantaXlamda*(2*sum(R1[1:JX]+1)+2*RXJ)),
                                        -(1/n)*(-sum(1/(X-selfmiuliantaXlamda)^2)-selflamdaliantaXlamda*(2*sum(R1[1:JX]+1)+2*RXJ)),
                                        (1/n)*(sum((R1[1:JX]+1)*(X-selfmiuliantaXlamda)^2)+RXJ*(TimX-selfmiuliantaXlamda)^2)),
                                 2,2)
selfOmigaGXlamdaValue <- det(solve(selfOmigagXlamdamatrix))
selfliantaGXlamdaValue <- -selflamdamida

selfparliantaXmiu<- optim(c(1,0.1),method='L-BFGS-B',selfliantaXmiufunction,lower=c(0.01,0.01),
                          hessian=T,X=X,R=R1)$par  #第三种情况1 0.1  
selflamdaliantaXmiu <- selfparliantaXmiu[1]
selfmiuliantaXmiu <- selfparliantaXmiu[2]
selfmiumida<- selfliantaXmiufunction(selfparliantaXmiu,X,R1)
selfOmigagXmiumatrix <- matrix(data=c((1/n)*(JX+a+1)*(1/selflamdaliantaXmiu^2),
                                      -(1/n)*(-sum(1/(X-selfmiuliantaXmiu)^2)-selflamdaliantaXmiu*(2*sum(R1[1:JX]+1)+2*RXJ)),
                                      -(1/n)*(-sum(1/(X-selfmiuliantaXmiu)^2)-selflamdaliantaXmiu*(2*sum(R1[1:JX]+1)+2*RXJ)),
                                      1/(n*selfmiuliantaXmiu^2)+(1/n)*(sum((R1[1:JX]+1)*(X-selfmiuliantaXmiu)^2)+RXJ*(TimX-selfmiuliantaXmiu)^2)),
                               2,2)
selfOmigaGXmiuValue <- det(solve(selfOmigagXmiumatrix))
selfliantaGXmiuValue <- -selfmiumida  

selftklamda<-tkvalue(selfOmigaGXlamdaValue,OmigaValue,liantaValue,selfliantaGXlamdaValue)
selftkmiu<-tkvalue(selfOmigaGXmiuValue,OmigaValue,liantaValue,selfliantaGXmiuValue)

library(xlsx)
ff <- data.frame(MLElamda=NRvalue$root[1],MLEmiu=NRvalue$root[2],
                 EMlamda=EMvalue$para[2],EMmiu=EMvalue$para[1],
                 selfbayeslamda=selfbayeslamda,selfbayesmiu=selfbayesmiu,
                 bayeslamdahpdL=bayeslamdahpd[1],bayeslamdahpdU=bayeslamdahpd[2],
                 bayesmiuhpdL=bayesmiuhpd[1],bayesmiuhpdU=bayesmiuhpd[2],
                 bayeslamdaCIL=bayeslamdaCI[1],bayeslamdaCIU=bayeslamdaCI[2],
                 bayesmiuCIL=bayesmiuCI[1],bayesmiuCIU=bayesmiuCI[2],
                 selflindlylamda=selflindlylamda,selflindlymiu=selflindlymiu,
                 selftklamda=selftklamda,selftkmiu=selftkmiu,
                 row.names = '估计值')

path <- file.path("E://progress result1119","case2realdataEstimateIP.xls")
write.xlsx(ff,file = path)
