#Numerical simulation estimation results when T=1.4 under NIP.
library(HDInterval)
library(xlsx)
library(ars)
library(HDInterval)

lamda=0.8
miu=0.5

a=0
b=0

R1=rep(0,25)
R1[25]=5

R2=rep(0,25)
R2[25]=15

R3=rep(0,35)
R3[35]=5

R4=rep(0,40)
R4[40]=20

R5=rep(0,40)
R5[1]=20

R6=rep(0,40)
R6[1]=R6[40]=10


CS=list(R1=R1,R2=R1,R3=R2,R4=R2,R5=R3,R6=R3,R7=R4,R8=R4,R9=R5,R10=R5,R11=R6,R12=R6)
SD<- paste("R", 1:12, sep = "")

nt=c(30,30,40,40,40,40,60,60,60,60,60,60)
mt=c(25,25,25,25,35,35,40,40,40,40,40,40)
ct=c(15,20,15,20,25,30,30,35,30,35,30,35)
Timt=c(1.4,1.7)  

MLElamdamean<- matrix(ncol = 2,nrow=12)
MLEmiumean <- matrix(ncol = 2,nrow=12)
MLEintervallamda <- matrix(ncol = 2,nrow=12)
MLEintervalmiu<- matrix(ncol = 2,nrow=12)

selfbayeslamdamean <-matrix(ncol = 2,nrow=12)
selfbayesmiumean <-matrix(ncol = 2,nrow=12)

hpdbayesintervallamda <- matrix(ncol = 2,nrow=12)
hpdbayesintervalmiu<-matrix(ncol = 2,nrow=12)
bcibayesintervallamda <- matrix(ncol = 2,nrow=12)
bcibayesintervalmiu<-matrix(ncol = 2,nrow=12)

selflindlylamdamean <- matrix(ncol = 2,nrow=12)
selflindlymiumean <- matrix(ncol = 2,nrow=12)

selftklamdamean <- matrix(ncol = 2,nrow=12)
selftkmiumean <- matrix(ncol = 2,nrow=12)

Tim=Timt[1]
for (kx in 1:12) {
  R1 = unlist((CS[SD[kx]]))
  n = nt[kx]
  m= mt[kx]
  c=ct[kx]
  Nsim1=1000

  theta0 <- c(0.3,0.5)
  MLElamda <- array(dim=Nsim1)
  MLEmiu <- array(dim=Nsim1)
  
  Nsim2 <- 1000
  theta0 = c(0.3,0.5)
  bigmiuCI <- matrix(nrow=Nsim2,ncol=2)
  biglamdaCI <- matrix(nrow=Nsim2,ncol=2)
  
  Nsim5 <- 1000
  Nsim6 <- 1000
  selfbayesmiu <- array(dim=Nsim5)
  selfbayeslamda <- array(dim=Nsim5)
  
  bayeslamdahpd <- matrix(nrow = Nsim5,ncol=2)
  bayesmiuhpd <- matrix(nrow=Nsim5,ncol=2)
  bayeslamdabci<-matrix(nrow = Nsim5,ncol=2)
  bayesmiubci<-matrix(nrow = Nsim5,ncol=2)
  
  Nsim8 <- 1000
  
  selflindlymiu <- array(dim=Nsim8)
  selflindlylamda <- array(dim=Nsim8)
  
  Nsim9 <- 1000
  selftkmiu <- array(dim=Nsim9)
  selftklamda <- array(dim=Nsim9)
  
  set.seed(20)
  
  for (i in 1:Nsim1) {
    data1 <- Progressive2usual(m,R1)
    Datah<-hybriddata(Tim,m,c,data1,R1)
    X=Datah$X
    JX=Datah$JX
    RXJ=Datah$RXJ
    TimX=Datah$TimX
    
    tryCatch({
      MLE <- NewtonsX(funs,Datah,theta0,R1)
      MLEe <- MLE$root
      MLElamda[i] <- MLEe[1]
      MLEmiu[i] <- MLEe[2]
    },error=function(e){
      we=1
    })
    
    tryCatch({
      MLE <- NewtonsX(funs,Datah,theta0,R1)
      miuhat1<-(MLE$root)[2]
      lamdahat1 <- (MLE$root)[1]
      FisherI <- matrix(data=c(JX/(lamdahat1^2),
                               -(2*sum(X-miuhat1)+2*sum(R1[1:JX]*(X-miuhat1))+2*RXJ*(TimX-miuhat1)),
                               -(2*sum(X-miuhat1)+2*sum(R1[1:JX]*(X-miuhat1))+2*RXJ*(TimX-miuhat1)),
                               -(-sum(1/((X-miuhat1)^2))-2*lamdahat1*JX-2*lamdahat1*sum(R1[1:JX])-2*RXJ*lamdahat1)),
                        2,2)
      covMatrix<-solve(FisherI)
      varmiu <- abs(covMatrix[2,2])
      varlamda <- abs(covMatrix[1,1])
      bigmiuCI[i,]<-deltaCI(miuhat1,varmiu,alpha = 0.05)
      biglamdaCI[i,]<-deltaCI(lamdahat1,varlamda)
    },error=function(e){
      we=1
    })
    
    gibbsmiu <- array(dim=Nsim6)
    gibbslamda <- array(dim=Nsim6)
    for (j in 1:Nsim6) {
      
      tryCatch({
        gibbsmiu[j]<- ars(1,logPaimiu,dlogPaimiu,x=0.5,m=1,lb=TRUE,xlb=0.01,ub=TRUE,xub=X[1])
        Gmiu<- sum((R1[1:JX]+1)*(X-gibbsmiu[j])^2)+RXJ*(TimX-gibbsmiu[j])^2
        gibbslamda[j] <- rgamma(1,shape = (JX+a),rate = (b+Gmiu))
      },error=function(e){
        we=1
      })
    }
    selfbayesmiu[i] <- mean(gibbsmiu,na.rm = T)
    selfbayeslamda[i] <- mean(gibbslamda,na.rm=T)
    tryCatch({
      bayeslamdahpd[i,]<-hdi(gibbslamda[complete.cases(gibbslamda)],credMass = 0.95)
      bayesmiuhpd[i,]<-hdi(gibbsmiu[complete.cases(gibbsmiu)],credMass = 0.95)
      bayeslamdabci[i,]<-CI(sort(gibbslamda[complete.cases(gibbslamda)]),tau=0.05,B=length(gibbslamda[complete.cases(gibbslamda)]))
      bayesmiubci[i,]<-CI(sort(gibbsmiu[complete.cases(gibbsmiu)]),tau=0.05,B=length(gibbsmiu[complete.cases(gibbsmiu)]))
    },error=function(e){
      we=1
    })
    
    theta0 = c(0.3,0.5)
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
    selflindlylamda[i] <- lamdahat1+value1
    
    omiga1lagmiu1 <- matrix(data=c(0,1),1,2)
    omiga2lagmiu1 <- matrix(data=c(0,0,0,0),2,2,byrow = T)
    value2<- lindlyvalue(omiga1lagmiu1,omiga2lagmiu1,Taumatrix,Lthreematrix,Kmatrix)
    selflindlymiu[i] <- miuhat1+value2
    
    tryCatch({
      parlianta<- optim(c(0.3,0.5),method='L-BFGS-B',liantafunction,lower=c(0.01,0.01),
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
      
      selfparliantaXlamda<- optim(c(0.3,0.5),method='L-BFGS-B',selfliantaXlamdafunction,lower=c(0.01,0.01),
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
     
      selfparliantaXmiu<- optim(c(0.3,0.5),method='L-BFGS-B',selfliantaXmiufunction,lower=c(0.01,0.01),
                                hessian=T,X=X,R=R1)$par
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
      
      selftklamda[i] <- tkvalue(selfOmigaGXlamdaValue,OmigaValue,liantaValue,selfliantaGXlamdaValue)
      selftkmiu[i] <- tkvalue(selfOmigaGXmiuValue,OmigaValue,liantaValue,selfliantaGXmiuValue)
      },error=function(e){
        we=1
      })
  }
  
  MLElamdamean[kx,1]<-mean(MLElamda,na.rm =T)
  MLEmiumean[kx,1]<-mean(MLEmiu,na.rm =T)
  MLElamdamean[kx,2]<-MSE(MLElamda,0.8,Nsim=length(MLElamda))
  MLEmiumean[kx,2]<-MSE(MLEmiu,0.5,Nsim=length(MLEmiu))
  
  MLEintervallamda[kx,1]<-Lengths(biglamdaCI[complete.cases(biglamdaCI),])#1.89
  MLEintervalmiu[kx,1]<-Lengths(bigmiuCI[complete.cases(bigmiuCI),])#0.4642
  MLEintervallamda[kx,2]<-coverpercenge(biglamdaCI[complete.cases(biglamdaCI),],0.8)#0.8842
  MLEintervalmiu[kx,2]<-coverpercenge(bigmiuCI[complete.cases(bigmiuCI),],0.5)#0.7704
  
  selfbayeslamdamean[kx,1]<-mean(selfbayeslamda,na.rm = T)
  selfbayesmiumean[kx,1]<-mean(selfbayesmiu,na.rm = T)
  selfbayeslamdamean[kx,2]<-MSE(selfbayeslamda[complete.cases(selfbayeslamda)],0.8,Nsim = length(selfbayeslamda[complete.cases(selfbayeslamda)]))
  selfbayesmiumean[kx,2]<-MSE(selfbayesmiu[complete.cases(selfbayesmiu)],0.5,Nsim = length(selfbayesmiu[complete.cases(selfbayesmiu)]))
  
  hpdbayesintervallamda[kx,1]<-Lengths(bayeslamdahpd[complete.cases(bayeslamdahpd),])
  hpdbayesintervalmiu[kx,1]<-Lengths(bayesmiuhpd[complete.cases(bayesmiuhpd),])
  hpdbayesintervallamda[kx,2]<-coverpercenge(bayeslamdahpd[complete.cases(bayeslamdahpd),],0.8)
  hpdbayesintervalmiu[kx,2]<-coverpercenge(bayesmiuhpd[complete.cases(bayesmiuhpd),],0.5)
  
  bcibayesintervallamda[kx,1]<-Lengths(bayeslamdabci[complete.cases(bayeslamdabci),])
  bcibayesintervalmiu[kx,1]<-Lengths(bayesmiubci[complete.cases(bayesmiubci),])
  bcibayesintervallamda[kx,2]<-coverpercenge(bayeslamdabci[complete.cases(bayeslamdabci),],0.8)
  bcibayesintervalmiu[kx,2]<-coverpercenge(bayesmiubci[complete.cases(bayesmiubci),],0.5)
  
  selflindlylamdamean[kx,1]<-mean(selflindlylamda,na.rm = T)
  selflindlymiumean[kx,1]<-mean(selflindlymiu,na.rm = T)
  selflindlylamdamean[kx,2]<-MSE(selflindlylamda[complete.cases(selflindlylamda)],0.8,Nsim = length(selflindlylamda[complete.cases(selflindlylamda)]))
  selflindlymiumean[kx,2] <- MSE(selflindlymiu[complete.cases(selflindlymiu)],0.5,Nsim = length(selflindlymiu[complete.cases(selflindlymiu)]))
  
  selftklamdamean[kx,1]<-mean(selftklamda,na.rm = T)
  selftkmiumean[kx,1]<-mean(selftkmiu,na.rm = T)
  selftklamdamean[kx,2]<-MSE(selftklamda[complete.cases(selftklamda)],0.8,Nsim = length(selftklamda[complete.cases(selftklamda)]))
  selftkmiumean[kx,2] <- MSE(selftkmiu[complete.cases(selftkmiu)],0.5,Nsim = length(selftkmiu[complete.cases(selftkmiu)]))
  
}

library(xlsx)
ff <- data.frame(MLElamdamean=MLElamdamean[,1],MLElamdaMSE=MLElamdamean[,2],
                 MLEmiumean=MLEmiumean[,1],MLEmiuMSE=MLEmiumean[,2],
                 
                 selfbayeslamdamean=selfbayeslamdamean[,1],selfbayeslamdaMSE=selfbayeslamdamean[,2],
                 selfbayesmiumean=selfbayesmiumean[,1],selfbayesmiuMSE=selfbayesmiumean[,2],

                 selflindlylamdamean=selflindlylamdamean[,1],selflindlylamdaMSE=selflindlylamdamean[,2],
                 selflindlymiumean=selflindlymiumean[,1],selflindlymiuMSE=selflindlymiumean[,2],

                 selftklamdamean=selftklamdamean[,1],selftklamdaMSE=selftklamdamean[,2],
                 selftkmiumean=selftkmiumean[,1],selftkmiuMSE=selftkmiumean[,2],

                 MLEintlengthamda=MLEintervallamda[,1],MLEintcovplamda=MLEintervallamda[,2],
                 MLEintlengthmiu=MLEintervalmiu[,1],MLEintcovpmiu=MLEintervalmiu[,2],
                 bcibayesintlengthlamda=bcibayesintervallamda[,1],bcibayesintcovplamda=bcibayesintervallamda[,2],
                 bcibayesintlengthmiu=bcibayesintervalmiu[,1],bcibayesintcovpmiu=bcibayesintervalmiu[,2],
                 hpdbayesintlengthlamda=hpdbayesintervallamda[,1],hpdbayesintcovplamda=hpdbayesintervallamda[,2],
                 hpdbayesintlengthmiu=hpdbayesintervalmiu[,1],hpdbayesintcovpmiu=hpdbayesintervalmiu[,2],
                 
                 row.names = paste('R',1:12,sep = ''))

path <- file.path("F://progress result","t14NIPlast.xls")
write.xlsx(ff,file = path)
