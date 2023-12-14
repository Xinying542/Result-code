####the EM algorithm is used to solve the maximum likelihood estimation when T=1.4.
library(HDInterval)
library(xlsx)
library(ars)

lamda=0.8
miu=0.5

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

EMlamdamean <- matrix(ncol = 2,nrow=12)
EMmiumean <- matrix(ncol = 2,nrow=12)

Tim=Timt[1]  
for (kx in 1:12) {
  R1 = unlist((CS[SD[kx]]))#将列表变为数组
  n = nt[kx]
  m= mt[kx]
  c=ct[kx]
  Nsim3=1000
  emmiuhat <- array(dim = Nsim3)
  emlamdahat <- array(dim = Nsim3)
  set.seed(20)
  for(i in 1:Nsim3){
    tryCatch({
      data1 <- Progressive2usual(m,R1)
      Datah<-hybriddata(Tim,m,c,data1,R1)
      X=Datah$X
      JX =Datah$JX
      RXJ=Datah$RXJ
      TimX=Datah$TimX
      re <- EM12(R1,0.5,0.5)
      a <- (re$para)[1]
      b <- (re$para)[2]
      emmiuhat[i]<-a
      emlamdahat[i]<-b
    },error=function(e){
      we=1
    })
  }
  
  EMlamdamean[kx,1]<-mean(emlamdahat,na.rm = T)
  EMmiumean[kx,1]<-mean(emmiuhat,na.rm = T)
  
  nanum <- sum(is.na(emmiuhat))
  EMlamdamean[kx,2]<-MSE(emlamdahat[complete.cases(emlamdahat)],0.8,(Nsim3-nanum))
  EMmiumean[kx,2]<-MSE(emmiuhat[complete.cases(emmiuhat)],0.5,(Nsim3-nanum))
  
}


library(xlsx)
ff <- data.frame(
  EMlamdamean=EMlamdamean[,1],EMlamdaMSE=EMlamdamean[,2],
  EMmiumean=EMmiumean[,1],EMmiuMSE=EMmiumean[,2],
  row.names = paste('R',1:12,sep = ''))

path <- file.path("F://progress result","emlast14x.xls")
write.xlsx(ff,file = path)
