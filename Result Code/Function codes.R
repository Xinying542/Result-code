#The functions utilized for running these codes
####Generate ordinary progressively Type-II censored data####
Progressive2usual<-function(m,R){
  Y <- numeric(m)
  W <- runif(m,0,1)
  U <- numeric(m)
  V <- numeric(m)
  O <- numeric(m)
  for(i in 1:m){
    V[i]<-W[i]^(1/(i+sum(R[(m-i+1):m])))
  }
  for (i in 1:m){
    U[i]<-1-prod(V[m:(m-i+1)])
  }
  O<-sort(U) 
  for (i in 1:m){
    Y[i]<-miu+sqrt(-(1/lamda)*log(1-O[i])) }
  return(Y)
}
####Generate generalized progressively hybrid censored data####
hybriddata <- function(Tim,m,c,Y,R){
  if(Tim < Y[c] && Y[c] <Y[m]){
    X=array(dim=c)
    X=Y[1:c]
    type=1
    JX=c
    RXJ=sum(R[(c+1):m]+1)
    lastR=R[1:JX]
    lastR[JX]=R[JX]+RXJ
    TimX=X[c]
    result=list(X=X,type=type,JX=JX,RXJ=RXJ,TimX=TimX,lastR=lastR)
  }
  if(Tim < Y[m] && Y[c] <Tim){
    for (i in 1:m-1){
      if(Y[i]<Tim && Tim<Y[i+1]){
        J<- i
      }
    }
    X=array(dim=J)
    X=Y[1:J]
    type=2
    JX=J
    RXJ=n-sum(R[1:J])-J
    lastR=R[1:JX]
    lastR[JX]=R[JX]+RXJ
    TimX=Tim
    result=list(X=X,type=type,JX=JX,RXJ=RXJ,TimX=TimX,lastR=lastR)
  }
  if(Y[c]< Y[m] && Y[m] <Tim){
    X=array(dim=m)
    X=Y[1:m]
    type=3
    JX=m
    RXJ=0
    lastR=R[1:JX]
    lastR[JX]=R[JX]+RXJ
    TimX=X[m]
    result=list(X=X,type=type,JX=JX,RXJ=RXJ,TimX=TimX,lastR=lastR)
  }
  return(result)
}
####Newton-Raphson algorithm####
NewtonsX = function(fun, data, theta,R,ep = 1e-5, it_max = 100){
  norm1 <- 0.1
  k=1
  while (norm1 > ep){
    theta1 = theta
    obj = fun(data,theta,R)
    theta = theta - solve(obj$Hes,obj$f)
    norm1 = sqrt((theta - theta1) %*% (theta - theta1))
    k=k+1
    # cat(k,'\n')
  }
  list(root = theta)
}
funs = function(data,theta,R){
  X=data$X
  JX=data$JX
  RXJ=data$RXJ
  TimX=data$TimX
  lamda=theta[1]
  miu=theta[2]
  gmiu <- sum((R[1:JX]+1)*(X-miu)^2)+RXJ*(TimX-miu)^2
  f=c(JX/lamda-gmiu,
      sum(-1/(X-miu))+2*lamda*(sum((X-miu))+sum(R[1:JX]*(X-miu))+RXJ*(TimX-miu)))
  ax=c(-JX/(lamda^2),
      2*sum(X-miu)+2*sum(R[1:JX]*(X-miu))+2*RXJ*(TimX-miu),
      2*sum(X-miu)+2*sum(R[1:JX]*(X-miu))+2*RXJ*(TimX-miu),
      -sum(1/((X-miu)^2))-2*lamda*JX-2*lamda*sum(R[1:JX])-2*RXJ*lamda)
  Hes=matrix(ax,2,2,byrow=T)
  list(f = f, Hes= Hes)
}

####Expectation-Maximization algorithm####
Estep <- function(miu,lamda){
  lamda=lamda
  miu=miu
  Z=matrix(nrow=JX,ncol = (n-m))
  ZTX=array(dim=RXJ)
  f <- function(x){
    2*lamda*x*(x-miu)*exp(-lamda*(x-miu)^2)
  }
  for(i in 1:JX){
    klist=integrate(f,X[i],Inf)
    k=klist$value 
    Z[i,]=k/exp(-lamda*(X[i]-miu)^2)
  }
  if(RXJ!=0){
    for(j in 1:RXJ){
      k2list=integrate(f,TimX,Inf)
      k2=k2list$value
      ZTX[j]=k2/exp(-lamda*(TimX-miu)^2)
    }
  }else{
    ZTX=0
  }
  Zvalue =list(Z=Z,ZTX=ZTX)
  return(Zvalue)
}
hfun <-function(miu,R,Z,ZTX){
  ax1=0
  for (i in 1:JX) {
    if(R[i]!=0){
      for(j in 1:R[i]){
        ax1=ax1+(Z[i,j]-miu)^2
      }
    }
  }
  bx1=ax1+sum((X-miu)^2)+sum((ZTX[1:RXJ]-miu)^2)
  SumR=JX+sum(R[1:JX])+RXJ
  lamdakmiu =SumR/bx1
  a2=0
  for(i in 1:JX){
    if(R[i]!=0){
      for(j in 1:R[i]){
        a2=a2+1/(Z[i,j]-miu)
      }
    }
  }
  a3=0
  for(i in 1:JX){
    if(R[i]!=0){
      for(j in 1:R[i]){
        a3=a3+Z[i,j]
      }
    }
  }
  Amiu=-sum(1/(X-miu))-a2-sum(1/(ZTX[1:RXJ]-miu))+2*lamdakmiu*(sum(X)+a3+sum(ZTX[1:RXJ]))
  h=Amiu/(2*lamdakmiu*SumR)
  return(h)  
}
fixpoint <- function(hfun, miu0, tol=1e-1, max.iter=1e4,R,Z,ZTX){
  R=R
  miu.old <- miu0
  miu.new <- miu0
  k=1
  while(k <= max.iter){
    miu.new <- hfun(miu.old,R,Z,ZTX)
    if(abs(miu.new - miu.old) < tol){
      return(miu.new)
      break
    }
    miu.old <- miu.new
    k=k+1
  }
}
EM12 <- function(R,miu0,lamda0,max_it=1000){
  
  miu.oldx <- miu0
  miu.newx <- miu0
  lamda.oldx <- lamda0
  lamda.newx <- lamda0
  
  K=1
  while (K <= max_it) {
    Zvalue=Estep(miu.oldx,lamda.oldx)
    Z =Zvalue$Z
    ZTX=Zvalue$ZTX
    miu.newx=fixpoint(hfun,miu.oldx,tol=1e-1, max.iter=1e4,R=R,Z=Z,ZTX=ZTX)
    ax1=0
    for (j in 1:JX) {
      if(R[j]!=0){
        for(k1 in 1:R[j]){
          ax1=ax1+(Z[j,k1]-miu.newx)^2
        }
      }
    }
    bx1=ax1+sum((X-miu.newx)^2)+sum((ZTX[1:RXJ]-miu.newx)^2)
    SumR=JX+sum(R[1:JX])+RXJ
    lamda.newx=SumR/bx1
    if(abs(miu.newx - miu.oldx) < 1e-1){
      res <- list(para=c(miu.newx,lamda.newx))
      return(res)
      break
    }else{
      miu.oldx <- miu.newx
      lamda.oldx <- lamda.newx
      K=K+1
    }
  }
}

####The functions used to calculate MSE, asymptotic confidence interval, CP, and AW.####
CI <- function(theta,tau=0.05,B=1000){
  L=B*tau/2
  R =B*(1-tau/2)
  Ci =c(theta[L],theta[R])
  return(Ci)
}
deltaCI <- function(theta,vartheta,alpha=0.05){
  a=exp(qnorm(1-alpha/2)*sqrt(vartheta/theta^2))
  leftpoint <- theta/a
  rightpoint <- theta*a
  interval <- c(leftpoint,rightpoint)
  return(interval)
}
Lengths <- function(MLEthetaInterval){
  mean(MLEthetaInterval[,2]-MLEthetaInterval[,1])
}
coverpercenge<-function(MLEthetaInterval,theta){
  num=0
  number=length(MLEthetaInterval)/2
  for(k in 1:number){
    a=MLEthetaInterval[k,1]
    b=MLEthetaInterval[k,2]
    if(theta >= a && theta<=b){
      num=num+1
    }
  }
  return(num/number)
}
MSE <- function(MLEtheta,theta,Nsim=10000){
  sum((MLEtheta-theta)^2)/Nsim
}

####The functions are used for MCMC algorithm. #####
Paimiu <- function(miu,A=2.056781e-22){
  Gmiu <- sum((R1[1:JX]+1)*(X-miu)^2)+RXJ*(TimX-miu)^2
  value = ((1/A)*prod(X-miu))/((b+Gmiu)^(JX+a))
  return(value)
}
logPaimiu <- function(x){
  Gmiu <- sum((R1[1:JX]+1)*(X-x)^2)+RXJ*(TimX-x)^2
  value = sum(log(X-x))-(JX+a)*log(b+Gmiu)
  return(value)
}
dlogPaimiu <- function(x){
  Gmiu <- sum((R1[1:JX]+1)*(X-x)^2)+RXJ*(TimX-x)^2
  Gmiup<- -2*sum((R1[1:JX]+1)*(X-x))-2*RXJ*(TimX-x)
  value <- sum(-1/(X-x))-(JX+a)*Gmiup/(b+Gmiu)
  return(value)
}

####The functions required for Lindley's approximation.####
lindlyvalue <- function(omiga1lag,omiga2lag,Taumatrix,Lthreematrix,Kmatrix){
  Amatrix <- matrix(ncol=2,nrow=2)
  B<- matrix(ncol=2,nrow=2)
  C<- matrix(ncol=2,nrow=2)
  A=0
  for (i in 1:2) {
    for (j in 1:2) {
      A=A+omiga2lag[i,j]*Taumatrix[i,j]
      Amatrix[i,j]=omiga1lag[i]*Taumatrix[i,i]+omiga1lag[j]*Taumatrix[j,i]
      B[i,j] =(omiga1lag[i]*Taumatrix[i,i]+omiga1lag[j]*Taumatrix[i,j])*Taumatrix[i,i]
      C[i,j]=3*omiga1lag[i]*Taumatrix[i,i]*Taumatrix[i,j]+omiga1lag[j]*(Taumatrix[i,i]*Taumatrix[j,j]+2*(Taumatrix[i,j])^2)
    }
  }
  value = 0.5*(A+Lthreematrix[1,1]*B[2,1]+Lthreematrix[2,2]*B[1,2]+Lthreematrix[1,2]*C[2,1]+Lthreematrix[2,1]*C[1,2])+
    Kmatrix[1]*Amatrix[1,2]+Kmatrix[2]*Amatrix[2,1]
  return(value)
}

####The functions required for TK approximation.####
liantafunction <- function(theta,X,R){
  lamda=theta[1]
  miu=theta[2]
  gmiu <- sum((R[1:JX]+1)*(X-miu)^2)+RXJ*(TimX-miu)^2
  lianta <- (1/n)*((JX+a-1)*log(lamda)+sum(log(X-miu))-lamda*(b+gmiu))
  return(-lianta)
}


selfliantaXlamdafunction <- function(theta,X,R){
  lamda=theta[1]
  miu=theta[2]
  gmiu <- sum((R[1:JX]+1)*(X-miu)^2)+RXJ*(TimX-miu)^2
  lianta <- (1/n)*log(lamda)+(1/n)*((JX+a-1)*log(lamda)+sum(log(X-miu))-lamda*(b+gmiu))
  return(-lianta)
}


selfliantaXmiufunction <- function(theta,X,R){
  lamda=theta[1]
  miu=theta[2]
  gmiu <- sum((R[1:JX]+1)*(X-miu)^2)+RXJ*(TimX-miu)^2
  lianta <- (1/n)*log(miu)+(1/n)*((JX+a-1)*log(lamda)+sum(log(X-miu))-lamda*(b+gmiu))
  return(-lianta)
}

tkvalue <- function(omigagx,omiga,lianta,liantagx){
  value <- sqrt(abs(omigagx/omiga))*exp(n*(liantagx-lianta))
  return(value)
}


