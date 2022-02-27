library(Rcpp)
library(RcppArmadillo)
library(roptim)


sourceCpp("RSest_trede_C.cpp")

Est_comp_C=function(U, reg, maxiter = 1000,eps=1e-08, 
                    ninit = 1,h=0){

  n = dim(U)[1]
  d = dim(U)[2]
  r = reg * d
  
  R0=array(0,dim=c(d,d,reg))
  #Partitioning of y into "reg" parts, for each part we evaluate corr matrix
  n0 = floor((n / reg))
  ind0 = 1:n0
  for (j in 1:reg) {
    ind = (j - 1) * n0 + ind0
    x = U[ind,]
    R0[, , j] = cor(x)
  }
  nu0=rep(4,reg)
  #Q and init initialization (h is given as input)
  Q0=matrix(1/(h+reg),reg,reg)
  diag(Q0)=rep((h+1)/(h+reg),reg) 
  init0=rep(1,reg)/reg
  
  start_time <- Sys.time()
  est=RSest(U,reg,R0,nu0,Q0,init0,eps=eps, maxiter=maxiter,ninit=ninit)
  end_time = Sys.time()
  elaps_time=end_time-start_time
  R=est$R
  nu=est$nu
  Q = est$Q
  init=est$init
  w=est$w
  z=est$z
  llk=est$llk
  iters=est$iter
  
  ind=rank(apply(R,3,det))
  temp_Q=matrix(0,reg,reg)
  temp_i=rep(0,reg)
  temp_R=array(0,dim=c(d,d,reg))
  temp_w=matrix(0,n,reg)
  for(i in 1:reg){
    temp_i[ind[i]]=init[i]
    temp_R[,,ind[i]]=R[,,i]
    temp_w[,ind[i]]=w[,i]
    for(j in 1:reg){
      temp_Q[ind[i],ind[j]]=Q[i,j]
    }
  }
  
  Q=temp_Q
  nu=sort(nu)
  init=temp_i
  R=temp_R
  w=temp_w
  
  #k=length(est$nu)
  K=(reg-1)+reg*(reg-1)+reg+d*(d-1)/2
  #llk=est$llk
  AIC=2*K-2*llk
  BIC=K*log(n)-2*llk
  
  return(list(
    R=R,
    nu=nu,
    Q = Q,
    init=init,
    elaps_time = elaps_time,
    iters=iters,
    w=w ,
    z=z,
    llk=llk,
    AIC=AIC,
    BIC=BIC
  ))
  
}


tailInd=function(est){
  k=length(est$nu)
  d=dim(est$R)[1]
  lambda=array(0,dim=c(d,d,k))
  for(i in 1:k){
    R=est$R[,,i]
    nu=est$nu[i]
    lambda[,,i]=2*dt(-(sqrt(1+nu)*sqrt(1-R))/sqrt(1+R),df=nu+1)
    diag(lambda[,,i])=1
  }
  return(lambda)
}

kendTau=function(est){
  k=length(est$nu)
  d=dim(est$R)[1]
  tau=array(0,dim=c(d,d,k))
  for(i in 1:k){
    R=est$R[,,i]
    #nu=est$nu[i]
    tau[,,i]=asin(R)*2/pi
    diag(tau[,,i])=1
  }
  return(tau)
}


# simulation --------------------------------------------------------------
copSim=function(n,d=5,R, nu, Q, init,seed){
  
  #markov chain simulation
  reg = dim(Q)[1]
  x <- numeric(n)
  set.seed(seed)
  x[1] <- sample(1:reg, 1, prob = init)
  for(i in 2:n){
    x[i] <- sample(1:reg, 1, prob = Q[x[i - 1], ])
  }
  # 
  Sim = matrix(0, n, d * reg)
  SimData = matrix(0, n, d)
  
  for (k in 1:reg) {
    u = copula::rCopula(n, copula::tCopula(param=P2p(R[,,k]), dim = d,df=nu[k],dispstr = "un"))
    Sim[, (d * k - d + 1):(d * k)] = u
  }
  
  for (i in 1:n) {
    k = x[i]
    SimData[i, ] = Sim[i, (d * k - d + 1):(d * k)]
  }
  return(SimData)
  
}

###
R1=matrix(c(1,.9,.7,
            .9, 1,.75,
            .7, .75, 1),3,3)

R2=matrix(c(1,.5,.3,
            .5, 1,-.4,
            .3, -.4, 1),3,3)

R3=matrix(c(1,.1,.15,
            .1,1,-.1,
            .15,-.1,1),3,3)

Rtrue=array(0,dim=c(3,3,3))
Rtrue[,,1]=R1;
Rtrue[,,2]=R2;
Rtrue[,,3]=R3;

nutrue=c(3,6,10)

Qtrue=matrix(c(.7,.2,.1,
               .4,.5,.1,
               .1,.1,.8),3,3,byrow=T)
inittrue=rep(1/3,3)

library(copula)

Us=copSim(1000,d=5,Rtrue, nutrue, Qtrue, inittrue,seed=1)
