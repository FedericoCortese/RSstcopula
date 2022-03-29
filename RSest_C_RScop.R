# This script contains the functions for the estimation and simulation of a 
# Regime Switching Student t Copula (RSStC) model. 
# The functions needed to perform local and global decoding are also given.

# The following packages are required to load the C++ functions that implement
# the Expectation-Maximization (EM) algorithm.
library(Rcpp)
library(RcppArmadillo)
library(roptim)

# estimation --------------------------------------------------------------

# The file "RSest_C.cpp" must be saved on the working directory and loaded 
# via the following code. 
sourceCpp("RSest_C.cpp")

Est_comp_C=function(U, reg, maxiter = 1000,eps=1e-08, 
                    ninit = 1,h=0){
  
  # Est_comp_C estimates a RSStC model. 
  
  # ARGUMENTS:
  # U is the matrix of uniform pseudo-observations of dimension Txd where T is the number of time observations and d is the number of variables 
  # reg is the number of regimes of the RSStC model
  # maxiter and eps define the maximum number of iterations and the tolerance level for the convergence of the EM algorithm
  # ninit sets the number of tuning iterations of the algorithm
  # h is the weight put on the main diagonal of the initial estimate for the transition matrix

  # VALUE
  # R is an array of dimension dxdxreg, consisting of the estimated matrices of dependence parameters, each of dimension dxd
  # nu is the vector of the estimated number of degrees-of-freedom of the Student-t copula, it is of length equal to reg
  # Q is the estimated matrix of transition probabilities, it is of dimension regxreg
  # init is the estimated vector of initial probabilities, it is of the same length of nu
  # elaps_time is the time required for estimation of the RSStC model's parameters
  # iters is the number of iterations performed by the EM algorithm
  # w is the matrix of estimated posterior probabilities. It has dimension Txreg and can be used to perform local decoding
  # llk is the maximum log-likelihood at convergence of the EM algorithm
  # AIC is the value of the Akaike Information Criterion 
  # BIC is the value of the Bayesian Information Criterion 
  
  # DETAILS:
  # Regimes are ordered with respect to the determinant of the estimated dependence matrices
  # i.e., regimes labelled with smaller numbers are those with stronger dependence among variables
  
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
    #z=z,
    llk=llk,
    AIC=AIC,
    BIC=BIC
  ))
  
}

# The dataset "simest5" is a simulated dataset consisting of 500 observations and 5 variables
x=read.table("simdata5.txt")

# It is better to consider objects of class "matrix"
x=as.matrix(x)

# Est_comp_C requires a matrix of uniform pseudo-marginals as input. One possibility is to consider
# the normalized ranks of the observations, computed as follows
Ur=apply(x,2,rank)/(dim(x)[1]+1)

# The following code estimates a RSStC model with 3 states considering the previously simulated dataset
# We are adopting a semi-parametric approach given that no parametric assumptions are made on the marginals
est=Est_comp_C(Ur,3)

# The estimated matrices of dependence parameters can be extracted from the "est" object via the following
est$R

# In a similar way we can extract the estimated number of degrees-of-freedom
est$nu

# or the estimated matrix of transition probabilities
est$Q


# decoding ----------------------------------------------------------------

# Local decoding can be performed considering the estimated posterior probabilities given in 
# the output "w" of the object computed via the function "Est_comp_C"

loc_dec=apply(est$w,1,which.max)
plot(loc_dec,type='l',col="red",ylab="Decoded State",xlab="Time",main="Local Decoding")

# The following package is required to compute the Student-t copula density
library(copula)

RScop.viterbi <- function(U, nc){
  # This function performs the Viterbi algorithm for the global decoding
  
  # ARGUMENTS:
  # U is the matrix of uniform pseudo-observations, it is of dimension Txd
  # nc is the output object of the function "Est_comp_C"
  
  # VALUES:
  # iv is a vector of integers consisting of the decoded states sequence
  
  init = nc$init
  Q = nc$Q
  R=nc$R
  nu=nc$nu
  
  r=dim(Q)[2]
  d=dim(U)[2] 
  n=dim(U)[1]
  dc = matrix(0, n, r)
  
  for(i in 1:r){
    cop=tCopula(par=P2p(R[,,i]),dim=d,dispstr = "un",df=nu[i])
    dc[,i]=dCopula(U,cop)
  }
  
  dc=t(dc)
  
  xi <- matrix(0, n, r)
  foo <- init * dc[, 1]
  xi[1, ] <- foo / sum(foo)
  for (i in 2:n) {
    foo <-
      apply(xi[i - 1, ] * Q, 2, max) * dc[, i]
    xi[i, ] <- foo / sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (i in (n - 1):1) {
    iv[i] <- which.max(Q[, iv[i + 1]] * xi[i, ])
  }
  return(iv)
}

glob_dec=RScop.viterbi(Ur,est)
plot(glob_dec,type='l',col="blue",xlab="Time",ylab="Decoded State",main="Global Decoding")

# simulation --------------------------------------------------------------

copSim=function(n,d=5,R, nu, Q, init,seed){
  # This function simulates a d-variate dataset of uniform pseudo-observations
  # distributed as a RSStC model
  
  # ARGUMENTS:
  # n, number of time observations
  # d, number of variables
  # R, matrices of dependence parameters dxd, which is an array of dimension dxdxreg. See the example below for details on how to define R
  # nu, vector of the number of degrees-of-freedom, with length equal to the number of desired hidden states
  # Q, matrix of transition probabilities of dimension regxreg
  # init, vector of initial probabilities, it is of the same length of nu
  # seed, integer argument to set.seed for reproducible results
  
  # VALUES:
  # SimData, a matrix of dimension nxd consisting of the simulated uniform pseudo-marginals
  
  #markov chain simulation
  reg = dim(Q)[1]
  x <- numeric(n)
  set.seed(seed)
  x[1] <- sample(1:reg, 1, prob = init)
  for(i in 2:n){
    x[i] <- sample(1:reg, 1, prob = Q[x[i - 1], ])
  }
  
  Sim = matrix(0, n, d * reg)
  SimData = matrix(0, n, d)
  
  #pseudo-observations simulation
  for (k in 1:reg) {
    u = rCopula(n, copula::tCopula(param=P2p(R[,,k]), dim = d,df=nu[k],dispstr = "un"))
    Sim[, (d * k - d + 1):(d * k)] = u
  }
  
  for (i in 1:n) {
    k = x[i]
    SimData[i, ] = Sim[i, (d * k - d + 1):(d * k)]
  }
  return(SimData)
  
}

# The following is an example on how to simulate from a 3-states RSStC model with 3 variables

#R1 is the dependence matrix of the first states
R1=matrix(c(1,.9,.7,
            .9, 1,.75,
            .7, .75, 1),3,3)

#R2 is the dependence matrix of the second states
R2=matrix(c(1,.5,.3,
            .5, 1,-.4,
            .3, -.4, 1),3,3)

#R3 is the dependence matrix of the third states
R3=matrix(c(1,.1,.15,
            .1,1,-.1,
            .15,-.1,1),3,3)

#We collect R1 R2 and R3 in the array Rtrue
Rtrue=array(0,dim=c(3,3,3))
Rtrue[,,1]=R1;
Rtrue[,,2]=R2;
Rtrue[,,3]=R3;

nutrue=c(3,6,10)

Qtrue=matrix(c(.7,.2,.1,
               .4,.5,.1,
               .1,.1,.8),3,3,byrow=T)

inittrue=rep(1/3,3)


Us=copSim(1000,d=3,Rtrue, nutrue, Qtrue, inittrue,seed=1)




