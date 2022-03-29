# This script contains the R functions to perform the simulation studies presented in the paper 

# We make use of the function "mclapply" from the R package "parallel" to perform parallel computing
# We recall that it can be used only on Linux and Mac computers because it relies on forking and hence is not available on Windows


# First scenario ----------------------------------------------------------
# The first scenario refers to a 2-states RSStC model: the first state is characterized by
# strong dependencies, while the second is a more tranquil market regime

# In the following, we set the 2-states RSStC parameters 
Qtrue2=matrix(c(.8,.2,
                .2,.8),2,2,byrow=T)
nutrue2=c(5,15)
R1=matrix(c(1,0.9,0.9,0.9,0.9,
            0.9,1,0.9,0.9,0.9,
            0.9,0.9,1,0.9,0.9,
            0.9,0.9,0.9,1,0.9,
            0.9,0.9,0.9,0.9,1),5,5,byrow = T)
R2=matrix(c(1.000,  0.2 , 0 , 0.1 ,0,  
            0.2,  1.000, 0 , 0.2,0.2,
            0 , 0 , 1.000  ,0.1  ,0.2,
            0.1 , 0.2  ,0.1, 1.000,  0,
            0 , 0.2 , 0.2 , 0 , 1.000),5,5,byrow=T)
Rtrue2=array(0,dim=c(5,5,2))
Rtrue2[,,1]=R1;
Rtrue2[,,2]=R2;
inittrue2=rep(.5,2)

# The following function is passed to "mclapply" which performs parallel computing
stime_trede <- function(seed) {
  Us=copSim(1500,5,Rtrue2,nutrue2,Qtrue2,inittrue2,seed)
  
  est = Est_comp_C(Us, reg=2, maxiter = 1000,eps=1e-08, 
                   ninit = 1,h=0)
  return(est)
}

#D is the number of simulated datasets
#D=1000
D=2
start = Sys.time()
sim_est_25 <- parallel::mclapply(X=1:D, stime_trede, mc.cores = parallel::detectCores()-1)
end = Sys.time()
elapsed = end - start

# Biases of the estimates are computed via the following code
# The mean absolute difference between the average of the simulated estimates and the true
# parameters is considered as a measure of the bias

Qsim = sapply(sim_est_25, function(x) x$Q)
Qmean=matrix(rowMeans(Qsim),2,2)
abs(Qmean-Qtrue2)

Rsim=sapply(sim_est_25, function(x) x$R)
Rmean=array(rowMeans(Rsim),dim=c(5,5,2))
abs(Rmean-Rtrue2)

initsim=sapply(sim_est_25, function(x) x$init)
initmean=c(rowMeans(initsim))
abs(initmean-inittrue2)

nusim=sapply(sim_est_25, function(x) x$nu)
numean=c(rowMeans(nusim))
abs(numean-nutrue2)


# Second scenario ---------------------------------------------------------

# The second scenario refers to a 3-states RSStC model: the first state is characterized by
# strong dependencies, the second by moderate dependence, the third is a more tranquil market regime

# In the following, we set the 3-states RSStC parameters 
R1=matrix(c(1,.9,.7,.8,.8,
            .9, 1,.75,.9,.8,
            .7, .75, 1,.7,.8,
            .8,.9,.7,1,.8,
            .8,.8,.8,.8,1),5,5)

R2=matrix(c(1,.5,.3,.5,.4,
            .5, 1,.4,.4,.5,
            .3, .4, 1,.4,.5,
            .5,.4,.4,1,.3,
            .4,.5,.5,.3,1),5,5)

R3=matrix(c(1,.1,.15,.05,.05,
            .1,1,-.1,.1,-0.05,
            .15,-.1,1,.05,.1,
            .05,.1,.05,1,-.01,
            .05,-.05,.1,-.01,1),5,5)

Rtrue=array(0,dim=c(5,5,3))
Rtrue[,,1]=R1;
Rtrue[,,2]=R2;
Rtrue[,,3]=R3;

nutrue=c(3,6,10)

Qtrue=matrix(c(.7,.2,.1,
               .3,.6,.1,
               .1,.1,.8),3,3,byrow=T)
inittrue=rep(1/3,3)

stime_trede <- function(seed) {
  Us=copSim(1500,5,Rtrue,nutrue,Qtrue,inittrue,seed)
  
  est = Est_comp_C(Us, reg=3, maxiter = 1000,eps=1e-08, 
                   ninit = 1,h=0)
  return(est)
}

#D=1000
D=2
start = Sys.time()
sim_est_35 <- parallel::mclapply(X=1:D, stime_trede, mc.cores = parallel::detectCores()-1)
end = Sys.time()
elapsed = end - start

# Biases of the estimates are computed via the following code
Qsim = sapply(sim_est_35, function(x) x$Q)
Qmean=matrix(rowMeans(Qsim),3,3)
abs(Qmean-Qtrue)

Rsim=sapply(sim_est_35, function(x) x$R)
Rmean=array(rowMeans(Rsim),dim=c(5,5,3))
abs(Rmean-Rtrue)

initsim=sapply(sim_est_35, function(x) x$init)
initmean=c(rowMeans(initsim))
abs(initmean-inittrue)

nusim=sapply(sim_est_35, function(x) x$nu)
numean=c(rowMeans(nusim))
abs(numean-nutrue)


# Varying the number of time observations -----------------------------------------

# In this third scenario we focus on the computational time demanded for the estimation of 2- and 3-states
# RSStC models for varying number of time observations and variables

# In this first step, we consider 2- and 3-states RSStC models varying the number 
# of time observations 

Qtrue2=matrix(c(.8,.2,
                .2,.8),2,2,byrow=T)
nutrue2=c(5,15)
R1=matrix(c(1,0.9,0.9,0.9,0.9,
            0.9,1,0.9,0.9,0.9,
            0.9,0.9,1,0.9,0.9,
            0.9,0.9,0.9,1,0.9,
            0.9,0.9,0.9,0.9,1),5,5,byrow = T)
R2=matrix(c(1.000,  0.2 , 0 , 0.1 ,0,  
            0.2,  1.000, 0 , 0.2,0.2,
            0 , 0 , 1.000  ,0.1  ,0.2,
            0.1 , 0.2  ,0.1, 1.000,  0,
            0 , 0.2 , 0.2 , 0 , 1.000),5,5,byrow=T)
Rtrue2=array(0,dim=c(5,5,2))
Rtrue2[,,1]=R1;
Rtrue2[,,2]=R2;

inittrue2=c(1/2,1/2)

Qtrue3=matrix(c(.7,.2,.1,
                .2,.7,.1,
                .1,.2,.7),3,3,byrow=T)
nutrue3=c(3,7,15)
R1=matrix(c(1,0.9,0.9,0.9,0.9,
            0.9,1,0.9,0.9,0.9,
            0.9,0.9,1,0.9,0.9,
            0.9,0.9,0.9,1,0.9,
            0.9,0.9,0.9,0.9,1),5,5,byrow = T)
R2=matrix(c(1.000, 0.4, 0.5, 0.5, 0.4,
            0.4, 1.000, 0.4, 0.5, 0.4,
            0.5,0.4 ,1.000, 0.5, 0.5,
            0.5, 0.5 ,0.5 ,1.000 ,0.4,
            0.4, 0.4, 0.5 ,0.4, 1.000),5,5,byrow=T)

R3=matrix(c(1.000,  0.2 , 0 , 0.1 ,0,  
            0.2,  1.000, 0 , 0.2,0.2,
            0 , 0 , 1.000  ,0.1  ,0.2,
            0.1 , 0.2  ,0.1, 1.000,  0,
            0 , 0.2 , 0.2 , 0 , 1.000),5,5,byrow=T)

Rtrue3=array(0,dim=c(5,5,3))
Rtrue3[,,1]=R1;
Rtrue3[,,2]=R2;
Rtrue3[,,3]=R3;

inittrue3=c(1/3,1/3,1/3)

# Four possible time lengths are considered: 500, 1000, 1500 and 2000
n_sim <- seq(500, 2000, by=500)
#D=1000
#n_seed <- 1:D
D=2
n_seed <- 1:D

simulazione <- expand.grid(n_seed=n_seed, n_sim=n_sim)

stime_trede <- function(seed,n_sim,reg,Rtrue,nutrue,Qtrue,inittrue) {
  Us=copSim(n_sim,d=5,Rtrue,nutrue,Qtrue,inittrue,seed)
  
  est = Est_comp_C(Us, reg=reg, maxiter = 1000,eps=1e-08, 
                   ninit = 1,h=0)
  return(est)
}

startd = Sys.time()
sim_est_T2 <- parallel::mclapply(X=1:nrow(simulazione),
                                function(x) stime_trede(seed=simulazione[x,]$n_seed,
                                                        n_sim=simulazione[x,]$n_sim,reg=2,
                                                        nutrue=nutrue2,Rtrue=Rtrue2,
                                                        Qtrue=Qtrue2,inittrue = inittrue2) ,
                                mc.cores = parallel::detectCores()-1)
endd = Sys.time()
elapsedd = endd - startd

startd = Sys.time()
sim_est_T3 <- parallel::mclapply(X=1:nrow(simulazione),
                                function(x) stime_trede(seed=simulazione[x,]$n_seed,
                                                        n_sim=simulazione[x,]$n_sim,reg=3,
                                                        nutrue=nutrue3,Rtrue=Rtrue3,
                                                        Qtrue=Qtrue3,inittrue = inittrue3) ,
                                mc.cores = parallel::detectCores()-1)
endd = Sys.time()
elapsedd = endd - startd

# Average computational times for the 2-states RSStC model are computed via the following
library(hms)
time_T2 = sapply(sim_est_T2, function(x) as_hms(x$elaps_time))
t2_500=mean(time_T2[1:D])
t2_1000=mean(time_T2[(D+1):(2*D)])
t2_1500=mean(time_T2[(2*D+1):(3*D)])
t2_2000=mean(time_T2[(3*D+1):(4*D)])

# Average computational times for the 3-states RSStC model are computed via the following
time_T3 = sapply(sim_est_T3, function(x) as_hms(x$elaps_time))
t3_500=mean(time_T3[1:D])
t3_1000=mean(time_T3[(D+1):(2*D)])
t3_1500=mean(time_T3[(2*D+1):(3*D)])
t3_2000=mean(time_T3[(3*D+1):(4*D)])


# Varying the number of variables -----------------------------------------

# In this second step, we consider 2- and 3-states RSStC models varying the number 
# of variables

# The followings are the parameters for simulating a 2-states RSStC model
Qtrue=matrix(c(.8,.2,
               .2,.8),2,2,byrow=T)
nutrue=c(5,15)
inittrue=c(1/2,1/2)

#D=1000
#n_seed <- 1:D
D=2
n_seed <- 1:D

# Three possible values are considered: 2,5 and 10
n_d=c(2,5,10)

simulazione <- expand.grid(n_seed=n_seed, n_d=n_d)

library(copula)
stime_trede <- function(seed,d,nutrue,Qtrue,inittrue) {
  R1=matrix(.9,d,d)
  diag(R1)=1
  R2=matrix(.1,d,d)
  diag(R2)=1
  Rtrue=array(0,dim=c(d,d,2))
  Rtrue[,,1]=R1
  Rtrue[,,2]=R2
  Us=copSim(1000,d=d,Rtrue,nutrue,Qtrue,inittrue,seed)
  
  est = Est_comp_C(Us, reg=2, maxiter = 1000,eps=1e-08, 
                   ninit = 1,h=0)
  return(est)
}


startd = Sys.time()
sim_est_d2 <- parallel::mclapply(X=1:nrow(simulazione), 
                                function(x) stime_trede(seed=simulazione[x,]$n_seed,
                                                        d=simulazione[x,]$n_d,
                                                        nutrue=nutrue,
                                                        Qtrue=Qtrue,inittrue = inittrue) ,
                                mc.cores = parallel::detectCores()-1)
endd = Sys.time()
elapsedd = endd - startd

# In the following, we compute the average computational time demanded for the estimation of the 
# 2-states RSStC model for vatying time length T
time_d2 = sapply(sim_est_d2, function(x) as_hms(x$elaps_time))
d2_2=mean(time_d2[1:D])
d2_5=mean(time_d2[(D+1):(2*D)])
d2_10=mean(time_d2[(2*D+1):(3*D)])

# The followings are the parameters for simulating a 3-states RSStC model
inittrue=rep(1/3,3)

Qtrue=matrix(c(.7,.2,.1,
                .2,.7,.1,
                .1,.2,.7),3,3,byrow=T)
nutrue=c(3,7,15)

#D=100
#n_seed <- 1:100
D=2
n_seed <- 1:D
n_d=c(2,5,10)

simulazione <- expand.grid(n_seed=n_seed, n_d=n_d)


stime_trede <- function(seed,d,nutrue,Qtrue,inittrue) {
  R1=matrix(.9,d,d)
  diag(R1)=1
  R2=matrix(.5,d,d)
  diag(R2)=1
  R3=matrix(.1,d,d)
  diag(R3)=1
  Rtrue=array(0,dim=c(d,d,3))
  Rtrue[,,1]=R1
  Rtrue[,,2]=R2
  Rtrue[,,3]=R3
  Us=copSim(1000,d=d,Rtrue,nutrue,Qtrue,inittrue,seed)
  
  est = Est_comp_C(Us, reg=3, maxiter = 1000,eps=1e-08, 
                   ninit = 1,h=0)
  return(est)
}


startd = Sys.time()
sim_est_d3 <- parallel::mclapply(X=1:nrow(simulazione), 
                                function(x) stime_trede(seed=simulazione[x,]$n_seed,
                                                        d=simulazione[x,]$n_d,
                                                        nutrue=nutrue,
                                                        Qtrue=Qtrue,inittrue = inittrue) ,
                                mc.cores = parallel::detectCores()-1)
endd = Sys.time()
elapsedd = endd - startd

# In the following, we compute the average computational time demanded for the estimation of the 
# 3-states RSStC model for vatying time length T
time_d3 = sapply(sim_est_d3, function(x) as_hms(x$elaps_time))
d3_2=mean(time_d3[1:D])
d3_5=mean(time_d3[(D+1):(2*D)])
d3_10=mean(time_d3[(2*D+1):(3*D)])

