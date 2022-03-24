#In this file the R functions to perform the simulation studies presented in the paper are given.

# first scenario ----------------------------------------------------------
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

stime_trede <- function(seed) {
  Us=copSim(1500,5,Rtrue,nutrue,Qtrue,inittrue,seed)
  
  est = Est_comp_C(Us, reg=3, maxiter = 1000,eps=1e-08, 
                   ninit = 1,h=0)
  return(est)
}

#D=1000
D=100
start = Sys.time()
sim_est_25 <- parallel::mclapply(X=1:D, stime_trede, mc.cores = parallel::detectCores())
end = Sys.time()
elapsed = end - start


# second scenario ---------------------------------------------------------
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
D=100
start = Sys.time()
sim_est_35 <- parallel::mclapply(X=1:D, stime_trede, mc.cores = parallel::detectCores())
end = Sys.time()
elapsed = end - start


# third scenario ----------------------------------------------------------
##varying T
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

n_sim <- seq(500, 2000, by=500)
n_seed <- 1:100

simulazione <- expand.grid(n_seed=n_seed, n_sim=n_sim)


stime_trede <- function(seed,n_sim,reg,Rtrue,nutrue,Qtrue,inittrue) {
  Us=copSim(n_sim,d=5,Rtrue,nutrue,Qtrue,inittrue,seed)
  
  est = Est_comp_C(Us, reg=reg, maxiter = 1000,eps=1e-08, 
                   ninit = 1,h=0)
  return(est)
}


start2 = Sys.time()
sim_est_n2 <- parallel::mclapply(X=1:nrow(simulazione), 
                                 function(x) stime_trede(seed=simulazione[x,]$n_seed,
                                                         n_sim=simulazione[x,]$n_sim,reg=2,
                                                         Rtrue=Rtrue2,nutrue=nutrue2,
                                                         Qtrue=Qtrue2,inittrue = inittrue2) ,
                                 mc.cores = parallel::detectCores())
end2 = Sys.time()
elapsed2 = end2 - start2

start3= Sys.time()
sim_est_n3 <- parallel::mclapply(X=1:nrow(simulazione), 
                                 function(x) stime_trede(seed=simulazione[x,]$n_seed,
                                                         n_sim=simulazione[x,]$n_sim,reg=3,
                                                         Rtrue=Rtrue3,nutrue=nutrue3,
                                                         Qtrue=Qtrue3,inittrue = inittrue3) ,
                                 mc.cores = parallel::detectCores())
end3= Sys.time()
elapsed3 = end3 - start3

##varying r

##2 states
Qtrue=matrix(c(.8,.2,
               .2,.8),2,2,byrow=T)
nutrue=c(5,15)
inittrue=c(1/2,1/2)

n_seed <- 1:100
n_d=c(2,5,10)

simulazione <- expand.grid(n_seed=n_seed, n_d=n_d)


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
sim_est_d <- parallel::mclapply(X=1:nrow(simulazione), 
                                function(x) stime_trede(seed=simulazione[x,]$n_seed,
                                                        d=simulazione[x,]$n_d,
                                                        nutrue=nutrue,
                                                        Qtrue=Qtrue,inittrue = inittrue) ,
                                mc.cores = parallel::detectCores())
endd = Sys.time()
elapsedd = endd - startd

##3 states
inittrue2=c(1/2,1/2)

Qtrue3=matrix(c(.7,.2,.1,
                .2,.7,.1,
                .1,.2,.7),3,3,byrow=T)
nutrue3=c(3,7,15)

n_seed <- 1:100
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
  Rtrue[,,2]=R3
  Us=copSim(1000,d=d,Rtrue,nutrue,Qtrue,inittrue,seed)
  
  est = Est_comp_C(Us, reg=3, maxiter = 1000,eps=1e-08, 
                   ninit = 1,h=0)
  return(est)
}


startd = Sys.time()
sim_est_d <- parallel::mclapply(X=1:nrow(simulazione), 
                                function(x) stime_trede(seed=simulazione[x,]$n_seed,
                                                        d=simulazione[x,]$n_d,
                                                        nutrue=nutrue,
                                                        Qtrue=Qtrue,inittrue = inittrue) ,
                                mc.cores = parallel::detectCores())
endd = Sys.time()
elapsedd = endd - startd

