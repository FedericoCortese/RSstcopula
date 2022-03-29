# This script contains the R functions needed to perform the block-bootstraps presented 
# in the paper.

# We make use of the function "mclapply" from the R package "parallel" to perform parallel computing
# We recall that it can be used only on Linux and Mac because it relies on forking and hence is not available on Windows

library(bootstrap)
library(boot)
library(parallel)
library(rugarch)
library(fGarch)
library(skewt)

x=read.table("simdata5.txt")
x=as.matrix(x)
d=dim(x)[2]

# semiparam__2-states -----------------------------------------------------

semiparametric_RScop_boot <- function(x) {
  
  # This function performs the bootstrap for the estimation of standard errors of RSStC model parameters
  # with 2 states, estimated adopting the semiparametric approach
  
  # ARGUMENTS:
  # x is the matrix of observations
  
  # VALUES:
  # cop.pars is the vector of RSStC copula parameters
  
  # marginals
  R=apply(x,2,rank)/(dim(x)[1]+1)
  
  # RS copula param
  est_rsc=Est_comp_C(R,2)
  Qvec=as.vector(est_rsc$Q)
  Rvec=as.vector(apply(est_rsc$R,3,function(x)x[upper.tri(x)]))
  cop.pars=c(est_rsc$init,Qvec,Rvec,est_rsc$nu)
  return(cop.pars)
}

# The average block-length for the block-bootstrap is the equal to l1 as suggested in Politis and Roman (1994)
l1=round(dim(x)[1]^(2/3))
#nboot=1000
nboot=2

set.seed(1)
startsp2=Sys.time()
param_boot_2sp <- tsboot(x, semiparametric_RScop_boot, R = nboot, l = l1, sim = "geom",
                         parallel = "multicore", ncpus = detectCores()-1, n.sim=dim(x)[1])

endsp2=Sys.time()
endsp2-startsp2

sp2=param_boot_2sp$t
SE_sp2=apply(sp2,2,sd)
reg=2

# Standard errors:

# Initial probabilities
round(SE_sp2[1:reg],3)

# Transition probabilities
round(SE_sp2[(reg+1):(reg*reg+reg)],3)

# Dependence matrices
round(SE_sp2[(reg*reg+reg+1):(length(SE_sp2)-reg)],3)

# Number of degrees-of-freedom
round(tail(SE_sp2,reg),3)

# semiparam__3-states -----------------------------------------------------

semiparametric_RScop_boot <- function(x) {
  # This function performs the bootstrap for the estimation of standard errors of RSStC model parameters
  # with 3 states, estimated adopting the semiparametric approach
  
  # ARGUMENTS:
  # x is the matrix of observations
  
  # VALUES:
  # cop.pars is the vector of RSStC copula parameters
  
  #marginals
  R=apply(x,2,rank)/(dim(x)[1]+1)
  
  #RS copula param
  est_rsc=Est_comp_C(R,3)
  Qvec=as.vector(est_rsc$Q)
  Rvec=as.vector(apply(est_rsc$R,3,function(x)x[upper.tri(x)]))
  #Rvec=as.vector(est_rsc$R)
  cop.pars=c(est_rsc$init,Qvec,Rvec,est_rsc$nu)
  return(cop.pars)
}

l1=round(dim(x)[1]^(2/3))
#nboot=1000
nboot=2

set.seed(1)
startsp3=Sys.time()
param_boot_3sp <- tsboot(x, semiparametric_RScop_boot, R = nboot, l = l1, sim = "geom",
                         parallel = "multicore", ncpus = detectCores()-1, n.sim=dim(x)[1])

endsp3=Sys.time()
endsp3-startsp3

sp3=param_boot_3sp$t
SE_sp3=apply(sp3,2,sd)
reg=3

# Standard errors:

# Initial probabilities
round(SE_sp3[1:reg],3)

# Transition probabilities
round(SE_sp3[(reg+1):(reg*reg+reg)],3)

# Dependence matrices
round(SE_sp3[(reg*reg+reg+1):(length(SE_sp3)-reg)],3)

# Number of degrees-of-freedom
round(tail(SE_sp3,reg),3)

# param__2-states ---------------------------------------------------------

parametric_RScop_boot <- function(x) {
  # This function performs the bootstrap for the estimation of standard errors of RSStC model parameters
  # with 2 states, estimated adopting the parametric approach and assuming an 
  # ARMA(1,1)-GARCH(1,1) with skewed generalized error distribution model for the marginals
  
  # ARGUMENTS:
  # x is the matrix of observations
  
  # VALUES:
  # cop.pars is the vector of RSStC copula parameters
  
  #marginal parameters
  btc = x[,1]
  eth=x[,2]
  xrp=x[,3]
  ltc=x[,4]
  bch=x[,5]

  model = "sGARCH"
  spec = ugarchspec(
    variance.model = list(model = model, garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(1, 1), include.mean = FALSE) ,
    distribution.model = "sged"
  )

  fit.btc = ugarchfit(spec, btc)
  fit.eth = ugarchfit(spec, eth)
  fit.ltc = ugarchfit(spec, ltc)
  fit.xrp = ugarchfit(spec, xrp)
  fit.bch = ugarchfit(spec, bch)

  marg.pars=c(fit.btc@fit$coef,
              fit.eth@fit$coef,
              fit.ltc@fit$coef,
              fit.xrp@fit$coef,
              fit.bch@fit$coef)

  nu.b = fit.btc@fit$coef["shape"]
  nu.e = fit.eth@fit$coef["shape"]
  nu.l = fit.ltc@fit$coef["shape"]
  nu.x = fit.xrp@fit$coef["shape"]
  nu.h = fit.bch@fit$coef["shape"]

  sk.b = fit.btc@fit$coef["skew"]
  sk.e = fit.eth@fit$coef["skew"]
  sk.l = fit.ltc@fit$coef["skew"]
  sk.x = fit.xrp@fit$coef["skew"]
  sk.h = fit.bch@fit$coef["skew"]

  #standardize residuals
  st.res.b = fit.btc@fit$residuals / fit.btc@fit$sigma
  st.res.e = fit.eth@fit$residuals / fit.eth@fit$sigma
  st.res.l = fit.ltc@fit$residuals / fit.ltc@fit$sigma
  st.res.x = fit.xrp@fit$residuals / fit.xrp@fit$sigma
  st.res.h = fit.bch@fit$residuals / fit.bch@fit$sigma

  #get integral transforms
  U1 = psged(st.res.b, nu = nu.b,xi=sk.b)
  U2 = psged(st.res.e, nu = nu.e,xi=sk.e)
  U3 = psged(st.res.l, nu = nu.l,xi=sk.l)
  U4 = psged(st.res.x, nu = nu.x,xi=sk.x)
  U5 = psged(st.res.h, nu = nu.h,xi=sk.h)

  UU=cbind(U1,U2,U3,U4,U5)
  U=apply(UU,2,rank)/(dim(UU)[1]+1)
  
  est_rsc=Est_comp_C(U,2)
  Qvec=as.vector(est_rsc$Q)
  Rvec=as.vector(apply(est_rsc$R,3,function(x)x[upper.tri(x)]))
  cop.pars=c(est_rsc$init,Qvec,Rvec,est_rsc$nu)
  return(cop.pars)
}

l1=round(dim(x)[1]^(2/3))
#nboot=1000
nboot=2

set.seed(1)
startp2=Sys.time()
param_boot_2p <- tsboot(x, parametric_RScop_boot, R = nboot, l = l1, sim = "geom",
                        parallel = "multicore", ncpus = detectCores()-1, n.sim=dim(x)[1])

endp2=Sys.time()
endp2-startp2

p2=param_boot_2p$t
SE_p2=apply(p2,2,sd)
reg=2

# Standard errors

# Initial probabilities
round(SE_p2[1:reg],3)

# Transition probabilities
round(SE_p2[(reg+1):(reg*reg+reg)],3)

# Dependence matrices
round(SE_p2[(reg*reg+reg+1):(length(SE_p2)-reg)],3)

# Number of degrees-of-freedom
round(tail(SE_p2,reg),3)


# param__3-states ---------------------------------------------------------

parametric_RScop_boot <- function(x) {
  # This function performs the bootstrap for the estimation of standard errors of RSStC model parameters
  # with 3 states, estimated adopting the parametric approach and assuming an 
  # ARMA(1,1)-GARCH(1,1) with skewed generalized error distribution model for the marginals
  
  # ARGUMENTS:
  # x is the matrix of observations
  
  # VALUES:
  # cop.pars is the vector of RSStC copula parameters
  
  #marginal parameters
  btc = x[,1]
  eth=x[,2]
  xrp=x[,3]
  ltc=x[,4]
  bch=x[,5]
  
  model = "sGARCH"
  spec = ugarchspec(
    variance.model = list(model = model, garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(1, 1), include.mean = FALSE) ,
    distribution.model = "sged"
  )
  
  fit.btc = ugarchfit(spec, btc)
  fit.eth = ugarchfit(spec, eth)
  fit.ltc = ugarchfit(spec, ltc)
  fit.xrp = ugarchfit(spec, xrp)
  fit.bch = ugarchfit(spec, bch)
  
  marg.pars=c(fit.btc@fit$coef,
              fit.eth@fit$coef,
              fit.ltc@fit$coef,
              fit.xrp@fit$coef,
              fit.bch@fit$coef)
  
  nu.b = fit.btc@fit$coef["shape"]
  nu.e = fit.eth@fit$coef["shape"]
  nu.l = fit.ltc@fit$coef["shape"]
  nu.x = fit.xrp@fit$coef["shape"]
  nu.h = fit.bch@fit$coef["shape"]
  
  sk.b = fit.btc@fit$coef["skew"]
  sk.e = fit.eth@fit$coef["skew"]
  sk.l = fit.ltc@fit$coef["skew"]
  sk.x = fit.xrp@fit$coef["skew"]
  sk.h = fit.bch@fit$coef["skew"]
  
  #standardize residuals
  st.res.b = fit.btc@fit$residuals / fit.btc@fit$sigma
  st.res.e = fit.eth@fit$residuals / fit.eth@fit$sigma
  st.res.l = fit.ltc@fit$residuals / fit.ltc@fit$sigma
  st.res.x = fit.xrp@fit$residuals / fit.xrp@fit$sigma
  st.res.h = fit.bch@fit$residuals / fit.bch@fit$sigma
  
  #get integral transforms
  U1 = psged(st.res.b, nu = nu.b,xi=sk.b)
  U2 = psged(st.res.e, nu = nu.e,xi=sk.e)
  U3 = psged(st.res.l, nu = nu.l,xi=sk.l)
  U4 = psged(st.res.x, nu = nu.x,xi=sk.x)
  U5 = psged(st.res.h, nu = nu.h,xi=sk.h)
  
  UU=cbind(U1,U2,U3,U4,U5)
  U=apply(UU,2,rank)/(dim(UU)[1]+1)
  
  #RS copula param
  est_rsc=Est_comp_C(U,3)
  Qvec=as.vector(est_rsc$Q)
  Rvec=as.vector(apply(est_rsc$R,3,function(x)x[upper.tri(x)]))
  #Rvec=as.vector(est_rsc$R)
  cop.pars=c(est_rsc$init,Qvec,Rvec,est_rsc$nu)
  return(cop.pars)
}

l1=round(dim(x)[1]^(2/3))
#nboot=1000
nboot=2

set.seed(1)
startp=Sys.time()
param_boot_3p <- tsboot(x, parametric_RScop_boot, R = nboot, l = l1, sim = "geom",
                        parallel = "multicore", ncpus = detectCores()-1, n.sim=dim(x)[1])

endp=Sys.time()
endp-startp

p3=param_boot_3p$t
SE_p3=apply(p3,2,sd)
reg=3

# Standard errors:

# Initial probabilities
round(SE_p3[1:reg],3)

# Transition probabilities
round(SE_p3[(reg+1):(reg*reg+reg)],3)

# Dependence matrices
round(SE_p3[(reg*reg+reg+1):(length(SE_p3)-reg)],3)

# Number of degrees-of-freedom
round(tail(SE_p3,reg),3)

