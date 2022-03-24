
# semiparam__2-states -----------------------------------------------------
library(bootstrap)
library(boot)
library(parallel)
library(rugarch)
library(fGarch)
library(skewt)


semiparametric_RScop_boot <- function(x) {
  #This function performs the bootstrap for the estimation of standard errors of RSStC model with 2 states
  #adopting the semiparametric approach
  #x is the matrix of observations
  
  #marginals
  R=apply(x,2,rank)/(dim(x)[1]+1)
  
  #RS copula param
  est_rsc=Est_comp_C(R,2)
  Qvec=as.vector(est_rsc$Q)
  Rvec=as.vector(apply(est_rsc$R,3,function(x)x[upper.tri(x)]))
  cop.pars=c(est_rsc$init,Qvec,Rvec,est_rsc$nu)
  return(cop.pars)
}

l1=round(dim(x)[1]^(2/3))
#nboot=1000
nboot=100
library(bootstrap)
library(boot)
startsp2=Sys.time()
param_boot_2sp <- tsboot(x, semiparametric_RScop_boot, R = nboot, l = l1, sim = "geom",
                         parallel = "multicore", ncpus = detectCores(), n.sim=dim(x)[1])

endsp2=Sys.time()
endsp2-startsp2

# semiparam__3-states -----------------------------------------------------

library(bootstrap)
#theta <- function(x){mean(x)} 
#results = bootstrap(testFile,100,theta) 
library(boot)
library(parallel)

semiparametric_RScop_boot <- function(x) {
  #This function performs the bootstrap for the estimation of standard errors of RSStC model with 3 states
  #adopting the semiparametric approach
  #x is the matrix of observations
  
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
nboot=1000
library(bootstrap)
library(boot)
startsp2=Sys.time()
param_boot_3sp <- tsboot(x, semiparametric_RScop_boot, R = nboot, l = l1, sim = "geom",
                         parallel = "multicore", ncpus = detectCores(), n.sim=dim(x)[1])

endsp2=Sys.time()
endsp2-startsp2


# param__2-states ---------------------------------------------------------
library(rugarch)
library(fGarch)
library(skewt)
parametric_RScop_boot <- function(x) {
  #This function performs the bootstrap for the estimation of standard errors of RSStC model with 2 states
  #adopting the parametric approach
  #x is the matrix of observations: marginal distribution are here specified in terms of an ARMA(1,1)-GARCH(1,1)
  #model with skewed GED for the error terms
  
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
  #sourceCpp("RSest_trede_C.cpp")
  est_rsc=Est_comp_C(U,2)
  Qvec=as.vector(est_rsc$Q)
  Rvec=as.vector(apply(est_rsc$R,3,function(x)x[upper.tri(x)]))
  #Rvec=as.vector(est_rsc$R)
  cop.pars=c(est_rsc$init,Qvec,Rvec,est_rsc$nu)
  return(c(cop.pars,marg.pars))
}

l1=round(dim(U)[1]^(2/3))
#nboot=1000
nboot=100

library(bootstrap)
library(boot)
startp2=Sys.time()
param_boot_2p <- tsboot(U, parametric_RScop_boot, R = nboot, l = l1, sim = "geom",
                        parallel = "multicore", ncpus = detectCores(), n.sim=dim(x)[1])

endp2=Sys.time()
endp2-startp2


# param__3-states ---------------------------------------------------------
library(rugarch)
library(fGarch)
library(skewt)
parametric_RScop_boot <- function(x) {
  #This function performs the bootstrap for the estimation of standard errors of RSStC model with 3 states
  #adopting the parametric approach
  #x is the matrix of observations: marginal distribution are here specified in terms of an ARMA(1,1)-GARCH(1,1)
  #model with skewed GED for the error terms
  
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
  return(c(cop.pars,marg.pars))
}

l1=round(dim(U)[1]^(2/3))
#nboot=1000
nboot=100

library(bootstrap)
library(boot)
startp=Sys.time()
param_boot_3p <- tsboot(U, parametric_RScop_boot, R = nboot, l = l1, sim = "geom",
                        parallel = "multicore", ncpus = detectCores(), n.sim=dim(x)[1])

endp=Sys.time()
endp-startp