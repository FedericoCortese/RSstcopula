library(copula)

RScop.viterbi <- function(U, nc){
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