#include <cmath>  // std::pow

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

using namespace roptim;

arma::vec ldtcC(arma::vec u, arma::mat R, double nu){
  //nu=exp(nu);
  int K=u.size();
  arma::vec z(K);
  for(int i=0;i<K;i++){
    z[i]=R::qt(u(i),nu,true,false);
  }
  arma::mat zz=trans(z);
  arma::vec c =-log(det(R))/2+ log(std::abs(tgamma((nu+K)/2)))
    -log(std::abs(tgamma(nu/2)))
    +K*(log(std::abs(tgamma(nu/2)))-log(std::abs(tgamma((nu+1)/2))))
    +(nu+1)/2*arma::accu(log(1+pow(z,2)/nu)) -(nu+K)/2*log(1+(zz*inv(R)*z)/nu) ;
    
    return c;
}

arma::vec dtcC(arma::vec u, arma::mat R, double nu){
  //nu=exp(nu);
  int K=u.size();
  arma::vec z(K);
  arma::vec den(K);
  for(int i=0;i<K;i++){
    z[i]=R::qt(u(i),nu,true,false);
    den[i]=pow(1+pow(z[i],2)/nu,-(nu+1)/2);
  }
  arma::mat zz=trans(z);
  arma::vec c=(tgamma((nu+K)/2)*pow(tgamma(nu/2),K-1))/(pow(det(R),.5)*pow(tgamma((nu+1)/2),K))*
    ((pow((1+(zz*inv(R)*z)/nu),-(nu+K)/2))/prod(den));
  return c;
}

double llktC(arma::mat U, arma::cube R, arma::vec nu, arma::mat Q, arma::vec init) {
  int r=Q.n_cols;
  int n=U.n_rows;
  int d=U.n_cols;
  
  double sv;
  
  arma::vec alpha(r);
  double L=0;
  int t, i;
  int h=0;
  int j=0;
  double temp;
  arma::mat dc(n,r);
  for(i=0;i<n;i++){
    for(j=0;j<r;j++){
      arma::mat Rt=R.slice(j);
      temp=dtcC(trans(U.row(i)),Rt,nu(j)).eval()(0,0);
      dc(i,j)=temp;
    }
  }
  alpha=init%trans(dc.row(0));
  sv=sum(alpha);
  L+=log(sv);
  alpha=alpha/sv;
  for(t=1;t<n;t++){
    alpha=Q*alpha%trans(dc.row(t));
    sv=sum(alpha);
    L+=log(sv);
    alpha=alpha/sv;
  }
  return L;
}


double obj_fun(double nu,arma::mat U,arma::mat R, arma::mat w, int k){
  int n=U.n_rows;
  double llk=0.0;
  int t;
  arma::vec temp;
  
  for(t=0;t<n;t++){
    temp=ldtcC(trans(U.row(t)),R,nu);
    llk+=temp.eval()(0,0)*w(t,k);
    //llk+=temp.eval()(0,0);
  }
  return -llk;
}

class Fun : public Functor {
public:
  Fun(arma::mat RR, arma::mat UU, arma::mat ww, int kk){
    U=UU;
    R=RR;
    w=ww;
    k=kk;
  }
  double operator()(const arma::vec &x) override {
    double nu=x(0);
    double llk = obj_fun(nu, U, R, w, k);
    return llk;
  }
private:
  arma::mat R;
  arma::mat U;
  arma::mat w;
  int k;
};

// [[Rcpp::export]]
Rcpp::List EMstepC(arma::mat U, arma::cube R, arma::vec nu, arma::mat Q, arma::vec init){
  
  int d=U.n_cols;
  int r=Q.n_cols;
  int n=U.n_rows;
  
  int i,j,k,h;
  
  //***************E step********************
  double tempd;
  arma::mat dc(n,r);
  for(i=0;i<n;i++){
    for(j=0;j<r;j++){
      arma::mat Rt=R.slice(j);
      tempd=dtcC(trans(U.row(i)),Rt,nu(j)).eval()(0,0);
      dc(i,j)=tempd;
    }
  }
  
  //fill last row of etabar
  arma::mat etabar(n,r);
  double r0=pow(r,-1);
  for(i=0; i<r;i++){
    etabar(n-1,i)=r0;
  }
  
  double sv;
  //fill etabar
  for(k=n-2;k>=0;k--){
    for(j=0;j<r;j++){
      etabar(k,j)=etabar(k+1,j)*dc(k+1,j);
    }
    etabar.row(k)=etabar.row(k)*Q.t();
    sv=sum(etabar.row(k));
    etabar.row(k)=etabar.row(k)*pow(sv,-1);
  }
  
  //fill first row of eta
  arma::mat eta(n,r);
  for(i=0;i<r;i++){
    eta(0,i)=init(i)*dc(0,i);
  }
  sv=sum(eta.row(0));
  eta.row(0)=eta.row(0)*pow(sv,-1);
  
  //fill eta
  arma::mat v(1,r);
  for(k=1;k<n;k++){
    v=eta.row(k-1)*Q;
    for(j=0;j<r;j++){
      eta(k,j)=v(j)*dc(k,j);
    }
    sv=sum(eta.row(k));
    eta.row(k)=eta.row(k)*pow(sv,-1);
  }
  
  //w computation
  arma::mat w(n,r);
  for(k=0;k<n;k++){
    for(j=0;j<r;j++){
      w(k,j)=eta(k,j)*etabar(k,j);
    }
    sv=sum(w.row(k));
    w.row(k)=w.row(k)*pow(sv,-1);
  }
  
  //z computation
  arma::cube z(r,r,n-1);
  arma::mat gc(n,r);
  for(k=0;k<n;k++){
    for(j=0;j<r;j++){
      gc(k,j)=etabar(k,j)*dc(k,j);
    }
  }
  arma::mat m(r,r);
  i=1;
  for(i=1;i<n;i++){
    m=trans(eta.row(i-1))*gc.row(i)%Q;
    double M=0.0;
    for(k=0;k<r;k++){
      for(j=0;j<r;j++){
        M+=m(k,j);
      }
    }
    z.slice(i-1)=m*pow(M,-1);
  }
  
  //***************M step********************
  arma::vec initnew(r);
  arma::mat Qnew(r,r);
  arma::vec qv(r);
  
  //initial probabilities
  initnew=trans(w.row(0));
  
  //transition probabilities
  for(i=0;i<r;i++){
    for(j=0;j<r;j++){
      double q=0.0;
      for(h=0;h<n-1;h++){
        q+=z(i,j,h);
      }
      qv(j)=q;
    }
    double sqv=sum(qv);
    Qnew.row(i)=pow(sqv,-1)*trans(qv);
  }
  
  //copula parameters
  arma::mat Z(n,d);
  arma::mat temp;
  arma::mat Rt(d,d);
  arma::cube Rnew(d,d,r);
  arma::vec nunew(r);
  
  for(k=0;k<r;k++){
    
    Rt=R.slice(k);
    arma::mat R0i=inv(Rt);
    for(i=0;i<n;i++){
      for(j=0;j<d;j++){
        Z(i,j)=R::qt(U(i,j),nu[k],true,false);
      }
    }
    double W=0;
    arma::mat A(d,d);
    for(i=0;i<n;i++){
      temp=Z.row(i)*R0i*trans(Z.row(i));
      double den =1+temp.eval()(0,0)/nu[k];
      A=A+w(i,k)*((trans(Z.row(i))*Z.row(i))*(1/den));
      W+=w(i,k);
    }
    A=A*(nu[k] + d) / (W * nu[k]) ;
    arma::vec a(d);
    for(i=0;i<d;i++){
      a(i)=1-A(i,i);
    }
    arma::mat Rtnew(d,d);
    Rtnew=A + Rt*diagmat(inv(pow(Rt,2))*a)*Rt;
    Rnew.slice(k)=Rtnew;
    
    
    arma::vec lower(1);
    arma::vec upper(1);
    // lower[0]=log(2);
    // upper[0]=log(100);
    //nu0=nu[k];
    lower[0]=2;
    upper[0]=50;
    Fun rb(Rtnew,U,w,k);
    Roptim<Fun> opt("L-BFGS-B");
    //Roptim<Fun> opt("SANN");
    opt.set_lower(lower);
    opt.set_upper(upper);
    opt.control.trace = 0;
    arma::vec x(1);
    //x(0)=log(nu[k]);
    x(0)=nu[k];
    opt.minimize(rb, x);
    arma::vec optnu=opt.par();
    //nunew[k]=exp(optnu(0));
    nunew[k]=optnu(0);
  }
  
  Rcpp::List out(5);
  out["R"]=Rnew;
  out["nu"]=nunew;
  out["Q"]=Qnew;
  out["init"]=initnew;
  out["w"]=w;
  out["z"]=z;
  
  return out;
}


double parmax(arma::cube R, arma::cube R0, arma::vec nu, arma::vec nu0, arma::mat Q,
              arma::mat Q0, arma::vec init, arma::vec init0) {
  
  arma::vec mxs(4);
  mxs[0]=max(abs(init-init0));
  mxs[1]=max(abs(nu-nu0));
  arma::vec temp=max(abs(R-R0));
  mxs[2]= max(temp);
  arma::vec temp2=trans(max(abs(Q-Q0)));
  mxs[3]= max(temp2);
    
  
  return max(mxs);
}

// [[Rcpp::export]]
Rcpp::List RSest(arma::mat U, int r, arma::cube R0, arma::vec nu0,
                 arma::mat Q0, arma::vec init0,
                 int maxiter=1000, int ninit=1, double eps=0.00000001){
  
  int n=U.n_rows;
  int d=U.n_cols;
  arma::mat Q;
  arma::vec init;
  arma::cube R;
  arma::vec nu;
  arma::mat w(n,r);
  arma::cube z(r,r,n-1);
  Rcpp::List est;
  double mx, llk, llk0;
  int i;
  
  for(i=0; i<ninit;i++){
    est=EMstepC(U, R0, nu0, Q0, init0);
    R=Rcpp::as<arma::cube>(est["R"]);
    nu=Rcpp::as<arma::vec>(est["nu"]);
    Q=Rcpp::as<arma::mat>(est["Q"]);
    init=Rcpp::as<arma::vec>(est["init"]);
    R0=R;
    nu0=nu;
    Q0=Q;
    init0=init;
    w=Rcpp::as<arma::mat>(est["w"]);
    z=Rcpp::as<arma::cube>(est["z"]);
  }
  llk0=llktC(U,R,nu,Q,init);
  
  for(i=0;i<maxiter;i++){
    est=EMstepC(U, R0, nu0, Q0, init0);
    R=Rcpp::as<arma::cube>(est["R"]);
    nu=Rcpp::as<arma::vec>(est["nu"]);
    Q=Rcpp::as<arma::mat>(est["Q"]);
    init=Rcpp::as<arma::vec>(est["init"]);
    w=Rcpp::as<arma::mat>(est["w"]);
    z=Rcpp::as<arma::cube>(est["z"]);
    mx=parmax(R,R0,nu,nu0,Q,Q0, init, init0);
    llk=llktC(U,R,nu,Q,init);
    if(mx<eps||((llk-llk0)/std::abs(llk0))<eps){
      break;
    }
    llk0=llk;
    R0=R;
    nu0=nu;
    Q0=Q;
    init0=init;
  }
  
  Rcpp::List out;
  out["R"]=R;
  out["nu"]=nu;
  out["Q"]=Q;
  out["init"]=init;
  out["iter"]=i+1;
  out["w"]=w;
  out["z"]=z;
  out["llk"]=llk;
  
  return out;
  
}
