#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data inputs
  DATA_VECTOR(y); // observed Y
  DATA_IVECTOR(t); // vector of indices for the observations
  // parameter inputs
  PARAMETER(beta);
  PARAMETER_VECTOR(u); // true log(Nt)
  PARAMETER(logsigma_proc);
  PARAMETER(thetaf);
  PARAMETER(logphi);
  
  
  // parameter transformations
  Type sigma_proc = exp(logsigma_proc);
  Type p = 1 + exp(thetaf)  / (1 + exp(thetaf));
  Type phi = exp(logphi);
  int nobs = y.size(); // number of observations
  int n = u.size(); // number of years to estimate
  
  
  // monitor key parmaters
  ADREPORT(sigma_proc)
  ADREPORT(p)
  ADREPORT(phi)
  
  // process model
  
  Type nll = 0.0;
  
  for(int i = 1; i < n; i++) {
    Type m = beta + u[i-1];
    nll -= dnorm(u[i], m, sigma_proc, true);
  }
  
  // observation model
  for (int i = 0; i<nobs; i++){
    int yearindex = t[i];
    Type m =  exp(u[yearindex]);
    nll -=dtweedie(y(i), m, phi, p, true);
  }
  
  return(nll);
}
  
