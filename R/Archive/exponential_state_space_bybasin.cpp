#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data inputs
  DATA_VECTOR(y); // observed Y
  DATA_IVECTOR(t); // index Y into u
  DATA_INTEGER(nstock); // number of stocklets to model
  DATA_IVECTOR(ustockid); // stock indices for process
  DATA_IVECTOR(ystockid); // stock indices for observations
  // parameter inputs
  PARAMETER_VECTOR(beta);
  PARAMETER(logsigma_beta);
  PARAMETER(mu_beta);
  PARAMETER_VECTOR(u); // true log(Nt)
  PARAMETER_VECTOR(logsigma_proc);
  PARAMETER_VECTOR(logphi);
  PARAMETER_VECTOR(thetaf);
  
  // Transform and setup parameters
  Type sigma_beta = exp(logsigma_beta);
  int n = u.size();
  int nobs = y.size();
  vector<Type> sigma_proc(nstock);
  vector<Type> p(nstock);
  vector<Type> phi(nstock);
  // Make arrays for sigma_proc, p, and phi
  
  // parameter transformations
  for(int i = 0; i < nstock; i++) {
  sigma_proc(i) = exp(logsigma_proc[i]);
  p(i) = 1 + exp(thetaf[i])  / (1 + exp(thetaf[i]));
  phi(i) = exp(logphi[i]);
  }
  
  
  // monitor key parmaters
  ADREPORT(sigma_beta)
  ADREPORT(sigma_proc)
  ADREPORT(p)
  ADREPORT(phi)
  
  
  
  // process model
  
  Type  nll=0;
  
  for(int i = 1; i < n; i++) {
    if (ustockid[i] == ustockid[i-1] ) {
    int id = ustockid[i];
    Type m = beta[id] + u[i-1];
    Type sigma2use = sigma_proc(id);
    nll -= dnorm(u[i], m, sigma2use, true);
    }
  }
  
  // observation model
  for(int i = 0; i < nobs; i++) {
    int id = ystockid[i];
    int u2use = t[i];
    Type m = exp(u[u2use]);
    Type phi2use = phi(id);
    Type p2use = p(id);
    nll -= dtweedie(y[i],m, phi2use, p2use, true);
  }
  
 // random effects of the betas
 for(int i = 0; i < nstock; i++) {
   nll -= dnorm(beta[i], mu_beta, sigma_beta, true);
 }
  return(nll);
}
  
