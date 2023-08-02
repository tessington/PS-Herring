# code to debug hierarchical model

#include <TMB.hpp>

y <- data$y
t <- data$t
nstock <- data$nstock
ustockid <- data$ustockid
ystockid <- data$ystockid

beta <- fixed[1:4]
logsigma_proc <- fixed[5:8]
logphi <- fixed[which(rownames(fixed)=="logphi")]
thetaf <- fixed[which(rownames(fixed)=="logphi")]


    # DATA_VECTOR(y); // observed Y
    # DATA_IVECTOR(t); // index Y into u
    # DATA_INTEGER(nstock); // number of stocklets to model
    # DATA_IVECTOR(ustockid); // stock indices for process
    # DATA_IVECTOR(ystockid); // stock indices for observations
    # // parameter inputs
    # PARAMETER_VECTOR(beta);
    # //PARAMETER(logsigma_beta);
    # //PARAMETER_VECTOR(eps_beta);
    # PARAMETER_VECTOR(u); // true log(Nt)
    # PARAMETER_VECTOR(logsigma_proc);
    # PARAMETER_VECTOR(logphi);
    # PARAMETER_VECTOR(thetaf);
    # 
#    // Transform and setup parameters
#    //Type sigma_beta = exp(logsigma_beta);
    #int n = u.size();
n <- length(u)
    #int nobs = y.size();
nobs <- length(y)

sigma_proc <- exp(logsigma_proc)
p <- phi <- sigma_proc <- numeric(nstock)

#    vector<Type> sigma_proc(nstock);
#    vector<Type> p(nstock);
#    vector<Type> phi(nstock);
#    // Make arrays for sigma_proc, p, and phi
    
#    // parameter transformations
    for( i in 1:nstock) {
      sigma_proc[i] = exp(logsigma_proc[i]);
      p[i] = 1 + exp(thetaf[i])  / (1 + exp(thetaf[i]));
      phi[i] = exp(logphi[i]);
    }
    
    
    # // monitor key parmaters
    # ADREPORT(sigma_proc)
    # ADREPORT(p)
    # ADREPORT(phi)
    # 
    # // process model
    # 
    #Type nll = 0.0;
    nllproc = numeric(n)
    m <- numeric(n)
    
    for (i in 2:n) {
      if(ustockid[i] == ustockid[i-1]) {
        id = ustockid[i]+1
        m[i] = beta[id] + u[i-1];
        sigma2use = sigma_proc[id];
        nllproc[i] = - dnorm(u[i], m[i], sigma2use, TRUE);
      }
    }
    nllobs <- numeric(nobs)
    m <- numeric(nobs)
    for ( i in 1:nobs) {
      id = ystockid[i]+1
      u2use <- t[i]+1
      m[i] = exp(u[u2use])
      phi2use = phi[id];
      p2use = p[id];
      nllobs[i] = -log(dtweedie(y = y[i], mu = m[i], phi = phi2use, power = p2use))
    }
    
      
    
    sum(nllproc)+sum(nllobs)
  }

