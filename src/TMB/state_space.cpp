#include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_INTEGER(n);
    DATA_INTEGER(nyears);
    DATA_VECTOR(y);
    DATA_MATRIX(Xij);
    DATA_MATRIX(Uij);
    DATA_IVECTOR(yearindex);
    
    
    PARAMETER_VECTOR(beta); // fixed effect coefficients
    PARAMETER_VECTOR(gamma); // random effect coefficients
    PARAMETER(logsigma_gamma); // standard deviation of random effects
    PARAMETER_VECTOR(vt); // the vector of vt's
    PARAMETER(logitrho);
    PARAMETER(logsigma_vt); // variation in r_t
    PARAMETER(logno); // initial log population size
    PARAMETER(logphi); // negative binomial dispersion parameters
    
    Type jnll = 0.0;
    
    // Transformed parameters
    
    Type sigma_gamma = exp(logsigma_gamma);
    Type sigma_v = exp(logsigma_vt);
    Type phi = exp(logphi);
    Type rho = 1 / (1 + exp( logitrho ) );
    Type nyears_minus_1 = nyears - 1;
    vector<Type> rt(nyears - 1);
    vector<Type> lognt(nyears);
    vector<Type> nt(nyears);
    
    ADREPORT(sigma_gamma);
    ADREPORT(sigma_v);
    ADREPORT(phi);
    ADREPORT(rho);
    
    rt(0) = vt(0);
    lognt(0) = logno;
    nt(0) = exp(lognt(0));
    
    for ( int i  = 1; i < nyears_minus_1; i++){ 
      rt(i) = rt(i - 1) * rho + pow(1 - pow(rho, 2), -2) * vt(i - 1);
      }
    
    for ( int i  = 1; i < nyears; i++){ 
      lognt(i) = lognt(i - 1) + rt(i - 1);
      nt(i) = exp(lognt(i));
    }
    ADREPORT(nt);
    ADREPORT(rt);
    
    
    // generate a vector of lognt that align with observations
    vector<Type> logntdata(n);
    for (int i =0; i < n; i++) {
      logntdata(i) = lognt(yearindex(i));
    }
    
    
    // random effects likelihood
    jnll -= sum(dnorm(gamma,0,sigma_gamma, true));
    jnll -= sum(dnorm(vt, 0, sigma_v, true)); 
     
    vector<Type>  mu(n);
    mu = exp( logntdata + Xij* beta + Uij * gamma) ;
    
    
    vector<Type> var(n);
    for(int i = 0; i < n; i++){ 
      var[i] = mu[i] + pow(mu[i],2) * pow (phi, -1);
      jnll -=  	dnbinom2(y[i], mu[i],var[i], true);
      
    }
    
    return jnll;
  }
