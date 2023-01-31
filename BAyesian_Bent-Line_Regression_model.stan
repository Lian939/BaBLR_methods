data {
	int<lower=0> N;								// number of observations
	int<lower=0> Npat;							// number of individuals	
	real y[N];								// outcome data
	real x[N];								// time data
	int<lower=1,upper=Npat> id[N];						// patient id for each observation
	real betacp_lower;							// lower bound for (prior on) mean change point
	int<lower=0> cp_mean;                                                  // (prior on) mean change point
	int<lower=0> cp_sd;                                                    // (prior on) sd change point
}

parameters {
	vector[2] beta;				                                // fixed effects, intercept and pre slopes
	real<upper=0> beta3;                                                    // fixed effects, diff slope <=0
	real<lower=betacp_lower> betacp;                                       // fixed effects, changepoint (bounding specified to help convergence)
	vector<lower=0>[4] beta_sd;                                            // level 2 error sd (sds of the fix effects)
	vector<lower=0>[4] u_sd;                                               // level 2 error sd (sds of the random effects)
        real<lower=0> y_sd;                                                    // level 1 error sd
        real u1[Npat];                                                         // random effects (level 2) errors
        real u2[Npat];                                                         // random effects (level 2) errors
        real<upper=0> u3[Npat];                                                // random effects (level 2) errors
        real u4[Npat];                                                         // random effects (level 2) errors
}

transformed parameters {	
  real y_mu[N];                                                               // mean parameter based on regression equation

  //==========================
  // calculate random effects
  //==========================
  real<upper=0> alpha3[Npat];
  real alpha1[Npat];
  real alpha2[Npat];
  real alpha4[Npat];
  vector[4] u[Npat];                                                            // random effects (level 2) errors
  vector[4] alpha[Npat];                                                        // random effects

  
  for (i in 1:Npat) {
    alpha1[i]= beta[1] + u1[i];
    alpha2[i]= beta[2] + u2[i];
    alpha3[i]= beta3 + u3[i];
    alpha4[i] = betacp + u4[i];
  }
  
  
  u[:,1] = u1;
  u[:,2] = u2;
  u[:,3] = u3;
  u[:,4] = u4;
    
  alpha[:,1] = alpha1;
  alpha[:,2] = alpha2;
  alpha[:,3] = alpha3;
  alpha[:,4] = alpha4;


  //=====================
  // regression equation
  //=====================
  
  for (j in 1:N) {
    if (x[j] < alpha[id[j],4]) 
      y_mu[j]= alpha[id[j],1] + alpha[id[j],2] * (x[j] - alpha[id[j],4]);
    else  
      y_mu[j] = alpha[id[j],1] + (alpha[id[j],2]+ alpha[id[j],3]) * (x[j] - alpha[id[j],4]);   
  } 
  
}


model {

matrix[4, 4] Sigma; //Defining covariance matrix
  //========
  // priors
  //========
                       
  beta_sd[1] ~ cauchy(0,10);                    // prior: fixed effect sd, intercept
  beta_sd[2] ~ cauchy(0,1);                     // prior: fixed effect sd, slope before change
  beta_sd[3] ~ cauchy(0,5);                     // prior: fixed effect sd, slope difference between before and after change
  beta_sd[4] ~ cauchy(0,10);                   // prior: fixed effect sd, change point
  
   
  beta[1] ~ normal(0,pow(beta_sd[1],2));       // prior: fixed effect, intercept
  beta[2] ~ normal(0,pow(beta_sd[2],2));       // prior: fixed effect, initial slope
  beta3 ~ normal(0,pow(beta_sd[3],2));         // prior: fixed effect, slope difference
  betacp ~ normal(cp_mean,pow(beta_sd[4],2));  // prior: fixed effect, change point
  
  u_sd[1] ~ cauchy(0,10);                      // prior: random effect sd, intercept
  u_sd[2] ~ cauchy(0,1);                       // prior: random effect sd, slope before change
  u_sd[3] ~ cauchy(0,5);                       // prior: random effect sd, slope after change
  u_sd[4] ~ cauchy(0,10);                     // prior: random effect sd, change point
  
  y_sd ~ cauchy(0,10);                        // prior: level 1 error sd
  
	
  //=============================
  // random effects distribution
  //=============================
  
  for (i in 1:Npat) {
     u[i,1] ~ normal(0,pow(u_sd[1],2));
     u[i,2] ~ normal(0,pow(u_sd[2],2));
     u[i,3] ~ normal(0,pow(u_sd[3],2));
     u[i,4] ~ normal(0,pow(u_sd[4],2)); 
  }
 
                             
  
  //==================
  // model likelihood
  //==================
  
  y ~ normal(y_mu, y_sd); // likelihood for the observed data
  
}
