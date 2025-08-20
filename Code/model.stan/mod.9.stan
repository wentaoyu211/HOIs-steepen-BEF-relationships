data{
  int ND; // number of total data points
  int NS; // number of species
  int NP; // number of plots in total
  int NY; // number of years
  int NN; // number of neighbours
  int N_spy; // number of unique species * year
  int N_spp; // number of unique species * plot
  real div_ave;
  
  // biomass data
  real focal_bm[ND];  // biomass of focal tree
  real next_bm[ND];   // the next year biomass of corresponding focal tree
  real neib_bm[ND, NN]; // biomass of neighboring trees
  
  // fixed parameters from Yu & Albert et al 2024
  real<lower=0> beta[NS];  // species specific growth parameter
  real theta; 			   //MTE exponent
  matrix[NS,NS] alpha;    // interaction strengths
  real<lower=0> b;         // general exponent of neighbor tree biomass  
  
  // indices
  int fdiv[ND];  // need to take log2(div) centered around log2(4) 
  int neib_sp_ID[ND, NN]; //species ID for neighboring trees  
  int fsp_ID[ND]; // identify species of each data point
  int year_ID[ND];   // identify the collecting year for each data point
  int spy[ND];    // species * year indices
  int fplot_ID[ND]; // identify plot for each data point
  int spp[ND];      // species * plot indices
}

parameters{
  matrix[NS,NS] gamma; // species-pair specific slope for each insteraction ~ div relationship
  
  // random effects
  vector[NP] eps_p_raw;   // model plots as random effect
  real<lower=0> sigma_eps_p;
  vector[N_spp] eps_spp_raw;   // model species * plots
  real<lower=0> sigma_eps_spp;   
  vector[NY] eps_y_raw;   // model year as random effect
  real<lower=0> sigma_eps_y;
  vector[N_spy] eps_spy_raw;   // model species * year 
  real<lower=0> sigma_eps_spy;
  
  real<lower=0> sigma;    // sigma is observation error   
}

transformed parameters{
  vector[NP] eps_p;
  vector[NY] eps_y; 
  vector[N_spy] eps_spy;
  vector[N_spp] eps_spp;
  
  //non-centred parameterization to improve efficiency and convergence of hierachical model
  eps_p = sigma_eps_p * eps_p_raw;
  eps_y = sigma_eps_y * eps_y_raw;
  eps_spy = sigma_eps_spy * eps_spy_raw;
  eps_spp = sigma_eps_spp * eps_spp_raw;
}

model{
  // for model simplification
  real interactions;  // interaction effect
  real predBM;        // predicted tree biomass
  
  // specify priors
  for (i in 1:NS){
	  for (j in 1:NS){
	  gamma[i, j] ~ normal(0,1);
	  }	
	}
  
  // random effect structure
  sigma_eps_p ~ normal(0, 0.1);
  sigma_eps_y ~ normal(0, 0.1);
  sigma_eps_spy ~ normal(0, 0.1);  
  sigma_eps_spp ~ normal(0, 0.1); 
  eps_p_raw ~ std_normal();
  eps_y_raw ~ std_normal();
  eps_spy_raw ~ std_normal();
  eps_spp_raw ~ std_normal();
   
  for(i in 1:ND){    // loop over all data points
        interactions = 0;
        for(j in 1:NN){    // sum over all neighbors
          interactions = interactions + alpha[fsp_ID[i], neib_sp_ID[i,j]]*((fdiv[i]/div_ave)^gamma[fsp_ID[i], neib_sp_ID[i,j]]) * neib_bm[i,j]^b;	
        }
        // prediction
        predBM = focal_bm[i] + max({beta[fsp_ID[i]]*focal_bm[i]^theta * max({1 + eps_p[fplot_ID[i]] + eps_spp[spp[i]] + eps_y[year_ID[i]] + eps_spy[spy[i]], 0}) + interactions, 0}); 
		
        // likelihood		
		next_bm[i] ~ lognormal( log(predBM), sigma );
  }
}

generated quantities {
  real interactions;
  real predBM; 
  vector[ND] log_lik;  // number of log likelihood for loo
  
  for (i in 1:ND){    // loop through focal individual tree
        interactions = 0;
        for(j in 1:NN){    // sum over all neighbors
        interactions = interactions + alpha[fsp_ID[i], neib_sp_ID[i,j]]*((fdiv[i]/div_ave)^gamma[fsp_ID[i], neib_sp_ID[i,j]]) * neib_bm[i,j]^b;	  
        }
      // prediction
      predBM = focal_bm[i] + max({beta[fsp_ID[i]]*focal_bm[i]^theta * max({1 + eps_p[fplot_ID[i]] + eps_spp[spp[i]] + eps_y[year_ID[i]] + eps_spy[spy[i]], 0}) + interactions, 0}); 
	  // pointwise log likelihood
	  log_lik[i] = lognormal_lpdf( next_bm[i] | log(predBM), sigma);
  }
}
