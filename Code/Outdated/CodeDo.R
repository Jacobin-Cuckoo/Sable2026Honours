library(rstan)
library(DynForest)

options(mc.cores = parallel::detectCores())

data("pbc2")
for (i in 1:nrow(pbc2)) {
	if (pbc2$event[i] == 1) {
		pbc2$event[i] <- 0
	}
	if (pbc2$event[i] == 2) {
		pbc2$event[i] <- 1
	}
}

jm <- {"

// block definging functions
functions {

	/* Deterministic model */

	real SLD(real t, row_vector eta, vector mu, int nb_param) {
		// t is time
		// eta is random effects row vector
		// mu is fixed effect vector
		// nb_param is number of parameters
		vector[nb_param] psi; // psi[1] is TS_0 (intitial no. cancer cells)
													// psi[2] is d (decay rate)
													// psi[3] is g (growth rate)
													// psi[4] is phi (ratio of respondant cells)
		
		real sld;
		
		for(j in 1:(nb_param - 1)) { // for TS_0, d, g
			psi[j] = exp(log(mu[j]) + eta[j]);
		}

		psi[nb_param] = inv_logit(logit(mu[nb_param]) + eta[nb_param]); // for phi

		sld = psi[1] * (psi[4] * exp(-psi[2] * t) + (1 - psi[4]) * exp(psi[3] * t));

		return sld;
	}


	/* baseline exponetial hazard model*/
	real h0(real lambda) {
		real hr;
		hr = 1 / lambda;
		return hr;
	}


	/* hazard function
		 defined as integrable function */

	real[] h_int(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
		// theta is input parameters
		real SLDtemp;
		real h;
		real dSdt[1]; // Why array? idk
		real lambda = theta[1];
		real beta = theta[2]; // link parameter
		int nb_param = x_i[1];

		vector[nb_param] mu; // fixed effects
		row_vector[nb_param] eta; //random effects

		for(i in 1:(nb_param)) {
			mu[i] = theta[i + 2]; // theta[3 to 6]
			eta[i] = theta[i + 6]; // theta[7 to 10]
		}

		
		h = h0(lambda);
		SLDtemp = fmin(SLD(t, eta, mu, nb_param), 8000);
		// above \"avoids potential issues\"
		dSdt[1] = h * exp(beta * SLDtemp);

		return(dSdt);
	}

}



data {

	// Longitudinal data
	int nb_obs; // number of total observations
	int<lower = 1> N; // number of patients
	int nb_param;
	int nb_meas[N]; // number of measurements for each patient
	int<lower = 0> nb_max_meas; // max of above
	vector[nb_obs] measures; // measurement data
	vector[nb_max_meas] day; // day of measurement(63 day intervals)

	// Survival data
	real ToE[N, 1]; // why 2-D array??? Censoring time or time of event
	vector[N] delta; // vector of censoring indicators
	
	// integration parameters
	real t_0;
	real y_0[1];

}



transformed data {
	// integration parameters
	real x_r[0];
	real y[0];
}


parameters {
	// population parameters (fixed effects)
	real <lower = 0> mu_1; // mu_BSLD
	real <lower = 0> mu_2; // mu_d
	real <lower = 0> mu_3; // mu_g
	real <lower = 0, upper = 1> mu_4; // mu_phi

	vector<lower = 0>[nb_param] Omega; // Omega diagonal vector
	real<lower = 0> sigma; // sd for error
	real<lower = 0> lambda; //exponetial survival parameter
	real beta; // link parameter


	// random effects
	matrix[N, nb_param] eta_tilde;
}


transformed parameters{
	vector[nb_param] mu;
	// non-centered reparameterisation for hierarchical models
	matrix[N, nb_param] eta;
	for(j in 1:nb_param) {
		eta[,j] = eta_tilde[,j] * Omega[j];
	}

	mu[1] = mu_1;
	mu[2] = mu_2;
	mu[3] = mu_3;
	mu[4] = mu_4;

}



model {
	int w;
	w = 0;
	real theta[10];
	real SLDtemp;
	int xi[1]; // why array???
	xi[1] = nb_param;
	

	// population parameters
	theta[1] = lambda;
	theta[2] = beta;
	theta[3] = mu[1];
	theta[4] = mu[2];
	theta[5] = mu[3];
	theta[6] = mu[4];


	// model log-likelihood
	for(i in 1:N) {  // For each patient

		// random effects of patient i
		theta[7] = eta[i, 1];
		theta[8] = eta[i, 2];
		theta[9] = eta[i, 3];
		theta[10] = eta[i, 4];

		for(t in 1:nb_meas[i]) { // for each measurement for patient i

			SLDtemp = SLD(day[t], eta[i,], mu, nb_param);

			// density of longitudinal process with observed event

				measures[w + t] ~ normal(SLDtemp, SLDtemp*sigma);
		}

		w = w + nb_meas[i];

		// prior on random effects
		eta_tilde[i, ] ~ normal(0, 1);

		
		// density of survival process for observed death
		if(delta[i] == 1) {
			target += log(h_int(ToE[i, 1], y, theta, x_r, xi)[1] * 
								exp(-integrate_ode_rk45(h_int, y_0, t_0, ToE[i, ], theta,
										x_r, xi)[1, 1]));
		}

		// density of survival process with censoring
		else {
			target += log(exp(-integrate_ode_rk45(h_int, y_0, t_0, ToE[i, ],
												theta, x_r, xi)[1, 1]));
		}

	}

	// priors
	mu_1 ~ lognormal(4.8, 1);
	mu_2 ~ beta(1, 100);
	mu_3 ~ beta(1, 100);
	mu_4 ~ beta(2, 4);
	Omega ~ lognormal(0, 1);
	sigma ~ lognormal(0, 1);
	lambda ~ lognormal(5, 1);
	beta ~ normal(0, 0.5);

}
"}

jm <- {"
data {
	int <lower = 1> nPat; // number of patients
	int <lower = 1> nObs; // number of observations
	int <lower = 1> P; // number of Fixed effects
	int <lower = 1> Q; // number of Random effects
	
	matrix[nObs, P] X;
	matrix[nObs, Q] Z;
	vector[nObs] y;
	array[nObs] int id;
