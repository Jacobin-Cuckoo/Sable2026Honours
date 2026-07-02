library(rstan)
options(mc.cores = parallel::detectCores())

data <- read.csv("data.csv")

jm <- "

// block definging functions
functions {

	/* Deterministic model */
	real SLD(real t, row_vector eta, vector mu, int nb_param) {
		vector[nb_param] psi; // psi[1] is TS_0 (intitial no. cancer cells)
													// psi[2] is d (decay rate)
													// psi[3] is g (growth rate)
													// psi[4] is phi (ratio of respondant cells)
		
		real sld;
		
		for(j in 1:(nb_param - 1)) {
			psi[j] = exp(log(mu[j]) + eta[j]);
		}

		psi[nb_param] = inv_logit(logit(mu[nb_param]) + eta[nb_param]);

		sld = psi[1] * (psi[4] * exp(-psi[2]*t) + (1 - psi[4]) * exp(psi[3] * t));

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
		real SLDtemp;
		real h;
		real dSdt[1]; // Why array? idk
		real lambda = theta[1];
		real beta = theta[2];

"

##### Formatting data

nb_obs <- nrow(data) # number of observations
N <- length(unique(data$ID)) # number of patients

nb_meas <- c() # number of measurements per patient
for (i in 1:N) {
	nb_meas[i] <- sum(data$ID == i)
}

nb_max_meas <- max(nb_meas) # max number of measurements for a patient

measures <- data$SLD # measurements
censor <- data$Censor_SLD # censoring indicator

day <- seq(0, max(data$Time), 63) # measuring days (measurements taken every 63 days)
t <- subset(data, Time == 0)$T # Time of event
t <- matrix(t, nrow = N, byrow = FALSE) # Make it a column matrix

delta <- subset(data, Time == 0)$delta # Event indicator
nb_param <- 4 # number of parameters
							#    

### integration parameters
t_0 <- 0
y_0 <- array(0)

##### list of stan function input parameters
dat <- list(nb_obs = nb_obs, N = N,
						nb_meas = nb_meas,
						nb_max_meas = nb_max_meas,
						measures = measures,
						censor = censor,
						day = day, t = t,
						delta = delta, t_0 = t_0,
						y_0 = y_0,
						nb_param = nb_param,
						c = 2.5)

# HMC algorithm parameters
nb_iter <- 400 # number of iterations
nb_warmup <- 200 # no. warm up iterations
nb_chains <- 3 #no. chains
nb_cores <- nb_chains # one core per chain

# control params
treedepth <- 12
adapt_delta <- 0.9

# initialization of population parameters and random effects
inits <- list(
	mu_1 = 60, mu_2 = 0.0055, mu_3 = 0.0015, mu_4 = 0.2,
	sigma = 0.18, Omega = c(0.7, 1, 1, 1.5), lambda = 580,
	beta = 0, eta = matrix(rep(0, N*nb_param), nrow = N, ncol = nb_param, byrow = FALSE)
)

