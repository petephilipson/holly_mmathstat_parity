library(rstan)
#Poisson model for HFD data
model_pois_69 <- "
data{
  int<lower=0> n;  // Sample size
  int<lower=0> y[n];  // Observed count
}
parameters {
  real<lower=0> lambda; 
}
model {
  y ~ poisson(lambda);
}
"
hfd_data_69 <- list(
  n = length(counts69vec),
  y = counts69vec
)

fit_pois_69 <- stan(
  model_code = model_pois_69,      # R program
  data = hfd_data_69,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 6000,            # total number of iterations per chain
  cores = 1
)

traceplot(fit_pois_69) #plot
mean(counts69vec) #true mean value 
get_posterior_mean(fit_pois_69) #model mean value

#New model to differentiate between countries
model_pois_2 <- "
data{
  int<lower=0> n;  // Sample size
  int<lower=0> y[n];  // Observed count
  int<lower=1> c; // Country by number
  int<lower=1> crep[n]; 
}
parameters {
  real<lower=0> lambda[c];
}
model {
  for (i in 1:n)
  y[i] ~ poisson(lambda[crep[i]]);
}
"
hfd_data_2 <- list(
  n = length(counts69vec),
  y = counts69vec,
  c = length(codenumbs),
  crep = codenumbsrep
)

fit_pois_2 <- stan(
  model_code = model_pois_2,      # R program
  data = hfd_data_2,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 6000,            # total number of iterations per chain
  cores = 1
)

get_posterior_mean(fit_pois_2) #these are close to the means we have directly from the data!! 
traceplot(fit_pois_2) #plot

#New model to allow for prior 
model_pois_p <- "
data{
  int<lower=0> n;  // Sample size
  int<lower=0> y[n];  // Observed count
  int<lower=1> c; // Country by number
  int<lower=1> crep[n]; // Repeated country numbers 
  int<lower=0> sdval; // Prior variance
  real<lower=0> meanval; // Overall mean from 1969 
}
parameters {
  real<lower=0> lambda[c];
}
model {
  lambda ~ lognormal(log(meanval), sdval);
  for (i in 1:n)
  y[i] ~ poisson(lambda[crep[i]]);
}
"
hfd_data_p <- list(
  n = length(counts69vec),
  y = counts69vec,
  c = length(codenumbs),
  crep = codenumbsrep,
  sdval = 1, #could multiply by 1.5
  meanval = countsmean1969
)

fit_pois_p <- stan(
  model_code = model_pois_p,      # R program
  data = hfd_data_p,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 6000,            # total number of iterations per chain
  cores = 1
)

get_posterior_mean(fit_pois_p)[1:23,1] #means only - close - some closer than without the prior, some worse (marginally)
traceplot(fit_pois_p, col="blue") #plot

priorfit = as.matrix(fit_pois_p) #converts to matrix 

# colMeans(priorfit) = get_posterior_mean

apply(priorfit,2,sd) #standard dev of each lambda (per country)


#### INITIAL MODEL FOR PRIOR 1969 ####
{
hfd_model <- "
data{
  int<lower=0> n;  // Sample size
  int<lower=0> y[n];  // Observed count
  int<lower=1> c; // Country by number
  int<lower=1> crep[n]; // Repeated country numbers 
  real<lower=0> prior_sd[c]; // Prior standard deviation
  real<lower=0> prior_mean[c]; // Prior mean  
}
parameters {
  real<lower=0> lambda[c];
}
model {
  for (j in 1:c)
  lambda[j] ~ lognormal(log(prior_mean[j]), prior_sd[j]);
  for (i in 1:n)
  y[i] ~ poisson(lambda[crep[i]]);
}
"
} 
hfd_data1 <- list(
  n = length(counts69vec),
  y = counts69vec,
  c = length(codenumbs),
  crep = codenumbsrep,
  prior_sd = rep(1, length(codenumbs)), #could multiply by 1.5
  prior_mean = rep(countsmean1969, length(codenumbs))
)

res <- NULL

initial_fit <- stan(
  model_code = hfd_model,      # R program
  data = hfd_data1,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1
)

fit_mat <- as.matrix(initial_fit)
prior_mean <- apply(fit_mat, 2, mean)[1:23]
prior_sd <- apply(fit_mat, 2, sd)[1:23]

res <- rbind(res, fit_mat) #store samples

#model for the next year - here 1970
hfd_data2 <- list(
  n = length(counts70vec),
  y = as.numeric(counts70vec),
  c = length(codenumbs),
  crep = codenumbsrep,
  prior_sd = prior_sd, #could multiply by 1.5
  prior_mean = prior_mean
)

hfd_fit2 <- stan(
  model_code = hfd_model, # R program
  data = hfd_data2,     # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1
)

fit_mat <- as.matrix(hfd_fit2)
prior_mean <- apply(fit_mat, 2, mean)[1:23]
prior_sd <- apply(fit_mat, 2, sd)[1:23]

res <- rbind(res, fit_mat) #store samples










#### INITIAL MODEL FOR FIRST YEAR OF DATA = 1918 ####
init_model <- "
data{
  int<lower=0> n;  // Sample size
  int<lower=0> y[n];  // Observed count
  int<lower=1> c; // Country by number
  int<lower=1> crep[n]; // Repeated country numbers 
  real<lower=0> prior_sd[c]; // Prior standard deviation
  real<lower=0> prior_mean[c]; // Prior mean  
}
parameters {
  real lambda[c];
}
model {
  for (j in 1:c)
  lambda[j] ~ normal(log(prior_mean[j]), prior_sd[j]);
  for (i in 1:n)
  y[i] ~ poisson_log(lambda[crep[i]]);
}
"

m <- mean(as.numeric(as.matrix(all_years_hfd_counts[, 3:1002])))
init_data <- list(
  n = length(as.vector(as.matrix(subset(countsyearcode, cohort == 1918)[, 3:1002]))),
  y = as.numeric(as.vector(as.matrix(subset(countsyearcode, cohort == 1918)[, 3:1002]))),
  c = 25,
  crep = rep(subset(countsyearcode, cohort==1918)[,1003],1000),
  prior_sd = rep(1, 25),
  prior_mean = rep(m, 25)
)

res_pois <- NULL

initial_fit <- stan(
  model_code = init_model,      # R program
  data = init_data,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1
)

fit_mat <- as.matrix(initial_fit)
prior_mean <- apply(fit_mat, 2, mean)[1:25]
prior_sd <- apply(fit_mat, 2, sd)[1:25]

res_pois <- rbind(res_pois, fit_mat) #store samples

#### MODEL FOR REST OF DATA ####
loop_model <- "
data{
  int<lower=0> n;  // Sample size
  int<lower=0> y[n];  // Observed count
  int<lower=1> c; // Country by number
  int<lower=1> crep[n]; // Repeated country numbers 
  real<lower=0> prior_sd[c]; // Prior standard deviation
  real prior_mean[c]; // Prior mean  
}
parameters {
  real lambda[c];
}
model {
  for (j in 1:c)
  lambda[j] ~ normal(prior_mean[j], prior_sd[j]);
  for (i in 1:n)
  y[i] ~ poisson_log(lambda[crep[i]]);
}
"
date()
for(i in 1:64){
loop_data <- list(
  n = length(as.vector(as.matrix(subset(countsyearcode, cohort == 1918+i)[, 3:1002]))),
  y = as.numeric(as.vector(as.matrix(subset(countsyearcode, cohort == 1918+i)[, 3:1002]))),
  c = 25,
  crep = rep(subset(countsyearcode, cohort==1918+i)[,1003], 1000),
  prior_sd = prior_sd,
  prior_mean = prior_mean
)

loop_fit <- stan(
  model_code = loop_model,      # R program
  data = loop_data,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1
)


new_fit_mat <- as.matrix(loop_fit)
prior_mean <- apply(new_fit_mat, 2, mean)[1:25]
prior_sd <- apply(new_fit_mat, 2, sd)[1:25]

res <- rbind(res, new_fit_mat) #store samples
res_ll <- rbind(res_ll, new_fit_mat[, 26] - sum(lfactorial(loop_data$y)))
}
date()

