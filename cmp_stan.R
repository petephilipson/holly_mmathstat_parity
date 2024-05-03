#### Conway-Maxwell-Poisson for HFD data ####
model_CMP_int <- '
functions {
real log_Z_com_poisson(real log_mu, real nu, vector unique, vector logfac) {
  real log_Z;
  if (nu == 1) {
    return exp(log_mu);
  }
  // nu == 0 or Inf will fail in this parameterization
  if (nu <= 0) {
    reject("nu must be positive");
  }
  if (nu == positive_infinity()) {
    reject("nu must be finite");
  }
  // direct computation of the truncated series
  log_Z = log_sum_exp(unique*log_mu - nu*logfac);
  return log_Z;
}
real com_poisson_log_lpmf(int y, real log_mu, real nu, vector unique, vector logfac) {
  if (nu == 1) return poisson_log_lpmf(y | log_mu);
  return (y * log_mu - nu*lgamma(y + 1)) - log_Z_com_poisson(log_mu, nu, unique, logfac);
}
real com_poisson_lpmf(int y, real mu, real nu, vector unique, vector logfac) {
  if (nu == 1) return poisson_lpmf(y | mu);
  return com_poisson_log_lpmf(y | log(mu), nu, unique, logfac);
}
}
data {
  int<lower=0> n;  // Sample size
  int<lower=0> y[n];  // Observed count
  vector[11] unique; // Unique values for count
  vector[11] logfac; // Log-factorial for each unique count
}
parameters {
  real alpha; 
  real<lower=0> beta; 
}
model {
  for (i in 1:n){
    y[i] ~ com_poisson_log(alpha, beta, unique, logfac);
  }
}
'

ii1969 <- all_years_hfd_counts$cohort == 1969
counts69m <- as.matrix(all_years_hfd_counts[ii1969, -(1:2)])
counts69vec <- as.numeric(counts69m)

hfd_data_69_cmp <- list(
  n = length(counts69vec),
  y = counts69vec,
  unique = sort(unique(counts69vec)),
  logfac = lfactorial(sort(unique(counts69vec)))
)

fit_CMP_int <- stan(
  model_code = model_CMP_int,      # R program
  data = hfd_data_69_cmp,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1
)

traceplot(fit_CMP_int) #plot

#### Now include countries ####
model_CMP_int_country <- '
functions {
real log_Z_com_poisson(real log_mu, real nu, vector unique, vector logfac) {
  real log_Z;
  if (nu == 1) {
    return exp(log_mu);
  }
  // nu == 0 or Inf will fail in this parameterization
  if (nu <= 0) {
    reject("nu must be positive");
  }
  if (nu == positive_infinity()) {
    reject("nu must be finite");
  }
  // direct computation of the truncated series
  log_Z = log_sum_exp(unique*log_mu - nu*logfac);
  return log_Z;
}
real com_poisson_log_lpmf(int y, real log_mu, real nu, vector unique, vector logfac) {
  if (nu == 1) return poisson_log_lpmf(y | log_mu);
  return (y * log_mu - nu*lgamma(y + 1)) - log_Z_com_poisson(log_mu, nu, unique, logfac);
}
real com_poisson_lpmf(int y, real mu, real nu, vector unique, vector logfac) {
  if (nu == 1) return poisson_lpmf(y | mu);
  return com_poisson_log_lpmf(y | log(mu), nu, unique, logfac);
}
}
data {
  int<lower=0> n;  // Sample size
  int<lower=0> nc;  // Number of countries
  int<lower=0> y[n];  // Observed count
  int<lower=1> cid[n];  // Country index
  vector[11] unique; // Unique values for count
  vector[11] logfac; // Log-factorial for each unique count
  real<lower=0> prior_lam_m; // Prior mean for lambda 
  real<lower=0> prior_lam_s; // Prior sd for lambda
  real prior_nu_m; // Prior mean for nu (on log-scale) 
  real<lower=0> prior_nu_s; // Prior sd for nu
}
parameters {
  real lambda[nc]; 
  real nu[nc]; 
}
model {
  lambda ~ normal(prior_lam_m, prior_lam_s);
  nu ~ normal(0, prior_nu_s);
  for (i in 1:n){
    y[i] ~ com_poisson_log(lambda[cid[i]], exp(nu[cid[i]]), unique, logfac);
  }
}
'

# Get number of countries in given year
nc <- length(table(all_years_hfd_counts$country[ii1969]))
# Get vector of country ids in given year
cid <- rep(1:nc, 1000)
hfd_data_69_cmp_2 <- list(
  n = length(counts69vec),
  y = counts69vec,
  unique = sort(unique(counts69vec)),
  logfac = lfactorial(sort(unique(counts69vec))),
  nc = nc,
  cid = cid,
  prior_lam_m = mean(counts69vec),
  prior_lam_s = 1,
  prior_nu_m = 1,
  prior_nu_s = 1
)

fit_CMP_int_country <- stan(
  model_code = model_CMP_int_country,      # R program
  data = hfd_data_69_cmp_2,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1
)

traceplot(fit_CMP_int_country, col = "navy") #plot
traceplot(fit_CMP_int_country, col = "navy", pars = "nu") # Dispersions

# Next the sequential part
# Borrow from Poisson model
# Model code with vectors for priors now
#### Now include countries ####
model_CMP_int_country <- '
functions {
real log_Z_com_poisson(real log_mu, real nu, vector unique, vector logfac) {
  real log_Z;
  if (nu == 1) {
    return exp(log_mu);
  }
  // nu == 0 or Inf will fail in this parameterization
  if (nu <= 0) {
    reject("nu must be positive");
  }
  if (nu == positive_infinity()) {
    reject("nu must be finite");
  }
  // direct computation of the truncated series
  log_Z = log_sum_exp(unique*log_mu - nu*logfac);
  return log_Z;
}
real com_poisson_log_lpmf(int y, real log_mu, real nu, vector unique, vector logfac) {
  if (nu == 1) return poisson_log_lpmf(y | log_mu);
  return (y * log_mu - nu*lgamma(y + 1)) - log_Z_com_poisson(log_mu, nu, unique, logfac);
}
real com_poisson_lpmf(int y, real mu, real nu, vector unique, vector logfac) {
  if (nu == 1) return poisson_lpmf(y | mu);
  return com_poisson_log_lpmf(y | log(mu), nu, unique, logfac);
}
}
data {
  int<lower=0> n;  // Sample size
  int<lower=0> nc;  // Number of countries
  int<lower=0> y[n];  // Observed count
  int<lower=1> cid[n];  // Country index
  vector[11] unique; // Unique values for count
  vector[11] logfac; // Log-factorial for each unique count
  real prior_lam_m[nc]; // Prior mean for lambda 
  real<lower=0> prior_lam_s[nc]; // Prior sd for lambda
  real prior_nu_m[nc]; // Prior mean for nu (on log-scale) 
  real<lower=0> prior_nu_s[nc]; // Prior sd for nu
}
parameters {
  real lambda[nc]; 
  real nu[nc]; 
}
model {
  lambda ~ normal(prior_lam_m, prior_lam_s);
  nu ~ normal(0, prior_nu_s);
  for (i in 1:n){
    y[i] ~ com_poisson_log(lambda[cid[i]], exp(nu[cid[i]]), unique, logfac);
  }
}
'

# Initialise priors (these are on log-scale)
prior_lam_m = rep(log(mean(as.matrix(all_years_hfd_counts[, 3:1002]))), 25)
prior_lam_s = rep(1, 25)
prior_nu_m = rep(0, 25)
prior_nu_s = rep(1, 25)
res_lam <- res_nu <- res_ll <- NULL

for(i in 53:65){
  ii <- all_years_hfd_counts$cohort == 1917 + i
  counts_m <- as.matrix(all_years_hfd_counts[ii, 3:1002])
  counts_vec <- as.numeric(counts_m)
    #sample(as.numeric(counts_m), 100*sum(ii), replace = T)
  # Get vector of country ids in given year
  # Subset country id variable
  cid <- rep(all_years_hfd_counts$countryID[ii], 1000)
  loop_data <- list(
    n = length(counts_vec),
    y = counts_vec,
    unique = 0:10,
    logfac = lfactorial(0:10),
    nc = 25,
    cid = cid,
    prior_lam_m = prior_lam_m,
    prior_lam_s = prior_lam_s,
    prior_nu_m = prior_nu_m,
    prior_nu_s = prior_nu_s
  )
  
  loop_fit <- stan(
    model_code = model_CMP_int_country,      # R program
    data = loop_data,    # named list of data
    chains = 1,             # number of Markov chains
    warmup = 1000,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = 1
  )
  
  # Need a way of subsetting results for both lambda and nu
  new_fit_mat <- as.matrix(loop_fit)
  prior_lam_m <- apply(new_fit_mat, 2, mean)[1:25]
  prior_lam_s <- rep(1, 25) #2*apply(new_fit_mat, 2, sd)[1:25]
  prior_nu_m <- apply(new_fit_mat, 2, mean)[26:50]
  prior_nu_s <- rep(1, 25) #2*apply(new_fit_mat, 2, sd)[26:50]
  
  # Store samples for both lambda and nu - needs to be updated
  res_lam <- rbind(res_lam, new_fit_mat[, 1:25]) 
  res_nu <- rbind(res_nu, new_fit_mat[, 26:50]) 
  res_ll <- rbind(res_ll, new_fit_mat[, 51])
}

# Convert to mu
# Note: lambda and nu both on log-scale in above loop
res_lam_store <- exp(res_lam)
res_nu_store <- exp(res_nu)
# Use the below as a basis
G <- rowSums(exp(tcrossprod(res_lam[, 1], 0:10) - 
                        tcrossprod(exp(res_nu[, 1]), lfactorial(0:10))))
mu <- colSums(t(exp(tcrossprod(res_lam[, 1], 0:10) - 
                    tcrossprod(exp(res_nu[, 1]), lfactorial(0:10))))*(0:10))/G
# Loop through all countries
res_mu <- NULL
for (i in 1:25){
  G <- rowSums(exp(tcrossprod(res_lam[, i], 0:10) - 
                     tcrossprod(exp(res_nu[, i]), lfactorial(0:10))))
  mu <- colSums(t(exp(tcrossprod(res_lam[, i], 0:10) - 
                        tcrossprod(exp(res_nu[, i]), 
                                   lfactorial(0:10))))*(0:10))/G
  res_mu <- cbind(res_mu, mu)
}
res_sig <- NULL
# Need res_lam to contain log_lam values and res_nu to contain log_nu values
for (i in 1:25){
  G <- rowSums(exp(tcrossprod(res_llam[, i], 0:10) - 
                     tcrossprod(exp(res_lnu[, i]), lfactorial(0:10))))
  mu <- colSums(t(exp(tcrossprod(res_llam[, i], 0:10) - 
                        tcrossprod(exp(res_lnu[, i]), 
                                   lfactorial(0:10))))*(0:10))/G
  ex2 <- colSums(t(exp(tcrossprod(res_llam[, i], 0:10) - 
                         tcrossprod(exp(res_lnu[, i]), 
                                    lfactorial(0:10))))*(0:10)^2)/G
  sig <- ex2 - mu^2
  res_sig <- cbind(res_sig, sig)
}


# Save results
save(res_mu, file = "res_mu.RData")
save(res_nu, file = "res_nu.RData")
save(res_lam, file = "res_lam.RData")
save(res_ll, file = "res_ll.RData")
save(res_sig, file = "res_sig.RData")

# Some plotting
muUSA <- res_mu[1:64000, 25]
library(coda)
df_USA_res <- data.frame(mu = muUSA, year = rep(1918:1981, each = 1000))
# Boxplot
USA_plot <- ggplot(data = df_USA_res, aes(x = year, y = mu)) + 
  ylab("Average parity") + xlab("Year") + geom_boxplot(aes(group = year)) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", fill = "white")) +
  theme(legend.key = element_rect(fill = "transparent")) + 
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12)) 
# Ribbon plot
muUSAmat <- matrix(muUSA, 1000, 64, byrow = F)
bounds_mu <- HPDinterval(as.mcmc(muUSAmat))
lower_mu <- bounds_mu[, 1]
upper_mu <- bounds_mu[, 2]
df_USA_res2 <- data.frame(mu = colMeans(muUSAmat), muL = lower_mu, muU = upper_mu, 
                          year = 1918:1981)
USA_plot <- ggplot(data = df_USA_res2, aes(x = year, y = mu)) + 
  ylab("Parity") + xlab("Year") + 
  geom_ribbon(data = df_USA_res2, aes(x = year, ymin = muL, ymax = muU),
                                      , fill = "gray60") + 
  geom_line(aes(x = year, y = mu)) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", fill = "white")) +
  theme(legend.key = element_rect(fill = "transparent")) + 
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12)) 
# Repeat for variance (or SD)
sigUSA <- res_sig[1:64000, 25]
sigUSAmat <- matrix(sigUSA, 1000, 64, byrow = F)
bounds_sig <- HPDinterval(as.mcmc(sigUSAmat))
lower_sig <- bounds_sig[, 1]
upper_sig <- bounds_sig[, 2]
df_USA_res3 <- data.frame(sig = colMeans(sigUSAmat), sigL = lower_sig, 
                          sigU = upper_sig, year = 1918:1981)
USA_sig_plot <- ggplot(data = df_USA_res3, aes(x = year, y = sig)) + 
  ylab("Parity") + xlab("Year") + 
  geom_ribbon(data = df_USA_res3, aes(x = year, ymin = sigL, ymax = sigU),
              , fill = "gray60") + 
  geom_line(aes(x = year, y = sig)) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", fill = "white")) +
  theme(legend.key = element_rect(fill = "transparent")) + 
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12)) 
