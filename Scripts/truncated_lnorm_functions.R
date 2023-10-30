require(tidyverse)
require(cubature)


# Fit truncated lognormal distribution ----------

# Following method of
# Pueyo S (2014) Algorithm for the maximum likelihood estimation of the parameters
# of the truncated normal and lognormal distributions 23–6 http://arxiv.org/abs/1407.6518

# User provides:
# --- patch_dist: patch size distribution as vector
# --- xmin: lower truncation value
# --- xmax: upper truncation value
# --- eta: convergence parameter; 0 < eta < 1

fit_truncated_ln <- function(patch_dist, xmin, xmax = Inf, eta = 0.1){
  
  # Filter patch size distribution to only include values >= xmin
  patch_dist <- patch_dist[patch_dist >= xmin]
  
  # Log-transform data
  y <- log(patch_dist)
  ymin <- log(xmin)
  ymax <- log(xmax)
  
  # Calculate sampling means
  y_bar <- mean(y)
  y2_bar <- mean(y^2)
  y3_bar <- mean(y^3)
  y4_bar <- mean(y^4)
  
  # Calculate h, a, b, c
  h <- y4_bar*(-y2_bar + y_bar^2) + y3_bar*(y3_bar - 2*y_bar*y2_bar) + y2_bar^3
  a <- (y4_bar - y2_bar^2) / h
  b <- (-y3_bar + y_bar*y2_bar) / h
  c <- (y2_bar - y_bar^2) / h
  
  # Starting values for alpha and psi
  alpha_j <- 1
  psi_j <- 0.1
  
  # Starting values for deltas
  delta_j_alpha <- delta_j_psi <- 1
  
  ##### Implement algorithm
  while (abs(delta_j_alpha) > eta/100000 & abs(delta_j_psi) > eta/100000) {

    # Normalization constant
    f_con <- function(u) {exp(-alpha_j*u - psi_j*u^2)}
    constant <- cubintegrate(f_con, lower = ymin, upper = ymax, method = "pcubature")$integral
    
    # Expectation of y
    f_Ey <- function(y) {y * (exp(-alpha_j*y - psi_j*y^2) / constant)}
    Ey <- cubintegrate(f_Ey, lower = ymin, upper = ymax, method = "pcubature")$integral
    
    # Expectation of y^2
    f_Ey2 <- function(y) {y^2 * (exp(-alpha_j*y - psi_j*y^2) / constant)}
    Ey2 <- cubintegrate(f_Ey2, lower = ymin, upper = ymax, method = "pcubature")$integral
    
    # Calculate deltas
    delta_j_alpha <- a*eta*(y_bar - Ey) + b*eta*(y2_bar - Ey2)
    delta_j_psi <-   b*eta*(y_bar - Ey) + c*eta*(y2_bar - Ey2)
    
    # Update of parameters
    alpha_j <- alpha_j + delta_j_alpha
    psi_j <- psi_j + delta_j_psi
  }
  
  # Return final parameters
  return( list(beta = 1 + alpha_j, psi = psi_j) )

}


# Plot truncated lognormal distribution ----------

# User provides:
# --- patch_dist: patch size distribution as vector
# --- xmin: lower truncation value
# --- xmax: upper truncation value
# --- beta: estimated beta parameter
# --- psi: estimated psi parameter

plot_truncated_ln <- function(patch_dist, xmin, xmax, beta, psi){
  
  # Filter patch size distribution to only include values >= xmin
  patch_dist <- patch_dist[patch_dist >= xmin]
  
  # Calculate normalization constant, a
  f_a <- function(x) { exp( -beta*log(x) - psi*log(x)^2 ) }
  a <- 1 / cubintegrate(f = f_a, lower = xmin, upper = xmax, method = "pcubature")$integral
  
  # Specify fitted lognormal density function
  f_ln <- function(x) { a * exp( -beta*log(x) - psi*log(x)^2 ) }
  
  # Loop through and calculate inverse cumulative probabilities
  # i.e., P(X >= x)
  x_range <- exp( seq(from = log(xmin), to = log(max(patch_dist)), by = 0.1) )
  P_predict <- c()
  
  for (i in 1:length(x_range)) {
    P <- cubintegrate(f = f_ln, lower = x_range[i], upper = xmax, method = "pcubature")$integral
    P_predict <- c(P_predict, P)
  }
  
  # Store predictions in tibble for plotting
  ln_predict <- tibble(x = x_range, P_x = P_predict)

  # Calculate inverse empirical cumulative density
  x_n <- length(patch_dist)
  x_ecdf <- c()
  for(x in 1:x_n){
    x_ecdf <- c(x_ecdf,
                sum(patch_dist >= patch_dist[x]) / x_n)
  }
  dist <- tibble(patch_area_ha = patch_dist, inv_ecdf = x_ecdf)
  
  # Plot inverse empirical cumulative density and fitted curve
  ggplot(data = dist, aes(x = patch_area_ha, y = inv_ecdf)) +
    geom_point(shape = 1, alpha = 0.5) +
    geom_line(data = ln_predict, aes(x = x, y = P_x)) +
    scale_x_log10(name = "Patch size (ha)") +
    scale_y_log10(name = "Inverse CDF [P(X >= x)]") 
}


# Extract residuals from truncated ln fit to empirical data --------

# User provides:
# --- patch_dist: patch size distribution as vector
# --- xmin: lower truncation value
# --- xmax: upper truncation value
# --- beta: estimated beta parameter
# --- psi: estimated psi parameter

resid_truncated_ln <- function(patch_dist, xmin, xmax, beta, psi){
  
  # Filter patch size distribution to only include values >= xmin
  patch_dist <- patch_dist[patch_dist >= xmin]
  
  # Calculate normalization constant, a
  f_a <- function(x) { exp( -beta*log(x) - psi*log(x)^2 ) }
  a <- 1 / cubintegrate(f = f_a, lower = xmin, upper = xmax, method = "pcubature")$integral
  
  # Specify fitted lognormal density function
  f_ln <- function(x) { a * exp( -beta*log(x) - psi*log(x)^2 ) }
  
  # Loop through and calculate inverse cumulative probabilities for observed data
  # i.e., P(X >= x)
  P_predict <- c()
  
  for (i in 1:length(patch_dist)) {
    P <- cubintegrate(f = f_ln, lower = patch_dist[i], upper = xmax, method = "pcubature")$integral
    P_predict <- c(P_predict, P)
  }
  
  # Store predictions in tibble
  ln_resid <- tibble(patch_area_ha = patch_dist, P_x = P_predict)
  
  # Calculate inverse empirical cumulative densities
  x_n <- length(patch_dist)
  x_ecdf <- c()
  for(x in 1:x_n){
    x_ecdf <- c(x_ecdf,
                sum(patch_dist >= patch_dist[x]) / x_n)
  }
  
  # Add inv_ecdf to tibble and calculate residuals
  ln_resid <- ln_resid %>% 
    mutate(inv_ecdf = x_ecdf) %>%
    mutate(resid = inv_ecdf - P_x)
  
  # Return tibble
  return(ln_resid)
}


# Functions for simulating from truncated lognormal ---------
# that is parameterized by xmin, xmax, beta, psi

# Integration helper function
# Calculates integral of exp( -beta*log(x) - psi*log(x)^2 )
# between given bounds, lwr and upr
int_fn <- function(lwr, upr, beta, psi){

  upper = sqrt(pi) * exp( ((beta-1)^2) / (4*psi) ) *
    pracma::erf( (2*psi*log(upr) + beta - 1) / (2*sqrt(psi)) ) / (2*sqrt(psi))
  lower = sqrt(pi) * exp( ((beta-1)^2) / (4*psi) ) *
    pracma::erf( (2*psi*log(lwr) + beta - 1) / (2*sqrt(psi)) ) / (2*sqrt(psi))

  return(upper - lower)
}

# Probability density function
d_truncated_ln <- function(x, xmin, xmax, beta, psi){
  
  # Check that x (either scalar or vector) > xmin and < xmax
  if (sum(x < xmin) > 0)
    stop("'x' must be greater than 'xmin'")
  if (sum(x > xmax) > 0)
    stop("'x' must be less than 'xmax'")
  
  # Calculate normalization constant, a
  # When psi > 0.01, use closed-form integral, otherwise use numerical integration methods
  if(psi > 0.01){
    a <- 1 / int_fn(lwr = xmin, upr = xmax, beta, psi)
  } else {
    f_a <- function(x) { exp( -beta*log(x) - psi*log(x)^2 ) }
    a <- 1 / cubintegrate(f = f_a, lower = xmin, upper = xmax, method = "pcubature")$integral
  }
  
  # Specify probability density function
  f_pdf <- function(x) { a * exp( -beta*log(x) - psi*log(x)^2 ) }
  
  return(f_pdf(x))
}

# Cumulative density function
p_truncated_ln <- function(x, xmin, xmax, beta, psi){

  # Check that x (either scalar or vector) > xmin and < xmax
  if (sum(x < xmin) > 0)
    stop("'x' must be greater than 'xmin'")
  if (sum(x > xmax) > 0)
    stop("'x' must be less than 'xmax'")

  # When psi > 0.01, use closed-form integral, otherwise use numerical integration methods
  if(psi > 0.01){
    
    # Calculate normalization constant, a
    a <- 1 / int_fn(lwr = xmin, upr = xmax, beta, psi)
    
    # Specify cumulative density function
    f_cdf <- function(x) { a * int_fn(lwr = xmin, upr = x, beta, psi) }
    
    return(f_cdf(x))
    
  } else {
    
    # Calculate normalization constant, a
    f_a <- function(x) { exp( -beta*log(x) - psi*log(x)^2 ) }
    a <- 1 / cubintegrate(f = f_a, lower = xmin, upper = xmax, method = "pcubature")$integral
    
    # Specify probability density function
    f_pdf <- function(x) { a * exp( -beta*log(x) - psi*log(x)^2 ) }
    
    # Loop through and numerically integrate to calculate cumulative probabilities
    # i.e., P(X <= x)
    CDF <- c()
    for (i in 1:length(x)) {
      P <- cubintegrate(f = f_pdf, lower = xmin, upper = x[i], method = "pcubature")$integral
      CDF <- c(CDF, P)
    }
    
    return(CDF)
  }
}


# Quantile function
# This is an approximate numerical solution -- not an exact solution
q_truncated_ln <- function(p, xmin, xmax, beta, psi){

  # Check that p is between zero and one
  if (any(p < 0) || any(p > 1))
    stop("'p' must be between 0 and 1")

  q <- numeric(length(p))

  # Set 0 and 1 to xmin and xmax, respectively
  p.low <- p == 0
  q[p.low] <- xmin
  p.high <- p == 1
  q[p.high] <- xmax

  # Loop through and calculate remaining quantiles using uniroot()
  # Approach is to subtract input 'p' from CDF and then find root
  # i.e., find value of x that makes function equal to zero

  # Root function -- CDF minus cumulative probability, p
  root_fn <- function(x, xmin, xmax, beta, psi, p){
    p_truncated_ln(x, xmin, xmax, beta, psi) - p
  }

  # Find roots for input 0 < p < 1
  if (any(index <- !(p.low | p.high))) {

    for (i in which(index)){
      q[i] <- uniroot(root_fn, c(xmin, xmax),
                      xmin = xmin, xmax = xmax,
                      beta = beta, psi = psi, p = p[i],
                      tol = 0.0001)$root
    }
  }
  return(q)
}


# Random sample generation
# Randomly selects cumulative probabilities using runif()
# and calculates associated quantiles
r_truncated_ln <- function (n, xmin, xmax, beta, psi){
  
  if (is.na(n) || n <= 0 || n != trunc(n) || length(n) > 1) 
    stop("'n' must be a positive integer")
  
  r <- q_truncated_ln(p = runif(n), xmin, xmax, beta, psi)
  return(r)
  
}


