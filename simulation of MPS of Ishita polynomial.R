#MPS SIMULATION FOR ISHITA–POLYNOMIAL DISTRIBUTION


rm(list = ls())

library(stats)


# NORMALIZING CONSTANT

C_m <- function(theta, m){
  
  part1 <- (theta^3 + 2) / theta^3
  
  if(m < 3)
    return(1 / part1)
  
  k <- 3:m
  
  sum_k <- sum(factorial(k) * theta^(m - k))
  
  denom <- part1 +
    (1 / theta^(m + 1)) * sum_k
  
  return(1 / denom)
}


#  CDF

F_m <- function(x, theta, m){
  
  Cm <- C_m(theta, m)
  
  ex <- exp(-theta * x)
  
  # Base part
  part1 <- (1 - ex) +
    (2 / theta^3) * (1 - ex) -
    (x^2 * ex) / theta -
    (2 * x * ex) / theta^2
  
  # Polynomial part
  sum_poly <- 0
  
  if(m >= 3){
    
    for(k in 3:m){
      
      j <- 0:k
      
      inner <- sum(
        factorial(k) /
          factorial(j) *
          theta^(j - k - 1) *
          x^j * ex
      )
      
      sum_poly <- sum_poly + inner
    }
  }
  
  cdf <- (Cm * part1) + 1 - (Cm * sum_poly)
  
  return(pmin(pmax(cdf, 1e-10), 1 - 1e-10))
}


# VECTORISED CDF

pIP_fast <- function(x, theta, m){
  
  sapply(x, F_m,
         theta = theta,
         m = m)
}


# LOG-SPACING FUNCTION

log_spacing_fn <- function(theta, x_ord, m){
  
  if(theta <= 0)
    return(-1e12)
  
  F_vals <- pIP_fast(x_ord, theta, m)
  
  D <- diff(c(0, F_vals, 1))
  
  D[D <= 1e-12] <- 1e-12
  
  sum(log(D))
}


#RANDOM GENERATION

generate_data <- function(n, theta, m){
  
  u <- runif(n)
  
  x <- numeric(n)
  
  for(i in 1:n){
    
    root <- try(
      uniroot(
        function(v)
          F_m(v, theta, m) - u[i],
        
        interval = c(0, 50)
      )$root,
      
      silent = TRUE
    )
    
    if(inherits(root, "try-error")){
      x[i] <- NA
    } else {
      x[i] <- root
    }
  }
  
  sort(na.omit(x))
}

#n,m,theta,N values

n_vals     <- c(50, 100, 200, 500)

m_vals     <- c(2, 3, 4, 5)

theta_vals <- c(1.0, 1.5, 2.0, 2.5, 3.0)

N_reps     <- 1000


# RESULT

results <- expand.grid(
  
  n = n_vals,
  
  m = m_vals,
  
  theta_true = theta_vals
)

results$Mean_Est <- NA

results$Bias <- NA

results$MSE <- NA

results$RMSE <- NA


# SIMULATION 

cat("Starting Fast Simulation...\n")

start_time <- Sys.time()

#LOOP

for(i in 1:nrow(results)){
  
  n_i <- results$n[i]
  
  m_i <- results$m[i]
  
  t_i <- results$theta_true[i]
  
  estimates <- numeric(N_reps)
  
  for(r in 1:N_reps){
    
    # Generate sample
    samp <- generate_data(n_i, t_i, m_i)
    
    # Skip failed samples
    if(length(samp) < n_i/2){
      estimates[r] <- NA
      next
    }
    
    # MPS Optimization
    opt <- optimize(
      
      f = log_spacing_fn,
      
      interval = c(0.05, 5),
      
      maximum = TRUE,
      
      x_ord = samp,
      
      m = m_i
    )
    
    estimates[r] <- opt$maximum
  }
  
  # Remove NA estimates
  estimates <- na.omit(estimates)
  
  # Performance Measures
  mean_est <- mean(estimates)
  
  bias_est <- mean_est - t_i
  
  mse_est <- mean((estimates - t_i)^2)
  
  rmse_est <- sqrt(mse_est)
  
  # Store
  results$Mean_Est[i] <- round(mean_est, 4)
  
  results$Bias[i] <- round(bias_est, 4)
  
  results$MSE[i] <- round(mse_est, 6)
  
  results$RMSE[i] <- round(rmse_est, 6)
  
  cat(
    "Completed:",
    i, "/",
    nrow(results),
    "\n"
  )
}


end_time <- Sys.time()

cat("\nTotal Time:\n")

print(end_time - start_time)


# FINAL RESULTS

cat("\nSimulation Results\n")

print(results)



