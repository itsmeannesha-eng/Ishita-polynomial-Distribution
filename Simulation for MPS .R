
# MPS Estimation for Ishita-Polynomial Distribution
# Simulation 

# C_m(theta)
C_m <- function(theta, m) {
  part1 <- (theta^3 + 2) / theta^3
  if (m < 3) return(1 / part1)
  
  sum_k <- 0
  for (k in 3:m) {
    sum_k <- sum_k + (factorial(k) * (theta^(m - k)))
  }
  denom <- part1 + (1 / theta^(m + 1)) * sum_k
  return(1 / denom)
}

# 2.CDF
F_m <- function(x, theta, m) {
  if (x <= 0) return(0)
  Cm <- C_m(theta, m)
  
  # Base Ishita-like components
  part1 <- (1 - exp(-theta * x)) +
    (2 / theta^3) * (1 - exp(-theta * x)) -
    (x^2 * exp(-theta * x)) / theta -
    (2 * x * exp(-theta * x)) / theta^2
  
  # Polynomial Summation 
  sum_poly <- 0
  if (m >= 3) {
    for (k in 3:m) {
      inner_j <- 0
      for (j in 0:k) {
        inner_j <- inner_j + (factorial(k) / factorial(j)) * (theta^(j - k - 1)) * (x^j) * exp(-theta * x)
      }
      sum_poly <- sum_poly + inner_j
    }
  }
  
  cdf <- (Cm * part1) + 1 - (Cm * sum_poly)
  return(min(max(cdf, 0), 1)) 
}

# 3. Log-Spacing
log_spacing_fn <- function(theta, x_ord, m) {
  if (theta <= 0) return(-1e10) 
  
  #  F(x) 
  F_vals <- sapply(x_ord, function(v) F_m(v, theta, m))
  
  # Di(theta) = F(xi) - F(xi-1)
  # F(x0)=0 and F(xn+1)=1
  D <- diff(c(0, F_vals, 1))
  
  # replace zero or negative spacings with small value
  D[D <= 0] <- 1e-15
  
  
  return(sum(log(D)))
}

# Inverse Transform Sampling
generate_data <- function(n, theta, m) {
  u <- runif(n)
  x_gen <- sapply(u, function(ui) {
    res <- try(uniroot(function(v) F_m(v, theta, m) - ui,
                       interval = c(0, 100), extendInt = "yes"), silent = TRUE)
    if(inherits(res, "try-error")) return(NA) else return(res$root)
  })
  return(sort(na.omit(x_gen)))
}

#vzlues of n,m, theta 
n_vals     <- c(20, 50, 100, 200, 500)      
m_vals     <- c(3, 4, 5, 6, 7)              
theta_vals <- c(0.5, 1.0, 1.5, 2.0, 2.5)    
N_reps     <- 500                           

# Results
results <- expand.grid(n = n_vals, m = m_vals, theta_true = theta_vals)
results$Mean_Est <- 0
results$Bias     <- 0
results$MSE      <- 0

cat("Starting Simulation Grid...\n")

for (i in 1:nrow(results)) {
  n_i <- results$n[i]
  m_i <- results$m[i]
  t_i <- results$theta_true[i]
  
  estimates <- numeric(N_reps)
  
  for (r in 1:N_reps) {
    # 1. Generate data
    samp <- generate_data(n_i, t_i, m_i)
    
    # 2. Maximize Log-Spacing
    
    opt <- optimize(f = log_spacing_fn,
                    interval = c(0.01, 15),
                    maximum = TRUE,
                    x_ord = samp,
                    m = m_i)
    estimates[r] <- opt$maximum
  }
  
  # 3. Calculate metrics
  results$Mean_Est[i] <- mean(estimates)
  results$Bias[i]     <- mean(estimates) - t_i
  results$MSE[i]      <- mean((estimates - t_i)^2)
  
  
  if(i %% 10 == 0) cat(sprintf("Progress: %d/%d combinations complete\n", i, nrow(results)))
}

# Final Table
print(results)

