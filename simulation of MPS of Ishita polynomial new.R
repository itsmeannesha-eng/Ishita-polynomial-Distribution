# MPS SIMULATION FOR ISHITA-POLYNOMIAL DISTRIBUTION
rm(list = ls())
library(stats)
# NORMALIZING CONSTANT
C_m <- function(theta, m) {
  
  part1 <- (theta^3 + 2) / theta^3
  
  if (m < 3)
    return(1 / part1)
  
  k <- 3:m
  sum_k <- sum(factorial(k) * theta^(m - k))
  
  denom <- part1 + (1 / theta^(m + 1)) * sum_k
  
  return(1 / denom)
}

# CDF
F_m <- function(x, theta, m) {
  
  Cm <- C_m(theta, m)
  ex <- exp(-theta * x)
  
  # Base part
  part1 <- (1 - ex) +
    (2 / theta^3) * (1 - ex) -
    (x^2 * ex) / theta -
    (2 * x * ex) / theta^2
  
  # Polynomial part
  sum_poly <- 0
  
  if (m >= 3) {
    for (k in 3:m) {
      j <- 0:k
      
      inner <- sum(
        factorial(k) / factorial(j) *
          theta^(j - k - 1) *
          x^j * ex
      )
      
      sum_poly <- sum_poly + inner
    }
  }
  
  cdf <- (Cm * part1) + 1 - (Cm * sum_poly)
  
  return(pmin(pmax(cdf, 1e-10), 1 - 1e-10))
}

# VECTORIZED CDF
pIP_fast <- function(x, theta, m) {
  sapply(x, F_m, theta = theta, m = m)
}

# LOG-SPACING FUNCTION
log_spacing_fn <- function(theta, x_ord, m) {
  
  if (theta <= 0)
    return(-1e12)
  
  F_vals <- pIP_fast(x_ord, theta, m)
  D <- diff(c(0, F_vals, 1))
  D[D <= 1e-12] <- 1e-12
  
  return(sum(log(D)))
}

# RANDOM GENERATION
generate_data <- function(n, theta, m) {
  
  u <- runif(n)
  x <- numeric(n)
  
  for (i in 1:n) {
    
    root <- try(
      uniroot(
        function(v) F_m(v, theta, m) - u[i],
        interval = c(0, 50)
      )$root,
      silent = TRUE
    )
    
    if (inherits(root, "try-error")) {
      x[i] <- NA
    } else {
      x[i] <- root
    }
  }
  
  sort(na.omit(x))
}

# PERFORMANCE MEASURES
calculate_metrics <- function(estimates, true_theta) {
  
  errors <- estimates - true_theta
  
  c(
    Mean_Est = mean(estimates),
    Bias     = mean(errors),
    Variance = var(estimates),
    MSE      = mean(errors^2),
    AAD      = mean(abs(errors)),
    RMSE     = sqrt(mean(errors^2)),
    SD       = sd(estimates)
  )
}

# SIMULATION SETTINGS
n_vals <- c(50, 100, 200, 500, 1000, 2000, 5000)
m_vals <- c(2, 3, 4, 5)
theta_vals <- c(1.0, 1.5, 2.0, 2.5, 3.0)
N_reps <- 1000

# RESULT STRUCTURE
results <- expand.grid(
  n = n_vals,
  m = m_vals,
  theta_true = theta_vals
)

results$Mean_Est <- NA_real_
results$Bias <- NA_real_
results$Variance <- NA_real_
results$MSE <- NA_real_
results$AAD <- NA_real_
results$RMSE <- NA_real_
results$SD <- NA_real_
results$Successful_Reps <- NA_integer_

# SIMULATION
cat("Starting MPS Simulation...\n")
start_time <- Sys.time()

for (i in 1:nrow(results)) {
  
  n_i <- results$n[i]
  m_i <- results$m[i]
  t_i <- results$theta_true[i]
  
  estimates <- rep(NA_real_, N_reps)
  
  for (r in 1:N_reps) {
    
    # Generate sample
    samp <- generate_data(n_i, t_i, m_i)
    
    # Skip failed samples
    if (length(samp) < n_i / 2) {
      next
    }
    
    # MPS Optimization
    opt <- try(
      optimize(
        f = log_spacing_fn,
        interval = c(0.05, 5),
        maximum = TRUE,
        x_ord = samp,
        m = m_i
      ),
      silent = TRUE
    )
    
    if (!inherits(opt, "try-error") &&
        is.finite(opt$maximum)) {
      estimates[r] <- opt$maximum
    }
  }
  
  # Remove failed estimates
  estimates <- estimates[is.finite(estimates)]
  
  # Store successful replication count
  results$Successful_Reps[i] <- length(estimates)
  
  # Calculate performance measures
  if (length(estimates) > 1) {
    
    metrics <- calculate_metrics(estimates, t_i)
    
    results$Mean_Est[i] <- round(metrics["Mean_Est"], 4)
    results$Bias[i] <- round(metrics["Bias"], 4)
    results$Variance[i] <- round(metrics["Variance"], 6)
    results$MSE[i] <- round(metrics["MSE"], 6)
    results$AAD[i] <- round(metrics["AAD"], 6)
    results$RMSE[i] <- round(metrics["RMSE"], 6)
    results$SD[i] <- round(metrics["SD"], 6)
  }
  
  cat(
    "Completed:", i, "/", nrow(results),
    "| m =", m_i,
    "| theta =", t_i,
    "| n =", n_i,
    "| Successful =", length(estimates),
    "\n"
  )
}

# COMPUTATION TIME
end_time <- Sys.time()

cat("\nTotal Simulation Time:\n")
print(end_time - start_time)

# FINAL RESULTS
cat("\nMPS Simulation Results:\n")
print(results)

# SAVE RESULTS
write.csv(
  results,
  "MPS_Ishita_Polynomial_Results.csv",
  row.names = FALSE
)

cat("\nResults saved to MPS_Ishita_Polynomial_Results.csv\n")