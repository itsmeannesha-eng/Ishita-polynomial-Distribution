# 1. Define Parameters
theta <- 0.5   
m <- 3         

# 2. Define Coefficients C_k 
C <- rep(0, m + 1)
C[1] <- theta      
C[2] <- 0          
C[3] <- 1          
C[4] <- 1          

# 3. Construct Matrix A 
A <- matrix(0, nrow = m + 1, ncol = m + 1)

for (i in 0:m) {
  for (j in 0:m) {
    #  formula: n! / a^(n+1) where a = 2*theta
    A[i + 1, j + 1] <- factorial(i + j) / (2 * theta)^(i + j + 1)
  }
}

# 4. Construct Vector b
b <- rep(0, m + 1)

for (i in 0:m) {
  sum_val <- 0
  # We sum over k from 0 to m to account for all components of f(x)
  for (k in 0:m) {
    multiplier <- if(k == 0) theta else 1
    
    term <- C[k + 1] * multiplier * factorial(i + k) / (2 * theta)^(i + k + 1)
    sum_val <- sum_val + term
  }
  b[i + 1] <- sum_val
}

# 5. Solve for Beta (beta = A^-1 * b)
beta_vector <- solve(A, b)

# 6. Display Results
cat("Matrix A:\n")
print(A)
cat("\nVector b:\n")
print(b)
cat("\nResulting Beta Coefficients:\n")
names(beta_vector) <- paste0("beta_", 0:m)
print(beta_vector)

# Check if A * beta = b
# print(A %*% beta_vector)

