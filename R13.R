# parameter
theta <- 0.36

#  x values
x <- seq(0, 20, length = 1000)

#  C_m(theta)
C_m <- function(m, theta){
  
  if(m == 1){
    # For m=1 (only theta term)
    val <- 1
  } else if(m == 2){
    val <- 1 + 2/theta^3
  } else {
    val <- 1 + 2/theta^3
    for(k in 3:m){
      val <- val + factorial(k)/theta^(k+1)
    }
  }
  
  return(1/val)
}

# Ishita-Polynomial pdf
f_m <- function(x, m, theta){
  
  poly_part <- theta + x^2
  
  if(m >= 3){
    for(k in 3:m){
      poly_part <- poly_part + x^k
    }
  }
  
  return(C_m(m, theta) * poly_part * exp(-theta * x))
}

# Plot setup
plot(x, f_m(x,1,theta), type="l", lwd=2, col=1,
     ylim=c(0, max(f_m(x,2,theta))),
     ylab = expression(f[m](x,~theta)),
     xlab="x",
     main="Ishita-Polynomial Distribution (theta = 0.36)")

#  curves for m = 2 to 8
for(m in 2:8){
  lines(x, f_m(x,m,theta), lwd=2, col=m)
}

legend("topright",
       legend=paste("m =", 1:8),
       col=1:8,
       lwd=2)

