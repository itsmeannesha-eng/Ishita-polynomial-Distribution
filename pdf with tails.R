rm(list = ls())


# Parameter

theta <- 2.5   

# x values
x <- seq(0, 20, length = 1000)
x_tail <- seq(5, 20, length = 1000)

# C_m(theta)
C_m <- function(m, theta){
  
  if(m == 1){
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

# Ishita polynomial pdf
f_m <- function(x, m, theta){
  
  poly_part <- theta + x^2
  
  if(m >= 3){
    for(k in 3:m){
      poly_part <- poly_part + x^k
    }
  }
  
  return(C_m(m, theta) * poly_part * exp(-theta * x))
}

# Distribution plot
plot(x, f_m(x,1,theta), type="l", lwd=2, col=1,
     ylim=c(0, max(f_m(x,2,theta))),
     ylab = expression(f[m](x,~theta)),
     xlab="x",
     main = bquote("Ishita-Polynomial Distribution ("*theta == .(theta)*")"))

for(m in 2:8){
  lines(x, f_m(x,m,theta), lwd=2, col=m)
}

legend("topright",
       legend = paste("m =", 1:8),
       col = 1:8,
       lwd = 2,
       cex = 0.6)

# tail grid plot
par(mfrow = c(2,4), mar = c(3,3,2,1))

for(m in 1:8){
  plot(x_tail, f_m(x_tail, m, theta),
       type = "l",
       lwd = 2,
       col = m,
       main = paste("Tail (m =", m, ")"),
       xlab = "x",
       ylab = "f(x)",
       cex.main = 0.9)
}

# Reset layout
par(mfrow = c(1,1))

rm(list = ls())


# Parameter

theta <- 0.36   

# x values
x <- seq(0, 20, length = 1000)
x_tail <- seq(5, 20, length = 1000)

# C_m(theta)
C_m <- function(m, theta){
  
  if(m == 1){
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

# Ishita polynomial pdf
f_m <- function(x, m, theta){
  
  poly_part <- theta + x^2
  
  if(m >= 3){
    for(k in 3:m){
      poly_part <- poly_part + x^k
    }
  }
  
  return(C_m(m, theta) * poly_part * exp(-theta * x))
}

# Distribution plot
plot(x, f_m(x,1,theta), type="l", lwd=2, col=1,
     ylim=c(0, max(f_m(x,2,theta))),
     ylab = expression(f[m](x,~theta)),
     xlab="x",
     main = bquote("Ishita-Polynomial Distribution ("*theta == .(theta)*")"))

for(m in 2:8){
  lines(x, f_m(x,m,theta), lwd=2, col=m)
}

legend("topright",
       legend = paste("m =", 1:8),
       col = 1:8,
       lwd = 2,
       cex = 0.6)

# tail grid plot
par(mfrow = c(2,4), mar = c(3,3,2,1))

for(m in 1:8){
  plot(x_tail, f_m(x_tail, m, theta),
       type = "l",
       lwd = 2,
       col = m,
       main = paste("Tail (m =", m, ")"),
       xlab = "x",
       ylab = "f(x)",
       cex.main = 0.9)
}

# Reset layout
par(mfrow = c(1,1))

