# ISHITA-POLYNOMIAL DISTRIBUTION
# COMPLETE REAL DATA ANALYSIS USING MOTOR FAILURE DATA
rm(list = ls())

library(boot)
library(ggplot2)
library(moments)


data(motor)
data <- motor$times
data <- data[data > 0]
n <- length(data)


Cm <- function(theta,m){
  
  val <- (theta^3+2)/theta^3 +
    sum(sapply(3:m,function(k)
      factorial(k)/theta^(k+1)))
  
  1/val
}


dIP <- function(x,theta,m){
  
  C <- Cm(theta,m)
  
  poly <- theta + x^2 +
    rowSums(sapply(3:m,function(k)x^k))
  
  C*poly*exp(-theta*x)
}


pIP <- function(q,theta,m){
  
  sapply(q,function(z)
    integrate(function(t)
      dIP(t,theta,m),
      0,z)$value)
}


hIP <- function(x,theta,m){
  
  dIP(x,theta,m)/(1-pIP(x,theta,m))
}

# INVERSE EXPONENTIAL DISTRIBUTION
dIE <- function(x,theta)
  (theta/x^2)*exp(-theta/x)


pIE <- function(x,theta)
  exp(-theta/x)

# TRANSMUTED ISHITA DISTRIBUTION
# Base: Standard Ishita (m = 2)

dTransIshita <- function(x,theta,lambda){
  
  f <- dIP(x,theta,2)
  
  F_val <- pIP(x,theta,2)
  
  f*(1-lambda+2*lambda*F_val)
}


pTransIshita <- function(q,theta,lambda){
  
  F_val <- pIP(q,theta,2)
  
  F_val*(1+lambda-lambda*F_val)
}

# LENGTH-BIASED WEIGHTED ISHITA DISTRIBUTION
# Base: Standard Ishita (m = 2)
#
# f_LB(x) = x f(x) / E(X)

mean_Ishita <- function(theta){
  
  integrate(
    function(x)
      x*dIP(x,theta,2),
    0,
    Inf
  )$value
}


dLengthBiasedIshita <- function(x,theta){
  
  mu <- mean_Ishita(theta)
  
  x*dIP(x,theta,2)/mu
}


pLengthBiasedIshita <- function(q,theta){
  
  sapply(q,function(z)
    integrate(
      function(x)
        dLengthBiasedIshita(x,theta),
      0,
      z
    )$value)
}

# UNIT ISHITA DISTRIBUTION
#
# Unit Ishita is defined on 0 < u < 1.
#
# f_U(u;theta) =
# theta^3 u^(theta-1)
# [theta + log^2(u)]/(theta^3+2)
#
# F_U(u;theta) =
# u^theta [
# theta^3 + theta^2 log^2(u)
# - 2 theta log(u) + 2
# ]/(theta^3+2)
#
# For motor-failure data x > 0, use
#
# u = x/(1+x)
#
# with
#
# du/dx = 1/(1+x)^2.
#
# Thus the density on the original x-scale is
#
# f_X(x) =
# f_U(x/(1+x)) /(1+x)^2
dUnitIshita <- function(x,theta){
  
  ifelse(
    x <= 0,
    0,
    {
      
      u <- x/(1+x)
      
      theta^3*
        u^(theta-1)*
        (theta+log(u)^2)/
        (theta^3+2)*
        1/(1+x)^2
    }
  )
}


pUnitIshita <- function(x,theta){
  
  sapply(x,function(z){
    
    if(z <= 0)
      return(0)
    
    u <- z/(1+z)
    
    u^theta*
      (
        theta^3+
          theta^2*log(u)^2-
          2*theta*log(u)+
          2
      )/
      (theta^3+2)
  })
}

# DESCRIPTIVE STATISTICS
Table1 <- data.frame(
  
  Mean=mean(data),
  Median=median(data),
  SD=sd(data),
  Skewness=skewness(data),
  Kurtosis=kurtosis(data),
  Minimum=min(data),
  Maximum=max(data),
  Sample_Size=n
)


cat(
  "\nTable 1: Descriptive Statistics of Motor Failure Data\n"
)

print(Table1)

# MPS ESTIMATION
logSpacing <- function(theta,data,m){
  
  F <- pIP(
    sort(data),
    theta,
    m
  )
  
  F <- c(0,F,1)
  
  D <- diff(F)
  
  if(any(D<=0))
    return(-Inf)
  
  sum(log(D))
}


estimate_theta <- function(data,m){
  
  optimize(
    function(th)
      -logSpacing(th,data,m),
    c(.001,2)
  )$minimum
}

# ORDER SELECTION
mgrid <- 1:6

fit_table <- data.frame()


for(m in mgrid){
  
  th <- estimate_theta(
    data,
    m
  )
  
  fit_table <- rbind(
    
    fit_table,
    
    data.frame(
      m=m,
      theta=th,
      logSpacing=
        logSpacing(
          theta=th,
          data=data,
          m=m
        )
    )
  )
}


best <- fit_table[
  which.max(
    fit_table$logSpacing
  ),
]


theta_hat <- best$theta
m_hat <- best$m

# PARAMETER ESTIMATION
Table2 <- fit_table


colnames(Table2) <-
  c(
    "Polynomial_Order_m",
    "Estimated_Theta",
    "Log_Spacing_Value"
  )


cat(
  "\nTable 2: MPS Parameter Estimates for Ishita Polynomial\n"
)

print(Table2)


# INVERSE EXPONENTIAL PARAMETER ESTIMATION
lambda_exp <- 1/mean(data)


theta_IE <- optimize(
  
  function(th)
    -sum(
      log(
        pmax(
          dIE(data,th),
          1e-10
        )
      )
    ),
  
  c(.001,10)
  
)$minimum

# TRANSMUTED ISHITA PARAMETER ESTIMATION
fit_trans <- optim(
  
  par=c(
    theta=0.15,
    lambda=0.1
  ),
  
  fn=function(p){
    
    th <- p[1]
    lam <- p[2]
    
    if(
      th <= 0 ||
      lam < -1 ||
      lam > 1
    )
      return(1e10)
    
    dens <- dTransIshita(
      data,
      th,
      lam
    )
    
    if(any(dens <= 0))
      return(1e10)
    
    -sum(log(dens))
  }
)


theta_trans <- fit_trans$par[1]
lambda_trans <- fit_trans$par[2]

# LENGTH-BIASED WEIGHTED ISHITA PARAMETER ESTIMATION
fit_LBI <- optimize(
  
  function(th){
    
    if(th <= 0)
      return(1e10)
    
    dens <- dLengthBiasedIshita(
      data,
      th
    )
    
    if(any(
      !is.finite(dens) |
      dens <= 0
    ))
      return(1e10)
    
    -sum(log(dens))
  },
  
  interval=c(.001,10)
)


theta_LBI <- fit_LBI$minimum

# UNIT ISHITA PARAMETER ESTIMATION
fit_UI <- optimize(
  
  function(th){
    
    if(th <= 0)
      return(1e10)
    
    dens <- dUnitIshita(
      data,
      th
    )
    
    if(any(
      !is.finite(dens) |
      dens <= 0
    ))
      return(1e10)
    
    -sum(log(dens))
  },
  
  interval=c(.001,10)
)


theta_UI <- fit_UI$minimum

# MULTIPLE m ANALYSIS (2 TO 5)

mgrid <- 2:5

fit_all <- data.frame()


for(m in mgrid){
  
  th <- estimate_theta(
    data,
    m
  )
  
  fit_all <- rbind(
    
    fit_all,
    
    data.frame(
      m=m,
      theta=th
    )
  )
}


print(fit_all)


# PREPARE PLOT DATA
xx <- seq(
  min(data),
  max(data),
  length=400
)


plot_data <- data.frame()


for(i in 1:nrow(fit_all)){
  
  m_val <- fit_all$m[i]
  th_val <- fit_all$theta[i]
  
  temp <- data.frame(
    
    x=xx,
    
    m=factor(m_val),
    
    pdf=dIP(
      xx,
      th_val,
      m_val
    ),
    
    cdf=pIP(
      xx,
      th_val,
      m_val
    ),
    
    hazard=hIP(
      xx,
      th_val,
      m_val
    )
  )
  
  plot_data <- rbind(
    plot_data,
    temp
  )
}


# PDF COMPARISON

ggplot(
  plot_data,
  aes(x,pdf,color=m)
)+
  
  geom_line()+
  
  labs(title="")



# CDF COMPARISON
ggplot(
  plot_data,
  aes(x,cdf,color=m)
)+
  
  geom_line()+
  
  labs(title="")


# HAZARD COMPARISON
ggplot(
  plot_data,
  aes(x,hazard,color=m)
)+
  
  geom_line()+
  
  labs(title="")


# HISTOGRAM + ALL ISHITA-POLYNOMIAL PDFs
ggplot(
  data.frame(data),
  aes(data)
)+
  
  geom_histogram(
    aes(y=..density..),
    bins=30,
    fill="grey70"
  )+
  
  geom_line(
    data=plot_data,
    aes(
      x=x,
      y=pdf,
      color=m
    )
  )


# ORDER SELECTION

ggplot(
  fit_table,
  aes(
    m,
    logSpacing
  )
)+
  
  geom_line(
    color="brown"
  )+
  
  geom_point(
    size=3
  )+
  
  labs(title="")

# MODEL COMPARISON
#
# Models included:
# 1. Ishita Polynomial
# 2. Transmuted Ishita
# 3. Inverse Exponential
# 4. Length-Biased Weighted Ishita
# 5. Unit Ishita
#
# Power Ishita and Beta-Exponentiated Ishita
# have been removed.
ggplot(
  data.frame(data),
  aes(x=data)
)+
  
  geom_histogram(
    aes(
      y=..density..,
      fill="Observed Data"
    ),
    bins=30,
    alpha=.6,
    color="black"
  )+
  
  
# Ishita Polynomial

stat_function(
  aes(
    color="Ishita Polynomial",
    linetype="Ishita Polynomial"
  ),
  
  fun=function(z)
    dIP(
      z,
      theta_hat,
      m_hat
    ),
  
  linewidth=1.4
)+
  
  
# Transmuted Ishita
  
stat_function(
  aes(
    color="Transmuted Ishita",
    linetype="Transmuted Ishita"
  ),
  
  fun=function(z)
    dTransIshita(
      z,
      theta_trans,
      lambda_trans
    ),
  
  linewidth=1.2
)+
  

# Inverse Exponential
stat_function(
  aes(
    color="Inverse Exponential",
    linetype="Inverse Exponential"
  ),
  
  fun=function(z)
    dIE(
      z,
      theta_IE
    ),
  
  linewidth=1.2
)+
  
  
# Length-Biased Weighted Ishita
stat_function(
  aes(
    color="Length-Biased Weighted Ishita",
    linetype="Length-Biased Weighted Ishita"
  ),
  
  fun=function(z)
    dLengthBiasedIshita(
      z,
      theta_LBI
    ),
  
  linewidth=1.2
)+
  

# Unit Ishita
  
stat_function(
  aes(
    color="Unit Ishita",
    linetype="Unit Ishita"
  ),
  
  fun=function(z)
    dUnitIshita(
      z,
      theta_UI
    ),
  
  linewidth=1.2
)+
  
 
# COLOUR SCALE

scale_color_manual(
  
  name="Fitted Models",
  
  values=c(
    
    "Ishita Polynomial"="red",
    
    "Transmuted Ishita"="darkgreen",
    
    "Inverse Exponential"="orange",
    
    "Length-Biased Weighted Ishita"="blue",
    
    "Unit Ishita"="brown"
  )
)+
  

# LINE TYPE SCALE

scale_linetype_manual(
  
  name="Fitted Models",
  
  values=c(
    
    "Ishita Polynomial"="solid",
    
    "Transmuted Ishita"="dotted",
    
    "Inverse Exponential"="solid",
    
    "Length-Biased Weighted Ishita"="dashed",
    
    "Unit Ishita"="dotdash"
  )
)+
  
  
  scale_fill_manual(
    
    values=c(
      "Observed Data"="grey80"
    )
  )+
  
  
  labs(
    title="",
    x="Failure Time",
    y="Density"
  )+
  
  theme_minimal()


# BOXPLOT
ggplot(
  data.frame(
    Model="Ishita Polynomial",
    data=data
  ),
  aes(
    x=Model,
    y=data
  )
)+
  
  geom_boxplot(
    fill="skyblue",
    color="black",
    width=.4,
    alpha=.7
  )+
  
  stat_summary(
    fun=mean,
    geom="point",
    shape=18,
    size=4,
    color="red"
  )+
  
  labs(
    title=
      "Boxplot under Ishita Polynomial Model",
    x="",
    y="Failure Time"
  )+
  
  theme_minimal(
    base_size=14
  )

# GOODNESS OF FIT COMPARISON
#
# Models included:
# 1. Ishita Polynomial
# 2. Transmuted Ishita
# 3. Inverse Exponential
# 4. Length-Biased Weighted Ishita
# 5. Unit Ishita

# LOG-LIKELIHOODS


# Ishita Polynomial
logLik_IP <- sum(
  log(
    pmax(
      dIP(
        data,
        theta_hat,
        m_hat
      ),
      1e-10
    )
  )
)

# Transmuted Ishita
logLik_Trans <- sum(
  log(
    pmax(
      dTransIshita(
        data,
        theta_trans,
        lambda_trans
      ),
      1e-10
    )
  )
)


# Inverse Exponential
logLik_IE <- sum(
  log(
    pmax(
      dIE(
        data,
        theta_IE
      ),
      1e-10
    )
  )
)


# Length-Biased Weighted Ishita
logLik_LBI <- sum(
  log(
    pmax(
      dLengthBiasedIshita(
        data,
        theta_LBI
      ),
      1e-10
    )
  )
)

# Unit Ishita

logLik_UI <- sum(
  log(
    pmax(
      dUnitIshita(
        data,
        theta_UI
      ),
      1e-10
    )
  )
)


# AIC AND BIC FUNCTIONS

AIC <- function(k,ll)
  -2*ll + 2*k


BIC <- function(k,ll,n)
  -2*ll + k*log(n)



# KOLMOGOROV-SMIRNOV TESTS

# Ishita Polynomial

KS_IP <- ks.test(
  data,
  function(q)
    pIP(
      q,
      theta_hat,
      m_hat
    )
)


# Transmuted Ishita

KS_Trans <- ks.test(
  data,
  function(q)
    pTransIshita(
      q,
      theta_trans,
      lambda_trans
    )
)


# Inverse Exponential

KS_IE <- ks.test(
  data,
  function(q)
    pIE(
      q,
      theta_IE
    )
)


# Length-Biased Weighted Ishita

KS_LBI <- ks.test(
  data,
  function(q)
    pLengthBiasedIshita(
      q,
      theta_LBI
    )
)


# Unit Ishita

KS_UI <- ks.test(
  data,
  function(q)
    pUnitIshita(
      q,
      theta_UI
    )
)


# FINAL GOODNESS-OF-FIT TABLE

Table3 <- data.frame(
  
  Model=c(
    
    "Ishita Polynomial",
    
    "Transmuted Ishita",
    
    "Inverse Exponential",
    
    "Length-Biased Weighted Ishita",
    
    "Unit Ishita"
  ),
  
  
  AIC=c(
    
    AIC(
      2,
      logLik_IP
    ),
    
    AIC(
      2,
      logLik_Trans
    ),
    
    AIC(
      1,
      logLik_IE
    ),
    
    AIC(
      1,
      logLik_LBI
    ),
    
    AIC(
      1,
      logLik_UI
    )
  ),
  
  
  BIC=c(
    
    BIC(
      2,
      logLik_IP,
      n
    ),
    
    BIC(
      2,
      logLik_Trans,
      n
    ),
    
    BIC(
      1,
      logLik_IE,
      n
    ),
    
    BIC(
      1,
      logLik_LBI,
      n
    ),
    
    BIC(
      1,
      logLik_UI,
      n
    )
  ),
  
  
  KS_Statistic=c(
    
    KS_IP$statistic,
    
    KS_Trans$statistic,
    
    KS_IE$statistic,
    
    KS_LBI$statistic,
    
    KS_UI$statistic
  ),
  
  
  KS_pvalue=c(
    
    KS_IP$p.value,
    
    KS_Trans$p.value,
    
    KS_IE$p.value,
    
    KS_LBI$p.value,
    
    KS_UI$p.value
  )
)

# PRINT FINAL GOODNESS-OF-FIT TABLE

cat(
  "\nTable 3: Goodness-of-Fit Comparison\n"
)

print(Table3)


# PRINT ESTIMATED PARAMETERS

cat(
  "\nTransmuted Ishita Parameters:\n"
)

cat(
  "Theta = ",
  theta_trans,
  "\n"
)

cat(
  "Lambda = ",
  lambda_trans,
  "\n"
)


cat(
  "\nLength-Biased Weighted Ishita Parameters:\n"
)

cat(
  "Theta = ",
  theta_LBI,
  "\n"
)


cat(
  "\nUnit Ishita Parameters:\n"
)

cat(
  "Theta = ",
  theta_UI,
  "\n"
)

