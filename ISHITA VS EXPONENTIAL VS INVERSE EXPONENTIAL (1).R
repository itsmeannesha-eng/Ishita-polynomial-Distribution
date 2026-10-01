# ISHITA–POLYNOMIAL DISTRIBUTION
# COMPLETE REAL DATA ANALYSIS USING MOTOR FAILURE DATA

rm(list = ls())

library(boot)
library(ggplot2)
library(moments)


# real data

data(motor)
data <- motor$times
data <- data[data > 0]
n <- length(data)

#Ishita-polynomial distribution

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
      dIP(t,theta,m),0,z)$value)
}

hIP <- function(x,theta,m){
  dIP(x,theta,m)/(1-pIP(x,theta,m))
}

# Inverse exponential distribution

dIE <- function(x,theta)
  (theta/x^2)*exp(-theta/x)

pIE <- function(x,theta)
  exp(-theta/x)

# Transmuted Ishita distribution (Base: Standard Ishita m=2)
dTransIshita <- function(x, theta, lambda) {
  f <- dIP(x, theta, 2)
  F_val <- pIP(x, theta, 2)
  f * (1 - lambda + 2 * lambda * F_val)
}

pTransIshita <- function(q, theta, lambda) {
  F_val <- pIP(q, theta, 2)
  F_val * (1 + lambda - lambda * F_val)
}

# Power Ishita distribution (Base: Standard Ishita m=2)
dPowerIshita <- function(x, theta, alpha) {
  f <- dIP(x, theta, 2)
  F_val <- pIP(x, theta, 2)
  alpha * f * (F_val^(alpha - 1))
}

pPowerIshita <- function(q, theta, alpha) {
  F_val <- pIP(q, theta, 2)
  F_val^alpha
}


# Descriptive Statistics

Table1 <- data.frame(
  Mean=mean(data),
  Median=median(data),
  SD=sd(data),
  Skewness=skewness(data),
  Kurtosis=kurtosis(data),
  Minimum=min(data),
  Maximum=max(data),
  Sample_Size=n)

cat("\nTable 1: Descriptive Statistics of Motor Failure Data\n")
print(Table1)

# MPS estimation


logSpacing <- function(theta,data,m){
  
  F <- pIP(sort(data),theta,m)
  F <- c(0,F,1)
  
  D <- diff(F)
  if(any(D<=0)) return(-Inf)
  
  sum(log(D))
}

estimate_theta <- function(data,m){
  optimize(function(th)
    -logSpacing(th,data,m),
    c(.001,2))$minimum
}

#order selection

mgrid <- 1:6
fit_table <- data.frame()

for(m in mgrid){
  
  th <- estimate_theta(data,m)
  
  fit_table <- rbind(fit_table,
                     data.frame(m=m,
                                theta=th,
                                logSpacing=
                                  logSpacing(theta=th,data=data,m=m)))
}

best <- fit_table[which.max(fit_table$logSpacing),]

theta_hat <- best$theta
m_hat <- best$m

#parameter estimation

Table2 <- fit_table
colnames(Table2) <-
  c("Polynomial_Order_m",
    "Estimated_Theta",
    "Log_Spacing_Value")

cat("\nTable 2: MPS Parameter Estimates for Ishita Polynomial\n")
print(Table2)

# model estimation  

lambda_exp <- 1/mean(data)

theta_IE <- optimize(
  function(th)
    -sum(log(pmax(dIE(data,th),1e-10))),
  c(.001,10))$minimum

# Standard Ishita parameter estimation (m=2)
theta_std <- optimize(
  function(th)
    -sum(log(pmax(dIP(data, th, 2), 1e-10))),
  c(.001, 10))$minimum

# Transmuted Ishita parameter estimation (optim 2D)
fit_trans <- optim(par = c(theta = 0.15, lambda = 0.1), fn = function(p) {
  th <- p[1]; lam <- p[2]
  if(th <= 0 || lam < -1 || lam > 1) return(1e10)
  dens <- dTransIshita(data, th, lam)
  if(any(dens <= 0)) return(1e10)
  -sum(log(dens))
})
theta_trans <- fit_trans$par[1]
lambda_trans <- fit_trans$par[2]

# Power Ishita parameter estimation (optim 2D)
fit_power <- optim(par = c(theta = 0.15, alpha = 1.0), fn = function(p) {
  th <- p[1]; alp <- p[2]
  if(th <= 0 || alp <= 0) return(1e10)
  dens <- dPowerIshita(data, th, alp)
  if(any(dens <= 0)) return(1e10)
  -sum(log(dens))
})
theta_power <- fit_power$par[1]
alpha_power <- fit_power$par[2]


xx <- seq(min(data),max(data),length=400)

# MULTIPLE m ANALYSIS (2 to 5)

mgrid <- 2:5
fit_all <- data.frame()

for(m in mgrid){
  th <- estimate_theta(data,m)
  
  fit_all <- rbind(fit_all,
                   data.frame(m=m,theta=th))
}

print(fit_all)

# Prepare plot data
plot_data <- data.frame()

for(i in 1:nrow(fit_all)){
  
  m_val <- fit_all$m[i]
  th_val <- fit_all$theta[i]
  
  temp <- data.frame(
    x = xx,
    m = factor(m_val),
    pdf = dIP(xx, th_val, m_val),
    cdf = pIP(xx, th_val, m_val),
    hazard = hIP(xx, th_val, m_val)
  )
  
  plot_data <- rbind(plot_data, temp)
}

# PDF comparison
ggplot(plot_data,aes(x,pdf,color=m))+
  geom_line()+
  labs(title="")

# CDF comparison
ggplot(plot_data,aes(x,cdf,color=m))+
  geom_line()+
  labs(title="")

# Hazard comparison
ggplot(plot_data,aes(x,hazard,color=m))+
  geom_line()+
  labs(title="")

# Histogram + all PDFs
ggplot(data.frame(data),aes(data))+
  geom_histogram(aes(y=..density..),
                 bins=30,fill="grey70")+
  geom_line(data=plot_data,
            aes(x=x,y=pdf,color=m))

# Order Selection
ggplot(fit_table,
       aes(m,logSpacing))+
  geom_line(color="brown")+
  geom_point(size=3)+
  labs(title="")

# Model comparison
ggplot(data.frame(data),
       aes(x=data))+
  
  geom_histogram(aes(y=..density..,
                     fill="Observed Data"),
                 bins=30,
                 alpha=.6,
                 color="black")+
  
  stat_function(
    aes(color="Ishita Polynomial"),
    fun=function(z)
      dIP(z,theta_hat,m_hat),
    linewidth=1.4)+
  
  stat_function(
    aes(color="Standard Ishita"),
    fun=function(z)
      dIP(z,theta_std,2),
    linewidth=1.2)+
  
  stat_function(
    aes(color="Transmuted Ishita"),
    fun=function(z)
      dTransIshita(z,theta_trans,lambda_trans),
    linewidth=1.2)+
  
  stat_function(
    aes(color="Power Ishita"),
    fun=function(z)
      dPowerIshita(z,theta_power,alpha_power),
    linewidth=1.2)+
  
  stat_function(
    aes(color="Inverse Exponential"),
    fun=function(z)
      dIE(z,theta_IE),
    linewidth=1.2)+
  
  scale_color_manual(
    name="Fitted Models",
    values=c("Ishita Polynomial"="red",
             "Standard Ishita"="blue",
             "Transmuted Ishita"="darkgreen",
             "Power Ishita"="purple",
             "Inverse Exponential"="orange"))+
  
  scale_fill_manual(values=c(
    "Observed Data"="grey80"))+
  
  labs(title="",
       x="Failure Time",
       y="Density")+
  theme_minimal()

# box plot
ggplot(data.frame(Model="Ishita Polynomial",
                  data=data),
       aes(x=Model,y=data))+
  
  geom_boxplot(fill="skyblue",
               color="black",
               width=.4,
               alpha=.7)+
  
  stat_summary(fun=mean,
               geom="point",
               shape=18,
               size=4,
               color="red")+
  
  labs(title=
         "Boxplot under Ishita Polynomial Model",
       x="",y="Failure Time")+
  theme_minimal(base_size=14)

# Goodness of fit comparison

logLik_IP <- sum(log(pmax(
  dIP(data,theta_hat,m_hat),1e-10)))

logLik_Std <- sum(log(pmax(
  dIP(data,theta_std,2),1e-10)))

logLik_Trans <- sum(log(pmax(
  dTransIshita(data,theta_trans,lambda_trans),1e-10)))

logLik_Power <- sum(log(pmax(
  dPowerIshita(data,theta_power,alpha_power),1e-10)))

logLik_IE <- sum(log(pmax(
  dIE(data,theta_IE),1e-10)))

AIC <- function(k,ll)-2*ll+2*k
BIC <- function(k,ll,n)-2*ll+k*log(n)

KS_IP  <- ks.test(data,
                  function(q)
                    pIP(q,theta_hat,m_hat))

KS_Std <- ks.test(data,
                  function(q)
                    pIP(q,theta_std,2))

KS_Trans <- ks.test(data,
                    function(q)
                      pTransIshita(q,theta_trans,lambda_trans))

KS_Power <- ks.test(data,
                    function(q)
                      pPowerIshita(q,theta_power,alpha_power))

KS_IE  <- ks.test(data,
                  function(q)
                    pIE(q,theta_IE))

Table3 <- data.frame(
  Model=c("Ishita Polynomial",
          "Standard Ishita",
          "Transmuted Ishita",
          "Power Ishita",
          "Inverse Exponential"),
  
  AIC=c(AIC(2,logLik_IP),
        AIC(1,logLik_Std),
        AIC(2,logLik_Trans),
        AIC(2,logLik_Power),
        AIC(1,logLik_IE)),
  
  BIC=c(BIC(2,logLik_IP,n),
        BIC(1,logLik_Std,n),
        BIC(2,logLik_Trans,n),
        BIC(2,logLik_Power,n),
        BIC(1,logLik_IE,n)),
  
  KS_Statistic=c(
    KS_IP$statistic,
    KS_Std$statistic,
    KS_Trans$statistic,
    KS_Power$statistic,
    KS_IE$statistic),
  
  KS_pvalue=c(
    KS_IP$p.value,
    KS_Std$p.value,
    KS_Trans$p.value,
    KS_Power$p.value,
    KS_IE$p.value)
)

cat("\nTable 3: Goodness-of-Fit Comparison\n")
print(Table3)

