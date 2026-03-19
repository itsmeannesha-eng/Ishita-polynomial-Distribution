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

mgrid <- 2:6
fit_table <- data.frame()

for(m in mgrid){
  
  th <- estimate_theta(data,m)
  
  fit_table <- rbind(fit_table,
                     data.frame(m=m,
                                theta=th,
                                logSpacing=
                                  logSpacing(th,data,m)))
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

xx <- seq(min(data),max(data),length=400)

#PDF

ggplot()+
  geom_line(aes(xx,
                dIP(xx,theta_hat,m_hat)),
            color="blue",linewidth=1.4)+
  labs(title="Estimated PDF of Ishita Polynomial",
       x="Failure Time",y="Density")

# CDF

ggplot()+
  geom_line(aes(xx,
                pIP(xx,theta_hat,m_hat)),
            color="red",linewidth=1.4)+
  labs(title="Estimated CDF",
       x="Failure Time",y="CDF")

# Hazard

ggplot()+
  geom_line(aes(xx,
                hIP(xx,theta_hat,m_hat)),
            color="darkgreen",linewidth=1.4)+
  labs(title="Hazard Rate Function",
       x="Failure Time",y="Hazard Rate")

# Histogram
ggplot(data.frame(data),
       aes(data))+
  geom_histogram(aes(y=..density..),
                 bins=30,
                 fill="skyblue",
                 color="black")+
  stat_function(
    fun=function(z)
      dIP(z,theta_hat,m_hat),
    color="red",
    linewidth=1.5)+
  labs(title="Histogram with Ishita Fit",
       x="Failure Time",y="Density")

# QQ plot

ggplot(data.frame(sample=data),
       aes(sample=sample))+
  stat_qq(color="purple")+
  stat_qq_line()+
  labs(title="Q–Q Plot")

# TTT

xs <- sort(data)
TTT <- cumsum(xs)/sum(xs)
r <- (1:n)/n

ggplot(data.frame(r,TTT),
       aes(r,TTT))+
  geom_line(color="orange",
            linewidth=1.4)+
  geom_abline(linetype=2)+
  labs(title="Total Time on Test Plot")

# Order Selection

ggplot(fit_table,
       aes(m,logSpacing))+
  geom_line(color="brown")+
  geom_point(size=3)+
  labs(title="Polynomial Order Selection")

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
    aes(color="Exponential Polynomial"),
    fun=function(z)
      dexp(z,lambda_exp),
    linewidth=1.3)+
  
  stat_function(
    aes(color="Inverse Exponential Polynomial"),
    fun=function(z)
      dIE(z,theta_IE),
    linewidth=1.3)+
  
  scale_color_manual(
    name="Fitted Models",
    values=c("Ishita Polynomial"="red",
             "Exponential Polynomial"="blue",
             "Inverse Exponential Polynomial"="darkgreen"))+
  
  scale_fill_manual(values=c(
    "Observed Data"="grey80"))+
  
  labs(title="Model Comparison",
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

logLik_EXP <- sum(log(pmax(
  dexp(data,lambda_exp),1e-10)))

logLik_IE <- sum(log(pmax(
  dIE(data,theta_IE),1e-10)))

AIC <- function(k,ll)-2*ll+2*k
BIC <- function(k,ll,n)-2*ll+k*log(n)

KS_IP  <- ks.test(data,
                  function(q)
                    pIP(q,theta_hat,m_hat))

KS_EXP <- ks.test(data,"pexp",
                  lambda_exp)

KS_IE  <- ks.test(data,
                  function(q)
                    pIE(q,theta_IE))

Table3 <- data.frame(
  Model=c("Ishita Polynomial",
          "Exponential",
          "Inverse Exponential"),
  
  AIC=c(AIC(2,logLik_IP),
        AIC(1,logLik_EXP),
        AIC(1,logLik_IE)),
  
  BIC=c(BIC(2,logLik_IP,n),
        BIC(1,logLik_EXP,n),
        BIC(1,logLik_IE,n)),
  
  KS_Statistic=c(
    KS_IP$statistic,
    KS_EXP$statistic,
    KS_IE$statistic),
  
  KS_pvalue=c(
    KS_IP$p.value,
    KS_EXP$p.value,
    KS_IE$p.value)
)

cat("\nTable 3: Goodness-of-Fit Comparison\n")
print(Table3)

