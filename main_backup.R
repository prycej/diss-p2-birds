library(data.table)
library(reshape2)
library(stats4)
library(rjags)

blackcap_data_raw <- fread("/Users/jacobpryce/Library/CloudStorage/OneDrive-UniversityofEdinburgh/diss-p2-birds/CR.blackcap_FixRing.csv", header = TRUE)

blackcap_data <- blackcap_data_raw[,5:11]

data_long <- melt(data)

closed_data <- fread("/Users/jacobpryce/Library/CloudStorage/OneDrive-UniversityofEdinburgh/diss-p2-birds/closed_sample.data", header = TRUE, drop = c(1))

closedlik <- function(theta, x, n, T) {

  # The data are stored in the array x;
  # n = number of observed individuals; T = number of capture occasions
  # Theta stores the set of parameter values - specified on the real line.
  # Define the parameter values in terms of capture probabs and population size.
  # Use the transformations: logit p[1:T] = theta[1:T]; log N = theta[T+1]

 p <- exp(theta[1:T])/(1+exp(theta[1:T]))
  unobs <- exp(theta[T+1])

  N <- unobs + n
  # Initialise the log-likelihood value:
  lik <- 0

  # Calculate the (log-)likelihood component for observed individual capture histories
  for (i in 1:n){
    for (t in 1:T){
      lik <- lik + x[i,t, with = FALSE]*log(p) + (1-x[i,t, with = FALSE])*log(1-p)
    }
  }

  # Calculate the (log) probability of not being observed within the study
  noprob <- sum(log(1-p))

  # Add the log-likelihood contribution of the probability of unobserved individuals
  lik <- lik + (N-n)*noprob

  # Add the Multinomial coefficient likelihood component:
  lik <- lik + lgamma(N+1) - lgamma(N-n+1)

  # Output the log-likelihood value:
  lik
}

result_test <- optim(par = c(rep(log(0.2), ncol(data)),log(100)), fn = closedlik, x = data, n = nrow(data), T = ncol(data), control=list(trace=1))

result <- optim(par = c(rep(log(0.2), ncol(blackcap_data)),log(100)), fn = closedlik, x = blackcap_data, n = nrow(blackcap_data), T = ncol(data), control=list(trace=1))



closedlik <- function(theta, x, n, T) {

  # The data are stored in the array x;
  # n = number of observed individuals; T = number of capture occasions
  # Theta stores the set of parameter values - specified on the real line.
  # Define the parameter values in terms of capture probabs and population size.
  # Use the transformations: logit p[1:T] = theta[1:T]; log N = theta[T+1]

  p <- exp(theta[1:T])/(1+exp(theta[1:T]))
  #  p <- exp(theta[1])
  unobs <- exp(theta[T+1])
  # unobs <- exp(theta[2])
  N <- unobs + n
  # Initialise the log-likelihood value:
  lik <- 0
  # Calculate the (log-)likelihood component for observed individual capture histories
  for (i in 1:n){
    for (t in 1:T){
      lik <- lik + x[i,t, with = FALSE]*log(p) + (1-x[i,t, with = FALSE])*log(1-p)
    }
  }
  # Calculate the (log) probability of not being observed within the study
  noprob <- sum(log(1-p))
  # Add the log-likelihood contribution of the probability of unobserved individuals
  lik <- lik + (N-n)*noprob
  # Add the Multinomial coefficient likelihood component:
  lik <- lik + lgamma(N+1) - lgamma(N-n+1)
  # Output the log-likelihood value:
  lik
}




