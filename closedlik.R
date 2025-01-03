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
      lik <- lik + x[i,t]*log(p[t]) + (1-x[i,t])*log(1-p[t])
    }
  }
  
  # Calculate the (log) probability of not being observed within the study  
  
  noprob <- sum(log(1-p[]))
  
  # Add the log-likelihood contribution of the probability of unobserved individuals
  
  lik <- lik + (N-n)*noprob
  
  # Add the Multinomial coefficient likelihood component:
  
  lik <- lik + lgamma(N+1) - lgamma(N-n+1)
  
  # Output the log-likelihood value:
  
  lik
}
