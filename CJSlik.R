CJSlik <- function(theta, x, f, l, n, T) {
  
  # The data are stored in the array x; 
  # n = number of observed individuals; T = number of capture occasions
  
  # f - array corresponding to first time each individual is observed
  # f can be calculated using e.g.: 
  #   for (i in 1:n){f[i] <- which(x[i,]==1)[1]}
  
  # l - array corresponding to last time each individual is observed
  # l can be calculated using e.g.: 
  #   for (i in 1:n){l[i] <- which(x[i,]==1)[length(which(x[i,]==1))]}
  
  # Theta stores the set of parameter values - specified on the real line.
  # Define the parameter values in terms of capture probabs and population size.
  # Use the transformations: logit phi[1:T-1] = theta[1:T-1]; logit p = theta[T]
  
  phi <- exp(theta[1:(T-1)])/(1+exp(theta[1:(T-1)]))
  p <- exp(theta[T])/(1+exp(theta[T]))
  
  # Calculate the chi terms: probability of not observed after time t
  # Initially set up chi to be an array of length T all elements equal to 1
  # chi[T] <- 1
  # Calculate chi[1:T-1] using recursive formula
  
  chi <- array(1,T)
  
  for (t in (T-1):1){
    chi[t] <- 1 - phi[t] + phi[t]*(1-p)*chi[t+1]
  }
  
  # Initialise the log-likelihood value for each individual:
  
  lik <- array(0,n)
  
  # Calculate the (log-)likelihood component for observed individual capture histories
  
  for (i in 1:n){
    
    # Consider the times between initial capture and final capture for each individual
    # Contribution only required when l[i] is not equal to f[i]
    
    if (f[i] != l[i]){
    
      for (t in f[i]:(l[i]-1)){
        lik[i] <- lik[i] + log(phi[t]) + x[i,t+1]*log(p) + (1-x[i,t+1])*log(1-p)
      }
    }
    
    # Note the above function assumes that we do not have any observed
    # histories with a single observation at the final capture occasion
    # I.e. we no not have any histories such that f = l = T. 
    # For such histories there is no information (and so can be omitted).
    
    # Add chi terms (probability of not being observed after time l[i])
    
    lik[i] <- lik[i] + log(chi[l[i]])
    
  }
  
# Calculate log-likelihood over all individuals:
  
  sumlik <- sum(lik)
  
  # Output the log-likelihood value:
  
  sumlik
}
