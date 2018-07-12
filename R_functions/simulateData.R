#### Function to generate simulated data and remove negative (nonsense) values
# Assume normal distribution
simulateData <- function(mean, SE, nsim = 10000, nmin = 1000){
  # nsim = number of simulations
  # nmin = minimum number of final simulations,
  # nmin used to eliminate negative values and maintain symmetry of normal distributions
  rawSim <- rnorm(nsim, mean = mean, sd = SE) # vector of parameter estimates
  simPos <- rawSim[rawSim >= 0] # remove estimates < 0
  simDiff <- length(rawSim) - length(simPos)
  # Remove an equal number of values in the right tail of the distribution
  simVec <- simPos[-order(simPos)[(length(simPos)-simDiff):length(simPos)]][1:nmin]
  return(simVec)
}
