#### Explorations of the SECI-Movement Model

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "deSolve")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/pdr1_epidemic_model_functions.R")


# Initial state variables
S0 <- 99; I0 <- 1; E0 <- 0; C0 <- 0
U0 <- 200; V0 <- 0

# Parameter vector
# Order: alpha_c = x[1], alpha_i = x[2], 
#   beta_c = x[3], beta_i = x[4],
#   phi_muc = x[5], phi_mui = x[6], phi_pc = x[7], phi_pi = x[8],
#   delta = x[9], gamma = x[10], 
#   nu = x[11], lambda = x[12], 
#   T = x[13])


alpha_c = 0.3; alpha_i = 0.3
beta_c = 0.6; beta_i = 0.6
phi_pc = 0.5; phi_pi = 0.5; phi_mui = 0; phi_muc = seq(0,1,0.2) 
delta = 0.66; gamma = 0.000833 
nu = 0.001; lambda = 0.02
T = 0.5

params <- cbind(alpha_c, alpha_i, beta_c, beta_i, 
                phi_pc, phi_pi, phi_muc, phi_mui,
                delta, gamma,
                nu, lambda, T)


# Exploring transient dynamics
test <- SECIMDynamics(params[1,])

matplot(test[,-1], type = "l",
        col = c("black", "pink", "green", "brown", "grey", "green", "brown"), 
        lty = c(1,1,1,1,2,2,2), lwd = 2)
legend("topright", c("S", "E", "C", "I", "U", "V_c", "V_i"), 
       col = c("black", "pink", "green", "brown", "grey", "green", "brown"), 
       lty = c(1,1,1,1,2,2,2), lwd = 2)


# Exploring sensitivity to various parameters
testSim <- params %>% apply(., 1, SECIMSimulations) %>% rbindlist() %>% as.data.frame()
testSim



#### Exploring the transmission probability term
pinoc <- function(x){
  x <- as.numeric(x)
  phi <- x[1]; beta <- x[2]; T <- x[3]
  b <- 1 - exp(-beta*T)
  inoc <- phi*b
  return(c(b, inoc))
}

inocPars <- data.frame(phi = 0.5,
                       beta = seq(0,10,0.5),
                       T = 0.5)
t(apply(inocPars, 1, pinoc))
