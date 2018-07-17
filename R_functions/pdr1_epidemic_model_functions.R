#### FUNCTIONS FOR SECI EPIDEMIC MODEL


#### Function to specify model equations
# Take the SECI model and then add two parameters determining vector movement and feeding
# phi[p,k] = attraction rate for k = C,I
# phi[mu,k] = leaving rate for k = C,I
# In Madden et al., phi is the number of plants visited per unit time -- is my interpretation of phi valid?
# T = time spent feeding per plant visit
# Host turnover changed to nu (used to be a)
# Vector turnover changed to lambda (used to be mu)


SECIMovementModel <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    # Vector-SECI model with preference decomposed to attraction and leaving
    dS <- nu*(E + C + I) - (phi_muc*(1-exp(-beta_c*T))*V_c*S)/(I + S + E + C) - (phi_mui*(1-exp(-beta_i*T))*V_i*S)/(I + S + E + C)
    dE <- (phi_muc*(1-exp(-beta_c*T))*V_c*S)/(I + S + E + C) + (phi_mui*(1-exp(-beta_i*T))*V_i*S)/(I + S + E + C) - (delta + nu)*E
    dC <- delta*E - nu*C - gamma*C
    dI <- gamma*C - nu*I
    dU <- lambda*(V_c + V_i) - (phi_pc*(1 - exp(-alpha_c*T))*C*U)/(I + S + E + C) - (phi_pi*(1 - exp(-alpha_i*T))*I*U)/(I + S + E + C)
    dV_c <- (phi_pc*(1 - exp(-alpha_c*T))*C*U)/(I + S + E + C) - lambda*V_c
    dV_i <- (phi_pi*(1 - exp(-alpha_i*T))*I*U)/(I + S + E + C) - lambda*V_i
    return(list(c(dS, dE, dC, dI, dU, dV_c, dV_i)))
  })
}



#### Function to return all State variables over all time steps for a single combination of parameters
SECIMDynamics <- function(x){
  x <- as.numeric(x)
  # Submit the parameter values as a row vector in the following order:
  Pars <- c(alpha_c = x[1], alpha_i = x[2], 
            beta_c = x[3], beta_i = x[4],
            phi_pc = x[5], phi_pi = x[6], phi_muc = x[7], phi_mui = x[8], 
            delta = x[9], gamma = x[10], 
            nu = x[11], lambda = x[12], 
            T = x[13])
  State <- c(S = S0, E = E0, C = C0, I = I0, U = U0, V_c = V0, V_i = V0) # Define starting values for states
  Time <- seq(0, 5000, by = 1) # time steps
  model.out <- as.data.frame(ode(func = SECIMovementModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  return(model.out)
}


#### Function to run simulations and output infective host (I) and infectious vector (V) densities
#### For multiple combinations of parameters
SECIMSimulations <- function(x){
  x <- as.numeric(x)
  # Submit the parameter values as a row vector in the following order:
  Pars <- c(alpha_c = x[1], alpha_i = x[2], 
            beta_c = x[3], beta_i = x[4],
            phi_pc = x[5], phi_pi = x[6], phi_muc = x[7], phi_mui = x[8], 
            delta = x[9], gamma = x[10], 
            nu = x[11], lambda = x[12], 
            T = x[13])
  State <- c(S = S0, E = E0, C = C0, I = I0, U = U0, V_c = V0, V_i = V0) # Define starting values for states
  Time <- seq(0, 5000, by = 1) # time steps
  model.out <- as.data.frame(ode(func = SECIMovementModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  # Need to remove some time steps because under/over flow results in NAs
  model.nona <- model.out[!is.na(model.out$I),]
  model.dat <- model.nona[nrow(model.nona),c("time", "C", "I", "V_c", "V_i")]
  #model.dat$trt <- x[14]
  return(model.dat)
}


#####################################################################################################################
#### 2-patch PdR1 model
#####################################################################################################################

#### Define model equations
SECIMPatchModel <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    # Vector-SECI model with preference decomposed to attraction and leaving
    # Resistant patch -- PdR1 patch
    dSr <- nu*(Er + Cr + Ir) - (phi_mucr*(1-exp(-beta_cr*T))*Vr*Sr)/(Ir + Sr + Er + Cr) - (phi_muir*(1-exp(-beta_ir*T))*Vr*Sr)/(Ir + Sr + Er + Cr)
    dEr <- (phi_mucr*(1-exp(-beta_cr*T))*Vr*Sr)/(Ir + Sr + Er + Cr) + (phi_muir*(1-exp(-beta_ir*T))*Vr*Sr)/(Ir + Sr + Er + Cr) - deltar*Er - nu*Er
    dCr <- deltar*Er - nu*Cr - gammar*Cr
    dIr <- gammar*Cr - nu*Ir
    dUr <- lambda*Vr - (phi_pcr*(1 - exp(-alpha_cr*T))*Cr*Ur)/(Ir + Sr + Er + Cr) - (phi_pir*(1 - exp(-alpha_ir*T))*Ir*Ur)/(Ir + Sr + Er + Cr) - mr*Ur + ms*Us + (Ss0/(Sr0+Ss0))*Mu
    dVr <- (phi_pcr*(1 - exp(-alpha_cr*T))*Cr*Ur)/(Ir + Sr + Er + Cr) + (phi_pir*(1 - exp(-alpha_ir*T))*Ir*Ur)/(Ir + Sr + Er + Cr) - lambda*Vr - mr*Vr + ms*Vs + (Sr0/(Sr0+Ss0))*Mv
    # Susceptible patch
    dSs <- nu*(Es + Cs + Is) - (phi_mucs*(1-exp(-beta_cs*T))*Vs*Ss)/(Is + Ss + Es + Cs) - (phi_muis*(1-exp(-beta_is*T))*Vs*Ss)/(Is + Ss + Es + Cs)
    dEs <- (phi_mucs*(1-exp(-beta_cs*T))*Vs*Ss)/(Is + Ss + Es + Cs) + (phi_muis*(1-exp(-beta_is*T))*Vs*Ss)/(Is + Ss + Es + Cs) - deltas*Es - nu*Es
    dCs <- deltas*Es - nu*Cs - gammas*Cs
    dIs <- gammas*Cs - nu*Is
    dUs <- lambda*Vs - (phi_pcs*(1 - exp(-alpha_cs*T))*Cs*Us)/(Is + Ss + Es + Cs) - (phi_pis*(1 - exp(-alpha_is*T))*Is*Us)/(Is + Ss + Es + Cs) - ms*Us + mr*Ur + (Ss0/(Sr0+Ss0))*Mu
    dVs <- (phi_pcs*(1 - exp(-alpha_cs*T))*Cs*Us)/(Is + Ss + Es + Cs) + (phi_pis*(1 - exp(-alpha_is*T))*Is*Us)/(Is + Ss + Es + Cs) - lambda*Vs - ms*Vs + mr*Vr + (Ss0/(Sr0+Ss0))*Mv
    return(list(c(dSr, dEr, dCr, dIr, dUr, dVr, dSs, dEs, dCs, dIs, dUs, dVs)))
  })
}


#### Function to return all State variables over all time steps for a single combination of parameters
SECIMPatchDynamics <- function(x){
  x <- as.numeric(x)
  # Submit the parameter values as a row vector in the following order,
  # Resistant patch parameters first, then susceptible patch, than patch-invariant parameters last
  Pars <- c(alpha_cr = x[1], alpha_ir = x[2], beta_cr = x[3], beta_ir = x[4],
            phi_pcr = x[5], phi_pir = x[6], phi_mucr = x[7], phi_muir = x[8],
            deltar = x[9], gammar = x[10], 
            alpha_cs = x[11], alpha_is = x[12], beta_cs = x[13], beta_is = x[14],
            phi_pcs = x[15], phi_pis = x[16], phi_mucs = x[17], phi_muis = x[18],
            deltas = x[19], gammas = x[20],
            nu = x[21], lambda = x[22], T = x[23],
            mr = x[24], ms = x[25], Mu = x[26], Mv = x[27])
  # Define starting values for states
  State <- c(Sr = x[28], Er = 0, Cr = 0, Ir = 0, Ur = Ur0, Vr = Vr0,
             Ss = Ss0, Es = 0, Cs = 0, Is = 0, Us = Us0, Vs = Vs0) 
  Time <- seq(0, 5000, by = 1) # time steps
  model.out <- as.data.frame(ode(func = SECIMPatchModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  return(model.out)
}


#### Function to run simulations and output infective host (I) and infectious vector (V) densities
#### For multiple combinations of parameters
#### Looking at 
SECIMPatchSimulations <- function(x){
  x <- as.numeric(x)
  # Submit the parameter values as a row vector in the following order,
  # Resistant patch parameters first, then susceptible patch, than patch-invariant parameters last
  Pars <- c(alpha_cr = x[1], alpha_ir = x[2], beta_cr = x[3], beta_ir = x[4],
            phi_pcr = x[5], phi_pir = x[6], phi_mucr = x[7], phi_muir = x[8],
            deltar = x[9], gammar = x[10], 
            alpha_cs = x[11], alpha_is = x[12], beta_cs = x[13], beta_is = x[14],
            phi_pcs = x[15], phi_pis = x[16], phi_mucs = x[17], phi_muis = x[18],
            deltas = x[19], gammas = x[20],
            nu = x[21], lambda = x[22], T = x[23],
            mr = x[24], ms = x[25], Mu = x[26], Mv = x[27],
            Sr0 = x[28])
  # Define starting values for states
  State <- c(Sr = x[28], Er = 0, Cr = 0, Ir = 0, Ur = Ur0, Vr = Vr0,
             Ss = Ss0, Es = 0, Cs = 0, Is = 0, Us = Us0, Vs = Vs0) 
  Time <- seq(0, 500, by = 1) # time steps
  model.out <- as.data.frame(ode(func = SECIMPatchModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  # Need to remove some time steps because under/over flow results in NAs
  model.nona <- model.out[!is.na(model.out$Ir),]
  model.dat <- model.nona[nrow(model.nona),c("time", "Sr", "Er", "Cr", "Ir", "Vr", "Ss", "Es", "Cs", "Is", "Vs")]
  return(model.dat)
}
