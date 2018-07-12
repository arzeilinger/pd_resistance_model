#####################################################################################################
#### 2-Patch SECI Movement Model
#### Expanded on model described in "pdr1_SECI_model.R" script
#### Uses the same (or similar) parameter values as the above script
#####################################################################################################

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "deSolve", "ggplot2")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/pdr1_epidemic_model_functions.R")
source("R_functions/simulateData.R")
source("R_functions/factor2numeric.R")

#########################################################################################
#### Parameter values

# Number of total Monte Carlo simulations to start out with and minimum to reduce to 
nsim <- 50000
nmin <- 5000

#### Import saved parameter values
## These parameter values were calculated from 2016 PdR1 transmission experiments
## They were calculated in the R script R_analyses/pdr1_SECI_model.R
parList <- readRDS("output/pdr1_SECI_parameter_values.rds")

#### Structuring the parameter estimates, same estimates as single patch model
resistantParams <- parList$R[,1:10]
names(resistantParams) <- paste(names(resistantParams), "r", sep = "")
susceptibleParams <- parList$S[,1:10]
names(susceptibleParams) <- paste(names(susceptibleParams), "s", sep = "")
invariantParams <- parList$R[,11:13]
mrVec <- rep(0.001, nmin) # Emmigration rate from Resistant patch
msVec <- rep(0.001, nmin) # Emmigration rate from Susceptible patch
M <- rep(1000, nmin) # Total number of migrant BGSS
# Infectivity of migrants
# Purcell 1975 estimated that 30% of BGSS in vineyards in early season were infectious
propV <- 0.3
MuVec <- M*(1-propV)
MvVec <- M*propV 

patchParams <- cbind(resistantParams, susceptibleParams, invariantParams, mrVec, msVec, MuVec, MvVec)


#### Initial state variables
Ss0 <- 100;
Us0 <- 0; Vs0 <- 0
Ur0 <- 0; Vr0 <- 0


##############################################################################################################################
#### Run 2-patch model over range of Resistant patch sizes
patchParams <- patchParams[1:1000,]

Sr0Vec <- seq(1,200,length.out = 10) %>% round()

patchParList <- vector("list", length(Sr0Vec))

for(i in 1:length(Sr0Vec)){
  patchParams$Sr0Vec <- Sr0Vec[i]
  patchParList[[i]] <- patchParams
}

# Understanding what the primary spread term looks like
# sAreaProp <- Ss0/(Sr0Vec + Ss0)
# sM <- sAreaProp*MvVec
# # Number of infectious migrant vectors in susceptible patch as a function of resistant patch area
# plot(Sr0Vec, sM[1:length(Sr0Vec)])

# Running a bunch of simulations
ti <- proc.time()
patchSimList <- lapply(patchParList, function(x) apply(x, 1, SECIMPatchSimulations) %>% rbindlist() %>% as.data.frame())
tf <- proc.time()
# Time the simulations took in minutes:
(tf-ti)/60

#saveRDS(patchSimList, file = "output/simulation_results_2-patch_area.rds")

for(i in 1:length(patchSimList)){
  patchSimList[[i]]$Sr0 <- Sr0Vec[i]
  patchSimList[[i]]$TI <- with(patchSimList[[i]], Cr + Ir + Cs + Is) # Total plants infected over all patches
}
# Convert to data.frame, remove Vc and Vi, and transform to "long" format
patchSimData <- patchSimList %>% rbindlist() %>% as.data.frame() %>% dplyr::select(., Cr, Ir, Cs, Is, TI, Vs, Sr0) %>% 
  gather(., key = state, value = density, Cr, Ir, Cs, Is, TI, Vs)
patchSimData$patchAreaRatio <- patchSimData$Sr0/Ss0


# Summarize simulation results
summaryPatch <- patchSimData %>% group_by(patchAreaRatio, state) %>% summarise(mean = mean(density),
                                                                               median = median(density),
                                                                               sd = sd(density),
                                                                               cil = quantile(density, 0.025),
                                                                               ciu = quantile(density, 0.975),
                                                                               max = max(density))

# Remove Vs
summaryPatch[summaryPatch$state == "Vs",]
summaryPatch <- summaryPatch %>% dplyr::filter(., state != "Vs", state != "Cr", state != "Ir", state != "TI")

# For plotting
summaryPatch$state <- ifelse(summaryPatch$state == "Cs", "HC", "HI")

#### Plotting with ggplot2
#### Mean infected density of C, I, and V
patchAreaPlot <- ggplot(data=summaryPatch, aes(x=patchAreaRatio, y=mean, group=state, shape=state)) +
  geom_line(aes(linetype=state), size=1.25) +
  #geom_point(position = position_dodge(width = 0.05), aes(shape=state), size=3.5) +
  #geom_errorbar(aes(ymax=ciu, ymin=cil, width=0.2), position = position_dodge(width = 0.05)) +
  scale_x_continuous(name = "Ratio Resistant patch : Susceptible patch area") +
  #                    breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Percent infected hosts in susceptible patch") +
  # ylab("% insects on source plant") + 
  #ylim(c(0,100)) +
  # xlab("Genotype") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

patchAreaPlot

ggsave("results/figures/SECI_patch_area_plot.jpg", plot = patchAreaPlot,
       width = 7, height = 7, units = "in")


###############################################################################################################
#### Simulation with pdr1 2016 experiment derived parameter values

# Add constant Sr0Vec to parameter matrix
patchParams$Sr0Vec <- rep(100, nmin)


# Testing out the model
patchDynamicsOut <- SECIMPatchDynamics(patchParams[5,])
# Plot only hosts because of larger vector numbers
plotDynamics <- patchDynamicsOut %>% dplyr::select(., -time, -starts_with("V"), -starts_with("U"))
matplot(plotDynamics[1:50,], type = "l")

# Testing the simulation function
patchSim <- SECIMPatchSimulations(patchParams[1,])

# Running a bunch of simulations
ti <- proc.time()
patchSim <- apply(patchParams, 1, SECIMPatchSimulations) %>% rbindlist() %>% as.data.frame()
tf <- proc.time()
# Time the simulations took in minutes:
(tf-ti)/60

saveRDS(patchSim, file = "output/simulation_results_2_patch_model.rds")
