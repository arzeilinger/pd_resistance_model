#####################################################################################################
#### 2-Patch SECI Movement Model
#### Expanded on model described in "pdr1_SECI_model.R" script
#### Uses the same (or similar) parameter values as the above script
#####################################################################################################

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "deSolve", "ggplot2", "akima")
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
runParams <- patchParams[1:600,]
runParams <- colMeans(patchParams) %>% t() %>% as.data.frame() # Just use the mean value of each parameter


Sr0Vec <- seq(1,200,length.out = 20) %>% round()

patchParList <- vector("list", length(Sr0Vec))

for(i in 1:length(Sr0Vec)){
  runParams$Sr0Vec <- Sr0Vec[i]
  patchParList[[i]] <- runParams
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


## Add initial resistant plant population (Sr0) to each list element, calculate total infected plants, and total healthy plants
for(i in 1:length(patchSimList)){
  patchSimList[[i]]$Sr0 <- Sr0Vec[i]
  patchSimList[[i]]$TI <- with(patchSimList[[i]], Cr + Ir + Cs + Is) # Total plants infected over all patches
  patchSimList[[i]]$healthys <- with(patchSimList[[i]], Ss + Es + Cs) # Total healthy/asymptomatic plants in susceptible patch
  patchSimList[[i]]$healthyr <- with(patchSimList[[i]], Sr + Er + Cr) # Total healthy/asymptomatic plants in resistant patch
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
  scale_x_continuous(name = "Ratio Resistant patch : Susceptible patch area") +
  scale_y_continuous(name = "Percent infected hosts in susceptible patch") +
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
#### Bioeconomic model
#### First iteration of the model -- Yield = (c1*healthys + c2*healthyr)[harvest]

## Define relative economic value of total healthy susceptible (healthys = Ss + Es + Cs) and health resistant (healthyr = Sr + Er + Cr) grapevines at time of harvest
## Assumes that Yield is proporitional to healthy grapevines at the end of the numerical simulations (harvest time)
## If c1 > c2, then resistant grapevines have lower value
c1 <- 1
c2 <- 0.1


## Calculate Yield for each simulation, from patchSimList output 
for(i in 1:length(patchSimList)){
  patchSimList[[i]]$Yield <- with(patchSimList[[i]], c1*healthys + c2*healthyr)
}

## Convert to data.frame
YieldData <- patchSimList %>% rbindlist() %>% as.data.frame()
YieldData$patchAreaRatio <- YieldData$Sr0/Ss0

hist(YieldData[YieldData$patchAreaRatio == 2,]$Yield, breaks = c(0,4,by=0.1))

# Summarize simulation results
summaryYield <- YieldData %>% group_by(patchAreaRatio) %>% summarise(mean = mean(Yield),
                                                                     median = median(Yield),
                                                                     sd = sd(Yield),
                                                                     cil = quantile(Yield, 0.025),
                                                                     ciu = quantile(Yield, 0.975),
                                                                     max = max(Yield))

#### Plotting with ggplot2
#### Yield and 95% confidence intervals from total healthy hosts
yieldAreaPlot <- ggplot(data=summaryYield) +
  geom_line(aes(x=patchAreaRatio, y=median), size=1.25, linetype = 1) +
  # Include lines for 95% confidence interval
  geom_line(aes(x=patchAreaRatio, y=cil), size=1.25, linetype = 2) +
  geom_line(aes(x=patchAreaRatio, y=ciu), size=1.25, linetype = 2) +
  scale_x_continuous(name = "Ratio Resistant patch : Susceptible patch area") +
  scale_y_continuous(name = "Yield") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

yieldAreaPlot

ggsave("results/figures/yield_patch_area_plot.jpg", plot = yieldAreaPlot,
       width = 7, height = 7, units = "in")




###############################################################################################################
#### Run Bioeconomic model over a range of c2 values, and Sr and Ss outputs from the epidemic simulation

## Function to calculate Yield when varying c2, Sr, and Ss
## Sr and Ss outputs are in list object patchSimList
# x[1] = healthy hosts in resistant patch
# x[2] = value of yield from resistant vines relative to susceptible vines, assuming c1 = 1
# x[3] = healthy hosts in susceptible patch
simpleYield <- function(x){
  Yield <- c1*x[3] + x[2]*x[1]
}

## Define a range of c2 values
totalHealthy <- patchSimList %>% rbindlist %>% as.data.frame %>% dplyr::select(healthys, healthyr, Sr0)
c2values <- seq(0.01, 10, length.out = 20)
## Create combinations of each set of values of totalHealthy and c2 and combine into a data.frame
yieldInputs <- expand.grid(healthyr = totalHealthy$healthyr, c2 = c2values) %>% 
  left_join(., totalHealthy, by = "healthyr")

# Set c1
c1 <- 1

## Calculate yield for each combination of totalHealthy hosts and c2, then clean it up for plotting
YieldData <- apply(yieldInputs, 1, simpleYield) %>% cbind(., yieldInputs$c2, yieldInputs$Sr0) %>% as.data.frame()
names(YieldData) <- c("yield", "c2", "Sr0")
YieldData$patchAreaRatio <- YieldData$Sr0/Ss0


#### Make a contourplot:
# resistant patch area (relative) on x-axis
# c2 value on y-axis
# yield as contours/colors/z-axis
zzyield <- interp(x = YieldData$patchAreaRatio, y = log10(YieldData$c2), z = YieldData$yield)
tiff("results/figures/contourplot_yield.tiff")
  filled.contour(zzyield, col = topo.colors(24),  
                 xlab = "Ratio Resistant patch : Susceptible patch area",
                 ylab = "Relative value of Resistant cultivar (log10(c2))",
                 cex.lab = 1.3, cex.axis = 1.2, pty = "s")
dev.off()


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
