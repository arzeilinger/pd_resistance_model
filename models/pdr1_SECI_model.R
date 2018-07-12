#### Pierce's Disease SECI epidemic Model
#### Developed by Adam Zeilinger, based on Zeilinger and Daugherty (2014)

rm(list = ls())
# Load packages
my.packages <- c("tidyr", "dplyr", "data.table", "deSolve", "ggplot2")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/pdr1_epidemic_model_functions.R")
source("R_functions/simulateData.R")
source("R_functions/factor2numeric.R")

##########################################################################################
#### Compiling parameter estimates
##########################################################################################
#### Loading data sets
# Rate parameters
rateParams <- read.csv("results/pdr1_cmm_rate_parameter_estimates.csv", header = TRUE)
rateParams

# Transmission parameters
transdata <- readRDS("output/pdr1_transmission_preference_dataset.rds")
str(transdata)

# Construct source plant infection status column
transdata$source.cfu.per.g <- as.numeric(transdata$source.cfu.per.g)
transdata$source.infection <- transdata$source.cfu.per.g
transdata$source.infection[transdata$source.infection > 0] <- 1

# Transmission summary table
transdata$test.plant.infection <- as.numeric(transdata$test.plant.infection)
transSummary <- transdata %>% group_by(., week, trt) %>% 
  summarise(propInfected = sum(test.plant.infection, na.rm = TRUE)/sum(!is.na(test.plant.infection)),
            sepropI = sqrt(propInfected*(1 - propInfected)/8),
            meanPD = mean(pd_index, na.rm = TRUE),
            sePD = sd(pd_index, na.rm = TRUE)/sqrt(sum(!is.na(pd_index))),
            meancfuSource = mean(log10(source.cfu.per.g+1)),
            secfuSource = sd(log10(source.cfu.per.g+1))/sqrt(sum(!is.na(source.cfu.per.g))),
            sourcepInfected = sum(source.infection, na.rm = TRUE)/sum(!is.na(source.infection)),
            meancfuVector = mean(log10(cagecfu+1), na.rm = TRUE),
            secfuVector = sd(log10(cagecfu+1), na.rm = TRUE)/sqrt(sum(!is.na(cagecfu), na.rm = TRUE)),
            meanPropVectorInfected = mean(propVectorInfected, na.rm = TRUE),
            sePropVectorInfected = sd(propVectorInfected, na.rm = TRUE)/sqrt(sum(!is.na(propVectorInfected), na.rm = TRUE)))
transSummary

# Empty list of length 2 to hold the simulated parameter values
parList <- vector("list", 2)

# Number of total Monte Carlo simulations to start out with and minimum to reduce to 
nsim <- 50000
nmin <- 5000

for(i in 1:length(levels(rateParams$trt))){
  trt.i <- levels(rateParams$trt)[i]
  
  #### Rate parameters
  # phi_pc = 8S p1
  # phi_pi = 12S p1
  # phi_muc = 8S mu1
  # phi_mui = 12S mu1
  # Multiply by 24 to get rate in day-1
  
  phi_pcBar <- rateParams %>% dplyr::filter(., parameter == "p1" & week == 3 & trt == trt.i) %>% dplyr::select(., estimate, se)*24
  phi_pcVec <- simulateData(phi_pcBar$estimate, phi_pcBar$se, nsim = nsim, nmin = nmin)
  
  phi_piBar <- rateParams %>% dplyr::filter(., parameter == "p1" & week == 12 & trt == trt.i) %>% dplyr::select(., estimate, se)*24
  phi_piVec <- simulateData(phi_piBar$estimate, phi_piBar$se, nsim = nsim, nmin = nmin)
  
  phi_mucBar <- rateParams %>% dplyr::filter(., parameter == "mu1" & week == 3 & trt == trt.i) %>% dplyr::select(., estimate, se)*24
  phi_mucVec <- simulateData(phi_mucBar$estimate, phi_mucBar$se, nsim = nsim, nmin = nmin)
  
  phi_muiBar <- rateParams %>% dplyr::filter(., parameter == "mu1" & week == 12 & trt == trt.i) %>% dplyr::select(., estimate, se)*24
  phi_muiVec <- simulateData(phi_muiBar$estimate, phi_muiBar$se, nsim = nsim, nmin = nmin)
  
  #### Transmission parameters
  # beta_c ~ P(trans) 3
  # beta_i ~ P(trans) 12
  # Divide P(trans) by 8 to give the inoculation per vector
  # Assume that inoculation is 60% of acquisition, because I don't have acquisition data yet
  # Inoculation
  beta_cBar <- transSummary %>% dplyr::filter(., week == 3 & trt == trt.i) %>% dplyr::select(., propInfected, sepropI)
  beta_cVec <- simulateData(mean = beta_cBar$propInfected, SE = beta_cBar$sepropI, nsim = nsim, nmin = nmin)/8
  
  beta_iBar <- transSummary %>% dplyr::filter(., week == 12 & trt == trt.i) %>% dplyr::select(., propInfected, sepropI)
  beta_iVec <- simulateData(mean = beta_iBar$propInfected, SE = beta_iBar$sepropI, nsim = nsim, nmin = nmin)/8
  
  # Acquisition
  alpha_cBar <- transSummary %>% dplyr::filter(., week == 3 & trt == trt.i) %>% dplyr::select(., meanPropVectorInfected, sePropVectorInfected)
  alpha_cVec <- simulateData(mean = alpha_cBar$meanPropVectorInfected, SE = alpha_cBar$sePropVectorInfected, nsim = nsim, nmin = nmin)/8
  
  alpha_iBar <- transSummary %>% dplyr::filter(., week == 12 & trt == trt.i) %>% dplyr::select(., meanPropVectorInfected, sePropVectorInfected)
  alpha_iVec <- simulateData(mean = alpha_iBar$meanPropVectorInfected, SE = alpha_iBar$sePropVectorInfected, nsim = nsim, nmin = nmin)/8
  
  #### Latent and incubation period
  # Latent period: take as uniform distribution between 4 and 21 days (or as a rate, 1/4 - 1/21),
  # 4 is based on Hill and Purcell 1997
  # 21 is based on significant transmission at 3 weeks
  deltaVec <- runif(n = nmin, min = 1/21, max = 1/4)
  
  # Incubation period: 
  # I know the minimum length of incubation (for S, min = 8 - 3 = 5 weeks; for R, min = 12 - 3 = 9 weeks),
  # but for both lines, the uncertainty is around the end, or max, of the incubation period. Need Monte Carlo for that.
  # For susceptible line, assume incubation period is between 5 weeks and 9 weeks after latent period,
  # because I have no evidence of change in preference or in transmission at 8 weeks but I do have changes at 12 weeks
  # For resistant line, assume incubation is between 9 weeks and 16 - 3 = 13 weeks
  # saw no decline in preference at 12 weeks, but many plants looked very bad at 16 weeks, so use this as cutoff
  if(trt.i == "R"){
    gammaVec <- runif(n = nmin, min = 1/(13*7), max = 1/(9*7))
  } else {
    gammaVec <- runif(n = nmin, min = 1/(9*7), max = 1/(5*7))
    }
  
  #### Vector and host recovery
  # Take vector recovery from Fabien's experiment with WT plants
  #lambdaVec <- simulateData(mean = 0.083, SE = 0.042, nsim = nsim, nmin = nmin)
  # Alternatively, use vector recovery rate from Oikos paper
  lambdaVec <- 10
  
  # Take host recovery as the change in proportion of source plants that were negative between 3 and 12 weeks
  # sourcepInfected for 12S - sourcepInfected for 3S
  # then divide the number of days to get recovery per day
  # prop.i <- as.numeric(transSummary[transSummary$week == 3 & transSummary$trt == trt.i, "sourcepInfected"])
  # prop.f <- as.numeric(transSummary[transSummary$week == 12 & transSummary$trt == trt.i, "sourcepInfected"])
  # nuVec <- ifelse(prop.i >= prop.f, prop.i - prop.f, 0)/((12-3)*7)
  nuVec <- 0.0001
  
  #### Feeding time 
  # Estimate of T is taken from Almeida and Backus 2004, Table 4, WDEI calculation for waveform C1 (ingestion)
  # Duration originaly given in minutes, convert to fraction of a day
  TVec <- simulateData(mean = 48.1/60, SE = 12.2/60, nsim = nsim, nmin = nmin)
  
  #### Combine all parameter values
  # Parameter vector order for simulation function
  # alpha_c = x[1], alpha_i = x[2], 
  # beta_c = x[3], beta_i = x[4],
  # phi_pc = x[5], phi_pi = x[6], phi_muc = x[7], phi_mui = x[8], 
  # delta = x[9], gamma = x[10], 
  # nu = x[11], lambda = x[12], 
  # T = x[13])
  par.i <- as.data.frame(cbind(alpha_cVec, alpha_iVec,
                               beta_cVec, beta_iVec,
                               phi_pcVec, phi_piVec, phi_mucVec, phi_muiVec,
                               deltaVec, gammaVec,
                               nuVec, lambdaVec,
                               TVec))
  #### Include column of genotype
  #par.i$trt <- trt.i
  parList[[i]] <- par.i
}
names(parList) <- levels(rateParams$trt)

# Look at distribution of parameters
lapply(1:length(parList), function(x) summary(parList[[x]]))

# Save parameter values
saveRDS(parList, "output/pdr1_SECI_parameter_values.rds")

###########################################################################################
#### Numerical Simulations
###########################################################################################

# Initial state variables
S0 <- 99; I0 <- 1; E0 <- 0; C0 <- 0
U0 <- 200; V0 <- 0

ti <- proc.time()
simList <- lapply(parList, function(x) apply(x, 1, SECIMSimulations) %>% rbindlist() %>% as.data.frame())
tf <- proc.time()
# Time the simulations took:
(tf-ti)/60
# Remove the "time" column
simList <- lapply(1:length(simList), function(x) simList[[x]] %>% dplyr::select(., -time))

# Save simulation output
saveRDS(simList, file = "output/pdr1_SECIM_simulation_output.rds")
# Load simulation output
simList <- readRDS("output/pdr1_SECIM_simulation_output.rds")

# Get V = V_c + V_i
for(i in 1:length(simList)){
  simList[[i]]$V <- with(simList[[i]], V_c + V_i)
  simList[[i]]$trt <- levels(rateParams$trt)[i]
}
# Convert to data.frame, remove V_c and V_i, and transform to "long" format
dataSim <- simList %>% rbindlist() %>% as.data.frame() %>% dplyr::select(., -V_c, -V_i) %>% gather(., key = state, value = density, C, I, V)
dataSim$genotype <- ifelse(dataSim$trt == "R", "Resistant", "Susceptible")
dataSim$propInfected <- with(dataSim, ifelse(state == "V", density/200, density/100))

# Summarize simulation results
summarySim <- dataSim %>% group_by(genotype, state) %>% summarise(mean = mean(propInfected),
                                                                  median = median(propInfected),
                                                                  sd = sd(propInfected),
                                                                  cil = quantile(propInfected, 0.025),
                                                                  ciu = quantile(propInfected, 0.975),
                                                                  max = max(propInfected))
summarySim

summarySim[summarySim$state != "V",] %>% group_by(genotype) %>% summarise(infsum = sum(mean))
# I + C total hosts infected are similar but Susceptible is higher:
# Resistant = 0.84
# Susceptible = 0.98

# Change state names to HC, HI, and V
summarySim$state[summarySim$state == "C"] <- "HC"
summarySim$state[summarySim$state == "I"] <- "HI"

#### Plotting with ggplot2
#### Mean infected density of C, I, and V
propInfectedPlot <- ggplot(data=summarySim, aes(x=genotype, y=median, group=state, shape=state)) +
  # geom_bar(position=position_dodge(), stat="identity", 
  #          aes(fill=trt)) +
  # geom_hline(aes(yintercept=50), linetype="dashed") +
  # geom_line(aes(linetype=trt), size=1.25) +
  geom_point(position = position_dodge(width = 0.5), aes(shape=state), size=3.5) +
  geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd, width=0.2), position = position_dodge(width = 0.5)) +
  # scale_x_continuous(name = "Weeks post inoculation", 
  #                    breaks = c(3,8,12)) + 
  scale_y_continuous(name = "Proportion infected") +
  # ylab("% insects on source plant") + 
  # ylim(c(0,100)) +
  xlab("Genotype") +
  theme_bw(base_size=18) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 

propInfectedPlot

ggsave("results/figures/SECI_prop_infected_plot.jpg", plot = propInfectedPlot,
       width = 7, height = 7, units = "in")






#####################################################################################################
#### Transient dynamics

dynS <- SECIMDynamics(parList[[2]][1,])
matplot(dynS[1:200,-1], type = "l",
        col = c("black", "pink", "green", "brown", "grey", "green", "brown"), 
        lty = c(1,1,1,1,3,3,3), lwd = 2)
legend("topright", c("S", "E", "C", "I", "U", "V_c", "V_i"), 
       col = c("black", "pink", "green", "brown", "grey", "green", "brown"), 
       lty = c(1,1,1,1,3,3,3), lwd = 2)



#####################################################################################################
#### Histograms of parameter simulations and results

parS <- parList[[2]]
for(i in 1:ncol(parS)){
  fileName <- paste("results/figures/histogram_seci_parameter_simulations_", names(parS)[i], ".jpg", sep="")
  jpeg(fileName)
  hist(parS[,i], main = "", xlab = "")
  dev.off()
}

jpeg("results/figures/histogram_seci_infected_vector_density.jpg")
  hist(simList[[2]]$V, main = "", xlab = "")
dev.off()


