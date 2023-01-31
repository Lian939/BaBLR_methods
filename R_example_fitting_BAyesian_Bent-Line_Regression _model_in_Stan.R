#============================================================================
# Project: BAyesian Bent-Line Regression model(BaBLR) in Stan
# File:    This file is to create a single simulated dataset and fit 
#          BAyesian Bent-Line Regression model(BaBLR) in Stan. 
# The assumptions similar to the preliminary results using our proposed BaBLR 
# model from the Wisconsin Registry for Alzheimer's Prevention (WRAP) data 
#============================================================================

#============================================================================
#load packages (if the packages are not installed, we should install the 
#packages from CRAN first)
library(rstan)   # package for running Stan from within R
library(dplyr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#============================================================================

library(bayesplot)
# library(bayestestR)

#============================================================================
# Function to create a single simulated dataset.
# Notes: The function returns an R list, containing within it two objects:
#        (i) a vector with the named elements of the vector providing the
#            true parameter values used in creating the simulated data, and
#        (ii) a list containing the simulated data, to be utilised by RStan.
#        Note that the dataset must be provided as an R list, since this is
#        the format required by RStan. If you are analysing your own dataset,
#        then it will need to be converted from its usual format (for
#        example a data frame) to a named R list.
#============================================================================
createdata <-
  function(seed                     = sample.int(.Machine$integer.max, 1),
           Npat                     = 100,
           agemean                  = 2.5,
           agesd                    = 7.2,
           timepointmin             = 3,
           timepointmax             = 7,
           agemin                   = -20,
           agemax                   = 25,
           intercept.truemean       = -0.00588,
           intercept.truesd         = 0.638,
           preslope.truemean        = -0.00516,
           preslope.truesd          = 0.012,
           diffslope.truemean       = -0.00851,
           diffslope.truesd         = 0.150,
           cp.truemean              = 10,
           cp.truesd                = 10,
           cp.lowerbound             = -20,
           error.truesd             = 0.30) {
    set.seed(seed)
    
    # Number of measurements per patient
    Ntimepoints      <- as.integer(runif(Npat, timepointmin, timepointmax))
    id.long          <- rep(1:Npat, times = Ntimepoints)
    
    NstartAgeRandom      <- as.integer(runif(Npat, agemin, (agemax-4))) 
    
    #Age
    age <- runif(n = length(id.long), min = agemin, max = agemax)    
    # Measurement times
    ageRan<- runif(n=length(id.long),
                   min=-0.2,
                   max = 0.2)
    age <- age+ageRan
    
    # Measurement errors
    error            <-
      rnorm(n = length(id.long),
            mean = 0,
            sd = error.truesd)
    # Long data (measurement times and error terms)
    dat1             <-
      data.frame(id = id.long, age = age, error = error)
    
    # True random effects (uncorrelated): intercept, preslope, diffslope, change point
    id.short         <- c(1:Npat)
    X1<- rnorm(n = Npat,
               mean = intercept.truemean,
               sd = intercept.truesd)
    X2<- rnorm(n = Npat,
               mean = preslope.truemean,
               sd = preslope.truesd)
    X3<- rnorm(n = Npat,
               mean = diffslope.truemean,
               sd = diffslope.truesd)
    X4<- rnorm(n = Npat,
               mean = cp.truemean,
               sd = cp.truesd)
    # Short data (true random effect parameters: intercept, preslope, diffslope, change point)
    dat2             <- data.frame(id = id.short, X1,X2,X3,X4)
    # Merge long and short data
    simdat           <- merge(dat1, dat2, by = "id")
    simdat           <- simdat[order(simdat$id, simdat$age), ]
    
    # Observed outcome
    simdat <- simdat %>%
      mutate(y= case_when(age< X4 ~ X1 + X2 *(age-X4) +error,
                          age > X4 ~X1 +  (X2 + X3) *(age-X4) +error),
             cp_y = X1,
             yideal = case_when(age< X4 ~ X1 + X2 *(age-X4),
                                age > X4 ~X1 +  (X2 + X3) *(age-X4))
      )
    # Return data for Stan
    return(list(
      trueparms = c(
        "beta[1]"       = intercept.truemean,
        "beta[2]"       = preslope.truemean,
        "beta3"       = diffslope.truemean,
        "betacp"        =  cp.truemean,
        "u_sd[1]"       = intercept.truesd,
        "u_sd[2]"       = preslope.truesd,
        "u_sd[3]"       = diffslope.truesd,
        "u_sd[4]"       = cp.truesd,
        "y_sd"          = error.truesd
      ),
      standata = list(
        "N"             = nrow(simdat),
        "Npat"          = Npat,
        "cp_mean"   = cp.truemean,
        "cp_sd"   = cp.truesd,
        "zeros4"        = rep(0, 4),
        "betacp_lower"  = cp.lowerbound,
        "id"            = simdat$id,
        "x"           = simdat$age,
        "y"             = simdat$y,
        "cp"      = simdat$X4,
        "yideal"             = simdat$yideal,
        "cp_y"      = simdat$cp_y,
        "agemin"        = agemin,
        "agemax"        = agemax
      )
    ))
    
  }

#============================================================================
# Create simulated data, and fit the BAyesian Bent-Line Regression model
#============================================================================
# Create a single simulated dataset with specified seed
seed_choice= sample(1:1000000, 100)
seed <- seed_choice[2]
Npat<-100

standata <- createdata(seed = seed,Npat=Npat)

# Specify the location of the Stan model code file
stancode.filepath <- "BAyesian_Bent-Line_Regression_model.stan"
stan_code <- readChar(stancode.filepath, file.info(stancode.filepath)$size)
cat(stan_code)
# Specify which parameters we wish to monitor the MCMC samples for
# Notes: "Monitoring" a parameter means that the MCMC samples related to 
#        that parameter are retained, so they can be used later for post-
#        processing and inference. The default is for RStan to monitor all 
#        model parameters.
#        However to save memory, we wish to only monitor those parameters
#        which will be used for inference or checking model diagnostics.
#        The parameter names are listed in a vector, and the parameters
#        names must align with those that are used in the Stan model code
#        file (ie, the file "BAyesian_Bent-Line_Regression_model.stan").
#        Here, we choose to monitor the log posterior (lp__), all key 
#        model parameters, as well as prediction and random effect estimates
#        for a subset of the patients (see the Stan model code and the data  
#        steps above for more detail).
stanmonitor <- c("beta", 
                 "beta3",
                 "betacp",
                 "beta_sd",
                 "y_sd", 
                 "u_sd",
                 "alpha", 
                 "lp__")

# Fit the BAyesian Bent-Line Regression model to the simulated data
# Notes: Here we are fitting 4 chains in parallel, using 8 PC cores. The 
#        number of cores may need to be adjusted according to the PC being
#        used to fit the model. The control option is used to adjust some
#        of the parameters used for the Hamiltonian Monte Carlo sampler; 
#        see the RStan vignette (http://mc-stan.org/interfaces/rstan) for 
#        a discussion of these.

n_chains= 4
n_cores <- 8
total_iterations <- 1000
warmups <- 500
max_treedepth = 18
adapt_delta = 0.99

# 1. Compile and fit model.
# Start the clock
ptm <- proc.time()
stanfit <- stan(file    = stancode.filepath, 
                data    = standata$standata,
                pars    = stanmonitor, 
                chains  = n_chains, 
                cores   = n_cores, 
                iter    = total_iterations, 
                warmup  = warmups,
                seed = seed,
                control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))

# Stop the clock
proc.time() - ptm

#============================================================================
# Some diagnostics for the fitted model
#============================================================================
#####################################################
## Create a general theme that can be applied to both plots
theme_DLL <- theme_bw() +
  theme(panel.grid.major=element_line(color="gray90",linetype="dashed"),
        panel.grid.minor=element_blank(),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        legend.text = element_text(face="bold", size=14),
        legend.position=c(0,1),
        legend.justification=c(-0.05,1.02),
        legend.title=element_blank())
####################################################
parsinfo<-c("beta[1]","beta[2]","beta3","betacp","u_sd[1]","u_sd[2]","u_sd[3]","u_sd[4]","y_sd")
#load("stanfitsim_n100_10000iteration_example.Rda")
stan_dens(stanfit, pars= parsinfo, include = TRUE, unconstrain = FALSE,
          inc_warmup = FALSE, 
          separate_chains = TRUE,alpha=0.1)+
  scale_fill_manual(values = c("red", "blue", "green", "black"))+
  theme_DLL

stan_trace(stanfit, pars= parsinfo, include = TRUE, unconstrain = FALSE,
           inc_warmup = FALSE, 
           separate_chains = FALSE,
           alpha=0.6)+theme_DLL

stan_ac(stanfit, pars= parsinfo, include = TRUE, unconstrain = FALSE,
        inc_warmup = FALSE, 
        separate_chains = FALSE,
        alpha=0.6)+theme_DLL
#rhat
stanfit %>%
  bayesplot::rhat() %>%
  mcmc_rhat()  

print(stanfit, pars = c("beta", 
                        "beta3",
                        "betacp", 
                        "y_sd",
                        "u_sd",
                        "alpha"),digits=5)

