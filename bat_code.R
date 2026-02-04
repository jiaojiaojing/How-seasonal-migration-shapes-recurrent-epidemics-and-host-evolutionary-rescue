# 1) Packages
library(deSolve)

# 2) Define parameters
parameters2 <- list(
  K       = 2000,   # carrying capacity
  r_SW    = 0.5,    # growth rate of Wild type of susceptible
  r_IW    = 0.05,   # growth rate of Wild type infected
  r_RW    = 0.4,    # growth rate of Wild type recovered
  r_ST    = 0.4,    # growth rate of Robust type susceptible
  r_IT    = 0.04,   # growth rate of Robust type infected
  r_RT    = 0.32,    # growth rate of Robust type recovered
  mu_1    = 0.1,    # natural mortality
  alpha_W = 0.05,   # disease-induced mortality in Wild type
  alpha_T = 0.01,   # disease-induced mortality in Robust type
  beta_WT = 0.00005, # transmission rate from infected W to susceptible R
  beta_TW = 0.00035,   # transmission rate from infected R to susceptible W
  beta_WW = 0.00035,   # transmission rate from infected W to susceptible W
  beta_TT = 0.00005, # transmission rate from infected R to susceptible R
  gamma_W = 0.01,
  gamma_R = 0.01,
  yita_S  = 20,      # coefficient of seasonal migration in susceptible
  yita_I  = 20,      # coefficient of seasonal migration in infected
  yita_R  = 20       # coefficient of seasonal migration in recovered
)


# 3) Define times and initial states
Timesteps <- 1000            # total time steps
times     <- seq(0, Timesteps, by = 1)

# Initial conditions: [S_W, S_T, I_W, I_T, R_W, R_T]
inits <- c(
  S_W = 100,
  S_T = 100,
  I_W = 10,
  I_T = 0,
  R_W = 0,
  R_T = 0
)


# 4) Define the signal function
y_vec <- sin(2 * pi * times / (365 / 1))

signal <- approxfun(x = times, y = y_vec, method = "linear", rule = 2)


# 5) Define a helper function check_neg
check_neg <- function(old_value, change) {
  if (old_value + change < 5) {
    return(0)
  } else {
    return(change)
  }
}

# 6) Define the ODE 
bat_function <- function(t, state, parameters) {
  
  # Unpack state variables
  S_W <- state[1]
  S_T <- state[2]
  I_W <- state[3]
  I_T <- state[4]
  R_W <- state[5]
  R_T <- state[6]
  
  # Unpack parameters
  with(as.list(parameters), {
    
    # Compute total population
    total_pop <- S_W + S_T + I_W + I_T + R_W + R_T
    
    # For safety, handle total_pop = 0 (avoid division by zero)
    if(total_pop <= 0) {
      return(list(c(0,0,0,0,0,0)))
    }
    
    # Calculate intermediate "birth" or growth terms (similar to your original approach)
    LS_WS_W <- max(S_W * S_W / total_pop * r_SW * (1 - total_pop / K), 0)
    LS_WS_R <- max(S_W * S_T / total_pop * (r_SW + r_ST)/2 * (1 - total_pop / K), 0)
    LS_RS_R <- max(S_T * S_T / total_pop * r_ST * (1 - total_pop / K), 0)
    
    LI_WI_W <- max(I_W * I_W / total_pop * r_IW * (1 - total_pop / K), 0)
    LI_WI_R <- max(I_W * I_T / total_pop * (r_IW + r_IT)/2 * (1 - total_pop / K), 0)
    LI_RI_R <- max(I_T * I_T / total_pop * r_IT * (1 - total_pop / K), 0)
    
    LR_WR_W <- max(R_W * R_W / total_pop * r_RW * (1 - total_pop / K), 0)
    LR_WR_R <- max(R_W * R_T / total_pop * (r_RW + r_RT)/2 * (1 - total_pop / K), 0)
    LR_RR_R <- max(R_T * R_T / total_pop * r_RT * (1 - total_pop / K), 0)
    
    LS_WI_W <- max(2 * S_W * I_W / total_pop * (r_SW + r_IW)/2 * (1 - total_pop / K), 0)
    LS_WI_R <- max(S_W * I_T / total_pop * (r_SW + r_IT)/2 * (1 - total_pop / K), 0)
    LS_RI_R <- max(2 * S_T * I_T / total_pop * (r_ST + r_IT)/2 * (1 - total_pop / K), 0)
    LS_RI_W <- max(S_T * I_W / total_pop * (r_ST + r_IW)/2 * (1 - total_pop / K), 0)
    
    LI_WR_W <- max(2 * I_W * R_W / total_pop * (r_IW + r_RW)/2 * (1 - total_pop / K), 0)
    LI_WR_R <- max(I_W * R_T / total_pop * (r_IW + r_RT)/2 * (1 - total_pop / K), 0)
    LI_RR_R <- max(2 * I_T * R_T / total_pop * (r_IT + r_RT)/2 * (1 - total_pop / K), 0)
    LI_RR_W <- max(I_T * R_W / total_pop * (r_IT + r_RW)/2 * (1 - total_pop / K), 0)
    
    LS_WR_W <- max(2 * S_W * R_W / total_pop * (r_SW + r_RW)/2 * (1 - total_pop / K), 0)
    LS_WR_R <- max(S_W * R_T / total_pop * (r_SW + r_RT)/2 * (1 - total_pop / K), 0)
    LS_RR_R <- max(2 * S_T * R_T / total_pop * (r_ST + r_RT)/2 * (1 - total_pop / K), 0)
    LS_RR_W <- max(S_T * R_W / total_pop * (r_ST + r_RW)/2 * (1 - total_pop / K), 0)
    
    # Sum terms for "birth" or net increase to each genotype
    B_W <- LS_WS_W + LS_WS_R + LI_WI_W + LI_WI_R + LR_WR_W + LR_WR_R +
      LS_WI_W + LS_WI_R + LS_RI_W + LI_WR_W + LI_WR_R + LI_RR_W +
      LS_WR_W + LS_WR_R + LS_RR_W
    
    B_T <- LS_RS_R + LS_WS_R + LI_RI_R + LI_WI_R + LR_RR_R + LR_WR_R +
      LS_RI_R + LS_WI_R + LS_RI_W + LI_RR_R + LI_WR_R + LI_RR_W +
      LS_RR_R + LS_WR_R + LS_RR_W
    
    # Disease transmission + migration + mortality
    # dS_W
    tmp1 <- B_W -
      beta_WW * S_W * I_W -
      beta_WT * S_W * I_T -
      mu_1 * S_W +
      yita_S * signal(t)
    
    dS_W <- check_neg(S_W, tmp1)
    
    # dS_R
    tmp2 <- B_T -
      beta_TT * S_T * I_T -
      beta_TW * S_T * I_W -
      mu_1 * S_T +
      yita_S * signal(t)
    
    dS_T <- check_neg(S_T, tmp2)
    
    # dI_W
    tmp3 <- beta_WW * S_W * I_W +
      beta_WT * S_W * I_T -
      (mu_1 + alpha_W) * I_W -
      gamma_W * I_W +
      yita_I * signal(t)
    
    dI_W <- check_neg(I_W, tmp3)
    
    # dI_T
    tmp4 <- beta_TT * S_T * I_T +
      beta_TW * S_T * I_W -
      (mu_1 + alpha_T) * I_T -
      gamma_R * I_T +
      yita_I * signal(t)
    
    dI_T <- check_neg(I_T, tmp4)
    
    # dR_W
    tmp5 <- gamma_W * I_W -
      mu_1 * R_W +
      yita_R * signal(t)
    
    dR_W <- check_neg(R_W, tmp5)
    
    # dR_T
    tmp6 <- gamma_R * I_T -
      mu_1 * R_T +
      yita_R * signal(t)
    
    dR_T <- check_neg(R_T, tmp6)
    
    # Return the list of derivatives
    list(c(dS_W, dS_T, dI_W, dI_T, dR_W, dR_T))
  })
}

# 7) Solve the system
out <- ode(
  y       = inits,
  times   = times,
  func    = bat_function,
  parms   = parameters2,
  method  = "ode45"   
)

# 8) Post-processing & Plotting
out_df <- as.data.frame(out)

# Compute totals
out_df$total_W <- out_df$S_W + out_df$I_W + out_df$R_W
out_df$total_R <- out_df$S_T + out_df$I_T + out_df$R_T
out_df$Prevalence <- ( out_df$I_W) / (out_df$S_W + out_df$S_T + out_df$I_W + out_df$I_T + out_df$R_W + out_df$R_T)

# Basic plots in base R
op <- par(mfrow = c(2, 3))  # 1x3 layout

# Susceptible
plot(out_df$time, out_df$S_W, type="l", col="red", lwd=2,
     xlab="Time", ylab="Number of Individuals", 
     main="Susceptible Individuals")
lines(out_df$time, out_df$S_T, col="blue", lwd=2)


# Infected
plot(out_df$time, out_df$I_W, type="l", col="red", lwd=2,
     xlab="", ylab="",
     main="",
     ylim = c(0, max(c(out_df$total_W, out_df$total_R), na.rm=TRUE)),
     cex.axis = 1.8, 
     cex.lab = 1.8,   
     cex.main = 2) 

lines(out_df$time, out_df$I_T, col="blue", lwd=2)


# Recovered
plot(out_df$time, out_df$R_W, type="l", col="red", lwd=2,
     xlab="Time", ylab="Number of Individuals",
     main="Recovered Individuals")
lines(out_df$time, out_df$R_T, col="blue", lwd=2)


# Total population
plot(out_df$time, out_df$total_W, type="l", col="red", lwd=2,
     xlab="Time", ylab="Number of Individuals",
     main="Total Population",
     ylim = c(0, 1500),
     cex.axis = 1.5,  # Makes axis numbers (e.g., 0, 200, 400) 50% larger
     cex.lab = 1.5,   # Makes axis labels ("Time", "Individuals") 50% larger
     cex.main = 2) # <- This line sets y-axis from 0 to 1000

lines(out_df$time, out_df$total_R, col="blue", lwd=2)




##################### Scenario 1 parameters#####################
parameters2 <- list(
  K       = 2000,   # carrying capacity
  r_SW    = 0.5,    # growth rate of Wild type of susceptible
  r_IW    = 0.05,   # growth rate of Wild type infected
  r_TW    = 0.4,    # growth rate of Wild type recovered
  r_ST    = 0.4,    # growth rate of Robust type susceptible
  r_IT    = 0.04,   # growth rate of Robust type infected
  r_RT    = 0.32,    # growth rate of Robust type recovered
  mu_1    = 0.1,    # natural mortality
  alpha_W = 0.03,   # disease-induced mortality in Wild type
  alpha_R = 0.01,   # disease-induced mortality in Robust type
  beta_WR = 0.00005, # transmission rate from infected W to susceptible R
  beta_TW = 0.00035,   # transmission rate from infected R to susceptible W
  beta_WW = 0.00035,   # transmission rate from infected W to susceptible W
  beta_TT = 0.00005, # transmission rate from infected R to susceptible R
  gamma_W = 0.01,    # recovery rate for Wild type
  gamma_R = 0.01,    # recovery rate for Wild type
  yita_S  = 0,      # coefficient of seasonal migration in susceptible
  yita_I  = 0,      # coefficient of seasonal migration in infected
  yita_R  = 0       # coefficient of seasonal migration in recovered
)
#### for scenario 1, yita_S/yita_I/yita_R have values equal to 0/10/25


##################### Scenario 2 parameters#####################
parameters2 <- list(
  K       = 2000,   # carrying capacity
  r_SW    = 0.5,    # growth rate of Wild type of susceptible
  r_IW    = 0.05,   # growth rate of Wild type infected
  r_TW    = 0.4,    # growth rate of Wild type recovered
  r_ST    = 0.4,    # growth rate of Robust type susceptible
  r_IT    = 0.04,   # growth rate of Robust type infected
  r_RT    = 0.32,    # growth rate of Robust type recovered
  mu_1    = 0.1,    # natural mortality
  alpha_W = 0.03,   # disease-induced mortality in Wild type
  alpha_ = 0.01,   # disease-induced mortality in Robust type
  beta_WR = 0.00005, # transmission rate from infected W to susceptible R
  beta_TW = 0.00035,   # transmission rate from infected R to susceptible W
  beta_WW = 0.00035,   # transmission rate from infected W to susceptible W
  beta_TT = 0.00005, # transmission rate from infected R to susceptible R
  gamma_W = 0.01,    # recovery rate for Wild type
  gamma_R = 0.01,    # recovery rate for Wild type
  yita_S  = 20,      # coefficient of seasonal migration in susceptible
  yita_I  = 20,      # coefficient of seasonal migration in infected
  yita_R  = 20       # coefficient of seasonal migration in recovered
)

### in scenario 2, yita is fixed to 20 but the step 4 have different values
y_vec <- sin(2 * pi * times / (365 / 1))
### based on the disease frequency, we tested 365/0.5; 365/1; 365/2


##################### Scenario 3 parameters#####################
### scenario 3 is a combined situation, we manipulated the parameters in scenario 1&2 to make a 3*3 grid.we set three different values for yita (5/15/30). And remains the disease frequency in scenario 2 (0.5/1/2). All the other parameters stay the same.

