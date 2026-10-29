# Load necessary package
library(deSolve)

# Define the SIR model function
SIR <- function(time, state, parameters) {
  S <- state[1]  # Susceptible individuals
  I <- state[2]  # Infected individuals
  R <- state[3]  # Recovered individuals
  
  beta <- parameters["beta"]  # Transmission rate
  v <- parameters["v"]        # Recovery rate
  
  # Differential equations
  dS <- -beta * S * I
  dI <- beta * S * I - v * I
  dR <- v * I
  
  return(list(c(dS, dI, dR)))
}

# Set initial conditions
N <- 1000     # Total population
I0 <- 1       # Initial infected
S0 <- N - I0  # Initial susceptible (everyone else)
R0 <- 0       # Initial recovered
initial_state <- c(S = S0/N, I = I0/N, R = R0/N)  # Normalize by population size

# Set parameters
R0_value <- 3       # Basic reproduction number
v <- 1/20           # Recovery rate (1/average duration of infection)
beta <- R0_value * v  # Transmission rate calculated from R0

parameters <- c(beta = beta, v = v)

# Set time range
times <- seq(0, 160, by = 1)  # Simulate for 160 days

# Solve the differential equations
sir_out <- ode(y = initial_state, times = times, func = SIR, parms = parameters)

# Convert output to dataframe
sir_df <- as.data.frame(sir_out)

# Plot the results
library(ggplot2)

ggplot(data = sir_df, aes(x = time)) +
  geom_line(aes(y = S, color = "Susceptible"), size = 1) +
  geom_line(aes(y = I, color = "Infected"), size = 1) +
  geom_line(aes(y = R, color = "Recovered"), size = 1) +
  labs(title = "SIR Model Simulation", x = "Time (days)", y = "Proportion of Population") +
  scale_color_manual(values = c("Susceptible" = "blue", "Infected" = "red", "Recovered" = "green")) +
  theme_minimal()
