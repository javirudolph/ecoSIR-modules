#WARNING: YES I did have AI code this script. This is the only script I have used AI to code since it is just a rote process

# Load package
library(deSolve)

set.seed(123)  # reproducibility

# ----------------------------
# 1. Define SIR model
# ----------------------------
sir_model <- function(t, state, parameters) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  with(as.list(parameters), {
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    
    list(c(dS, dI, dR))
  })
}

# ----------------------------
# 2. Parameters
# ----------------------------
N <- 1000
parameters <- c(beta = 1.2,
                gamma = 0.7,
                N = N)

# Initial state
init_state <- c(S = 999,
                I = 1,
                R = 0)

# ----------------------------
# 3. Integer-day timeline
# ----------------------------
times <- 0:35   # days 0,1,2,...,20

# Solve SIR
out <- ode(y = init_state,
           times = times,
           func = sir_model,
           parms = parameters)

out_df <- as.data.frame(out)

# True infectious counts
I_true <- out_df$I

# ----------------------------
# 4. Poisson observation noise
# ----------------------------
I_observed <- rpois(length(I_true), lambda = I_true)

# ----------------------------
# 5. Dataset
# ----------------------------
data <- data.frame(
  day = times,
  I_true = round(I_true, 2),
  I_observed = I_observed
)

print(data)
library(tidyverse)
write.csv(select(data,day,I_observed),"mle_example.csv")
