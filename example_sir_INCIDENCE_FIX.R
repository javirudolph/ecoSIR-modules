library(bbmle)
library(pracma)
library(dplyr)

# Importing our prevalence count data for Blackburn
data <- read.csv("https://raw.githubusercontent.com/javirudolph/ecoSIR-modules/refs/heads/main/data/measles_ex_data.csv")
data <- data[data$city=="Blackburn",]

# Extracting population size for later
pop_size <- data$population_size[1]

# Truncating the data to be only entries after the first infection appears
first_infect <- which(data$prevalence_count!=0)[1]
data_trunc <- data[3:38,]

# SIR model function
ModelFunc <- function(time,state,parameters){
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  beta <- parameters["beta"]
  gamma <- parameters["gamma"]
  N <- parameters["N"]
  
  dS <- -(beta*I*S/N)
  dI <- (beta*I*S/N) - (gamma*I)
  dR <- (gamma*I)
  
  return(c(dS,dI,dR))
}

# Outputs predicted prevalence counts
pred <- function(parameters,t){
  # S = population size - prevalence at first non-zero prevalence
  # I = first non-zero prevalence count
  # R = 0
  initial_vals <- c(S=pop_size - data_trunc$prevalence_count[1],
                    I=data_trunc$prevalence_count[1],
                    R=0)
  
  t_start <- min(t)
  t_end <- max(t)
  step_size <- (t_end - t_start)*500
  
  output <- rk4sys(f=function(t,y) ModelFunc(t,y,parameters),
                   a=t_start,
                   b=t_end,
                   y0=initial_vals,
                   n=step_size)
  
  # Extracting predicted prevalences
  output_df <- as.data.frame(output)
  prevalences <- filter(output_df,t_start <= x, t_end >= x) |> select(x,y.S,y.I)
  
  # Ensuring there are no NaNs or out of bounds numbers in the output
  prevalences[is.nan(prevalences$y.I) | prevalences$y.I < 0 ] <- NA

  return(prevalences)
}

find_incidence <- function(solutions){
  # Not estimating the first observation (week 6)
  # We do this because there are no values before it to integrate over
  incidence_times <- data_trunc$week[-1]
  incidence_times <- sapply(incidence_times,
                            incidence_func,
                            solutions = solutions)
  return(incidence_times)
}

incidence_func <- function(t_start,solutions){
  # Extracting the values needed to do trapezoidal integration
  times <- filter(solutions, x >= t_start-2, solutions$x <= t_start) |> pull(x)
  susceptibles <- filter(solutions, x >= t_start-2, solutions$x <= t_start) |> pull(y.S)
  infectious <- filter(solutions, x >= t_start-2, solutions$x <= t_start) |> pull(y.I)
  incidence_steps <- susceptibles*infectious/pop_size
  
  # Numerically integrating using the trapezoidal rule
  incidence <- trapz(times,incidence_steps)
  
  # Returning the incidence at t_start
  return(incidence)
}

# Likelihood function
poisson_ll <- function(b){
  times <- data_trunc$week
  parameters <- c(beta=b,gamma=1.3,N=pop_size)
  predicted <- pred(parameters,times)
  incidence <- find_incidence(predicted)
  print(b)
  print(incidence)
  sum(dpois(x=data_trunc$prevalence_count[-1],lambda=incidence,log=TRUE))
}

# Setting our initial guess as beta = 2
start <- list(b=4)

beta_seq <- seq(1,5,by=0.1)
neg_ll <- sapply(beta_seq,poisson_ll)
plot(beta_seq,neg_ll,type="l",xlab="beta",ylab="negative log-likelihood")

beta_seq[which.max(neg_ll)]

