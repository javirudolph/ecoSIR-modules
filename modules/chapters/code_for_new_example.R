pop_size <- 150000

data <- read.csv("./modules/data/disease_flu_2017.csv")
data<-mutate(data,
  total_sum = rowSums(select(data, -WEEK, -TOTAL.SPECIMENS,-YEAR))
)
data<-mutate(data,WEEK=((YEAR-2017)*52)+WEEK)
data <- select(data,WEEK,total_sum)
data <- mutate(data,WEEK=WEEK-40)
data <- rename(data,day=WEEK,prevalence_count=total_sum)
sum(data$prevalence_count)

ModelFunc <- function(time,state,parameters){
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  beta <- parameters["beta"]
  gamma <- parameters["gamma"]
  N <- parameters["N"]
  
  dS <- -(beta*I*S)/N
  dI <- (beta*I*S)/N - (gamma*I)
  dR <- (gamma*I)
  
  return(c(dS,dI,dR))
}

# Function Inputs: Named vector of model parameters and a vector of times
# Function Outputs: Predicted prevalences at input times
pred <- function(parameters,t){
  initial_vals <- c(S=parameters["N"]-data[1,2],I=data[1,2],R=0)
  t_start <- 0
  t_end <- max(t)
  step_size <- t_end - t_start
  
  # Calling rk4sys to solve our model equations
  output <- rk4sys(f=function(t,y) ModelFunc(t,y,parameters),
                   a=t_start,
                   b=t_end,
                   y0=initial_vals,
                   n=step_size)
  
  # Extracting predicted prevalences
  output_df <- as.data.frame(output)
  prevalences <- output_df$y.I[t]
  
  # Ensuring there are no NaNs or out of bounds numbers
  prevalences[is.nan(prevalences) | prevalences < 0 ] <- 0
  
  return(prevalences)
}

# Function Inputs: A named list containing a value for beta
# Function Outputs: The log-likelihood of the data given the input beta
binomial_ll <- function(b){
  times <- data$day
  parameters <- c(beta=b,gamma=1.4,N=pop_size)
  predicted <- pred(parameters,times)
  
  # Returns the negative log-likelihood
  -sum(dbinom(x=data$prevalence_count,
              size = parameters["N"], 
              prob = predicted/parameters["N"],
              log=TRUE))
}

start <- list(b=2)

# Using the L-BFGS-B algorithm / method for MLE 
# Setting the lower bound to be small but positive since beta is non-negative
result <- mle2(minuslogl=binomial_ll,
               start=start,
               method="L-BFGS-B",
               lower=c(b=0.1),upper=c(b=3))

# Printing our estimate for beta
Model <- function(time,state){
  # Extracting the state variables: susceptible, infected, and recovered (S,I,R)
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  # Specifying the model parameters
  # We only have two parameters - the transmission and recovery rate
  # These parameters will change depending on the model
  beta <- unname(coef(result))[1]
  gamma <- 1.4
  
  # Constructing our SIR model equations
  dS <- -(beta*I*S)/pop_size
  dI <- (beta*I*S)/pop_size - (gamma*I)
  dR <- (gamma*I)
  
  # Now we will return those model equations as a vector
  return(c(dS,dI,dR))
}


# The start and ending points of the solution interval
start <- 0
end <- 51

# Initial values to be passed into the parameter y
# Make sure values line up with the state variables in the model function
initial_vals <- c(S=pop_size-171,I=171,R=0)

# How many steps to take from a to b
step_size <- 100

library(pracma)

# Make sure the values we defined are being passed into the correct arguments
output <- rk4sys(f=Model,
                 a=start,
                 b=end,
                 y0=initial_vals,
                 n=step_size)

# Converting the output matrix to a data frame and printing it
output <- as.data.frame(output)

times <- output$x
S <- output$y.S
I <- output$y.I
R <- output$y.R
beta <- 1.15
gamma <- 1.4
View(output)
plot(times, I, type = "l", col = "steelblue", lwd = 2, xlab = "Time", ylab = "Proportion",
     main = paste0("SIR Epidemic (Runge-Kutta)  |  β=", beta,
                   "  γ=", gamma, "  R₀=", round(beta/gamma, 2)))

points(data[1:nrow(data),]$day,data[1:nrow(data),]$prevalence_count)

unname(coef(result))[1]

write.csv(data,"disease_flu_2017_MODIFIED.csv")
writ