library(tidyverse)
library(pracma)

euler_sir <- function(S0, I0, R0, beta, gamma, len, dt) {
  tbeg <- 0
  tend <- len/dt
  times <- seq(from=0, to=len, by=dt)
  S <- S0
  I <- I0
  R <- R0

  S.vect <- rep(0, length(times))
  I.vect <- rep(0, length(times))
  R.vect <- rep(0, length(times))
  
  for(tt in tbeg:tend){
    S.vect[(tt+1)] <- S
    I.vect[(tt+1)] <- I
    R.vect[(tt+1)] <- R
    
    dS <- -beta * S * I * dt
    dI <- beta * S * I * dt - gamma * I * dt
    dR <- gamma * I * dt
    S <- S + dS
    I <- I + dI
    R <- R + dR
  }

  df <- data.frame(times,S.vect,I.vect,R.vect) |> pivot_longer(cols=c(S.vect,I.vect,R.vect),names_to="n" ,values_to="v")
  
  View(df)
  
  ggplot(df,aes(x=times,y=v,color=n)) + geom_line(linewidth=1) + labs(title="SIR System Solved Using Euler's Method",x="Time (t)",y="Proportion n(t)",color="Population") + scale_color_discrete(labels = c("S.vect" = "Susceptible","I.vect" = "Infected","R.vect" = "Recovered"),limits=c("S.vect","I.vect","R.vect"))
}

rk4_sir <- function(S0, I0, R0, b, g, len, dt) {
  ModelFunc <- function(time,state,parameters){
    S <- state[1]
    I <- state[2]
    R <- state[3]
    
    beta <- parameters["beta"]
    gamma <- parameters["gamma"]
    
    dS <- -(beta*I*S)
    dI <- (beta*I*S) - (gamma*I)
    dR <- (gamma*I)
    
    return(c(dS,dI,dR))
  }
  
  start <- 0
  end <- len
  initial_vals <- c(S=S0,I=I0,R=R0)
  model_parameters <- c("beta"=b,"gamma"=g)
  step_size <- as.integer((end-start)/dt)
  output <- rk4sys(f = function (t,y) ModelFunc(t,y,model_parameters),
                   a=start,
                   b=end,
                   y0=initial_vals,
                   n=step_size)
  
  output <- as.data.frame(output)
  out <- output |> pivot_longer(cols=c("y.S","y.I","y.R"),names_to = "n",values_to = "v")
  out <- rename(out,time=x)
  
  ggplot(out,aes(x=time,y=v,color=n)) + geom_line(linewidth=1) + labs(title="SIR System Solved Using rk4sys",x="Time (t)",y="Proportion n(t)",color="Population") + scale_color_discrete(labels = c("y.S" = "Susceptible","y.I" = "Infected","y.R" = "Recovered"),limits=c("y.S","y.I","y.R"))
}

euler_sir_df <- function(S0, I0, R0, beta, gamma, len, dt) {
  tbeg <- 0
  tend <- len/dt
  time <- seq(from=0, to=len, by=dt)
  S <- S0
  I <- I0
  R <- R0
  #I should try to create an empty matrix or dataframe here instead of three individual vectors
  S.vect <- rep(0, length(time))
  I.vect <- rep(0, length(time))
  R.vect <- rep(0, length(time))
  
  for(tt in tbeg:tend){
    S.vect[(tt+1)] <- S
    I.vect[(tt+1)] <- I
    R.vect[(tt+1)] <- R
    
    dS <- -beta * S * I * dt
    dI <- beta * S * I * dt - gamma * I * dt
    dR <- gamma * I * dt
    S <- S + dS
    I <- I + dI
    R <- R + dR
  }
  
  df <- data.frame(time,S.vect,I.vect,R.vect) |> pivot_longer(cols=c(S.vect,I.vect,R.vect),names_to="n",values_to="v")
  df$method <- "Euler"
  return(df)
}

rk4_sir_df <- function(S0, I0, R0, b, g, len, dt) {
  ModelFunc <- function(time,state,parameters){
    S <- state[1]
    I <- state[2]
    R <- state[3]
    
    beta <- parameters["beta"]
    gamma <- parameters["gamma"]
    
    dS <- -(beta*I*S)
    dI <- (beta*I*S) - (gamma*I)
    dR <- (gamma*I)
    
    return(c(dS,dI,dR))
  }
  
  start <- 0
  end <- len
  initial_vals <- c(S=S0,I=I0,R=R0)
  model_parameters <- c("beta"=b,"gamma"=g)
  step_size <- as.integer((end-start)/dt)
  
  output <- rk4sys(f = function (t,y) ModelFunc(t,y,model_parameters),
                   a=start,
                   b=end,
                   y0=initial_vals,
                   n=step_size)
  
  output <- as.data.frame(output)
  out <- output |> pivot_longer(cols=c("y.S","y.I","y.R"),names_to = "n",values_to = "v")
  out <- rename(out,time=x)
  out$method <- "rk4sys"
  return(out)
}


euler_rk4_combined <- function(S0, I0, R0, beta, gamma, len, dt){
  euler <- euler_sir_df(S0, I0,R0,beta,gamma,len,dt)
  
  euler$n <- gsub("S\\.vect","S",euler$n)
  euler$n <- gsub("I\\.vect","I",euler$n)
  euler$n <- gsub("R\\.vect","R",euler$n)
  
  rk4 <- rk4_sir_df(S0, I0,R0,beta,gamma,len,dt)
  rk4$n <- gsub("y\\.S","S",rk4$n)
  rk4$n <- gsub("y\\.I","I",rk4$n)
  rk4$n <- gsub("y\\.R","R",rk4$n)
  
  combined_df <- rbind(euler,rk4)
  
  View(combined_df)
  
  ggplot(combined_df, aes(x = time, y = v, color = n, linetype = method)) + geom_line(linewidth = 1) + scale_color_manual(values = c("S" = "green", "I" = "red", "R" = "blue")) + scale_linetype_manual(values = c("Euler" = "dotted", "rk4sys" = "solid")) + labs(title = "Euler's method vs rk4sys", x = "Time (t)", y = "Proportion (n(t))", color = "Compartment", linetype = "Solving Method") + scale_color_discrete(labels = c("S" = "Susceptible","I" = "Infected","R" = "Recovered"),limits=c("S","I","R"))
  
}



