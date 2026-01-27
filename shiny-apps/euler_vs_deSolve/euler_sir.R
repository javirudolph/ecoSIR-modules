library(tidyverse)
library(deSolve)

euler_sir <- function(S0, I0, R0, beta, gamma, len, dt) {
  tbeg <- 0
  tend <- len/dt
  times <- seq(from=0, to=len, by=dt)
  S <- S0
  I <- I0
  R <- R0
  #I should try to create an empty matrix or dataframe here instead of three individual vectors
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

  df <- data.frame(times,S.vect,I.vect,R.vect) |> pivot_longer(cols=c(S.vect,I.vect,R.vect),names_to="n",values_to="v")
  ggplot(df,aes(x=times,y=v,color=n)) + geom_line(linewidth=1) + labs(title="SIR System Solved Using Euler's Method",x="Time (t)",y="Proportion n(t)",color="Population") + scale_color_discrete(labels = c("S.vect" = "Susceptible","I.vect" = "Infected","R.vect" = "Recovered"),limits=c("S.vect","I.vect","R.vect"))
}

desolve_sir <- function(S0, I0, R0, b, g, len, dt) {
  sir <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- -beta * S * I
      dI <- beta * S * I - gamma * I
      dR <- gamma * I
      
      return(list(c(dS, dI, dR)))
    })
  }
  init <- c(S=S0, I= I0, R= R0)
  ## Original parameters for the nice curve: parameters <- c(beta = 1.4247, gamma = 0.14286)
  parameters <- c(beta=b, gamma=g)
  times <- seq(0, len, by = dt)
  out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
  out <- out |> pivot_longer(cols=c("S","I","R"),names_to = "n",values_to = "v")
  ggplot(out,aes(x=time,y=v,color=n)) + geom_line(linewidth=1) + labs(title="SIR System Solved Using deSolve",x="Time (t)",y="Proportion n(t)",color="Population") + scale_color_discrete(labels = c("S" = "Susceptible","I" = "Infected","R" = "Recovered"),limits=c("S","I","R"))
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

desolve_sir_df <- function(S0, I0, R0, b, g, len, dt) {
  sir <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- -beta * S * I
      dI <- beta * S * I - gamma * I
      dR <- gamma * I
      
      return(list(c(dS, dI, dR)))
    })
  }
  init <- c(S=S0, I= I0, R= R0)
  ## Original parameters for the nice curve: parameters <- c(beta = 1.4247, gamma = 0.14286)
  parameters <- c(beta=b, gamma=g)
  times <- seq(0, len, by = dt)
  out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
  out <- out |> pivot_longer(cols=c("S","I","R"),names_to = "n",values_to = "v")
  out$method <- "deSolve"
  return(out)
}

euler_desolve_combined <- function(S0, I0, R0, beta, gamma, len, dt){
euler <- euler_sir_df(S0, I0,R0,beta,gamma,len,dt)

euler$n <- gsub("S\\.vect","S",euler$n)
euler$n <- gsub("I\\.vect","I",euler$n)
euler$n <- gsub("R\\.vect","R",euler$n)

desolve <- desolve_sir_df(S0, I0,R0,beta,gamma,len,dt)

combined_df <- rbind(euler,desolve)

ggplot(combined_df, aes(x = time, y = v, color = n, linetype = method)) + geom_line(linewidth = 1) + scale_color_manual(values = c("S" = "green", "I" = "red", "R" = "blue")) + scale_linetype_manual(values = c("Euler" = "dotted", "deSolve" = "solid")) + labs(title = "Euler vs deSolve", x = "Time (t)", y = "Proportion (n(t))", color = "Compartment", linetype = "Solving Method") + scale_color_discrete(labels = c("S" = "Susceptible","I" = "Infected","R" = "Recovered"),limits=c("S","I","R"))

}



