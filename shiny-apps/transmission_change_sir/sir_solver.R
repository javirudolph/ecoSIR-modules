#function to solve sir model using deSolve package
#beta not input, contact rate and probability of infection per contact instead

solve_sir <- function(S0, I0, R0, c, p, g, len, dt) {
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
  parameters <- c(beta=c*p, gamma=g)
  times <- seq(0, len, by = dt)
  out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
  out <- out |> pivot_longer(cols=c("S","I","R"),names_to = "n",values_to = "v")
  ggplot(out,aes(x=time,y=v,color=n)) + geom_line(linewidth=1) + labs(x="Time (t)",y="Proportion n(t)",color="Population") + scale_color_discrete(labels = c("S" = "Susceptible","I" = "Infected","R" = "Recovered"),limits=c("S","I","R"))
}