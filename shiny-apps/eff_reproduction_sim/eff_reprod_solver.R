library(deSolve)

solve_sir <- function(S0, I0, R0, b, g, len, dt) {
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
  out$Re <- (b/g)*out$S
  return(out)
}

