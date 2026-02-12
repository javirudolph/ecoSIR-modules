library(callr)
library(shiny)

apps <- list(
  list(path = "./shiny-apps/threshold_phenomenon/threshold_phenom_app.R", port = 4567),
  list(path = "./shiny-apps/eff_reproduction_sim/eff_reprod_app.R", port = 4568),
  list(path = "./shiny-apps/equilibrium/equilibrium_app.R", port = 4569),
  list(path = "./shiny-apps/euler_vs_deSolve/euler_desolve_app.R", port = 4570),
  list(path = "./shiny-apps/transmission_change_sir/trans_change_app.R", port = 4571)
)

bg_processes <- list()

for (app in apps) {
  cat("Launching:", app$path, "→ port", app$port, "\n")
  
  p <- r_bg(function(path, port) {
    shiny::runApp(path, port = port, launch.browser = FALSE)
  },
  args = list(path = app$path, port = app$port),
  supervise = TRUE)
  
  bg_processes[[as.character(app$port)]] <- p
}

