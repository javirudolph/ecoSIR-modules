library(callr)
library(here)

apps <- list(
  list(path = here("shiny-apps", "threshold_phenomenon", "threshold_phenom_app.R"), port = 4567),
  list(path = here("shiny-apps", "eff_reproduction_sim", "eff_reprod_app.R"), port = 4568),
  list(path = here("shiny-apps", "equilibrium", "equilibrium_app.R"), port = 4569),
  list(path = here("shiny-apps", "euler_vs_deSolve", "euler_desolve_app.R"), port = 4570),
  list(path = here("shiny-apps", "transmission_change_sir", "trans_change_app.R"), port = 4571)
)

for (app in apps) {
  log_file <- paste0("shiny_log_", gsub("[^a-zA-Z0-9_-]", "_", basename(app$path)), ".txt")
  file.create(log_file)  # create upfront
  
  p <- callr::r_bg(
    func = function(path, port, log_file) {
      # logging
      log <- file(log_file, open = "wt")
      sink(log)
      sink(log, type = "message")
      on.exit({ sink(type = "message"); sink(); close(log) })
      
      # required packages
      library(shiny)
      
      # absolute path
      path <- normalizePath(path, mustWork = TRUE)
      
      # run app
      shiny::runApp(path, port = port, launch.browser = FALSE)
    },
    args = list(app$path, app$port, log_file),
    supervise = TRUE
  )
  
  cat("Launched:", app$path, "→ logs at", log_file, "\n")
  cat("Process alive?", p$is_alive(), "\n")
}
