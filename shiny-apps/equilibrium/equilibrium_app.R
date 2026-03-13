library(shiny)
library(bslib)
library(tidyverse)
library(plotly)
source("equilibrium_solver.R")

ui <- page_sidebar(
  title = "Equilibrium Simulation",
  sidebar = sidebar(
    helpText(HTML("Hover over the curves to see more information. Double click on any of the legend labels to isolate that curve. <br> <br> <b>Note:</b> The x-axis limits are dynamic.")),
    textOutput("reprod"),
    sliderInput("beta",label="Transmission Rate",min=0.1,max=2,value=1.5,step=0.1),
    sliderInput("gamma",label="Recovery Rate",min=0.1,max=1,value=0.5,step=0.1),
    sliderInput("mu",label="Birth / Death Rate",min=0.1,max=1,value=0.5,step=0.1),
    sliderInput("initial_inf",label="Initial Proportion of Infected Individuals",min=0.01,max=1,value=0.5,step=0.01)
    ),
  plotlyOutput("plot_sir"),
)

server <- function(input, output) {

     output$plot_sir <- renderPlotly({
      
      df <- solve_sir(S0=1-input$initial_inf, I0 = input$initial_inf,R0=0,b=input$beta,g=input$gamma,m=input$mu,len=1000,dt=0.1)
      df_long <- pivot_longer(df, cols = c(S, I, R), names_to = "Compartment", values_to = "Value")
      df_long$Compartment <- factor(df_long$Compartment,levels = c("S", "I", "R"),labels = c("Susceptible", "Infected", "Recovered"))
      
      df$I_diff <- c(0, diff(df$I) / diff(df$time))
      df$S_diff <- c(0, diff(df$S) / diff(df$time))
      df$R_diff <- c(0, diff(df$R) / diff(df$time))
      
      threshold <- 1e-5
      t_min <- min(df$time[abs(df$I_diff) > threshold | abs(df$S_diff) > threshold | abs(df$R_diff) > threshold])
      t_max <- max(df$time[abs(df$I_diff) > threshold | abs(df$S_diff) > threshold | abs(df$R_diff) > threshold])
      
      p <- ggplot(df_long,aes(x = time, y = Value, color = Compartment,group=Compartment,text = paste0("Time: ", round(time, 2),"<br>", Compartment, ": ", round(Value,6),"<br>Effective Reproduction Number: ", round(Re, 3)))) +geom_line(linewidth = 0.6) +labs(x = "Time", y = "Proportion", color = "Compartment")+ggtitle(paste0("Endemic SIR (R₀ = ",round(input$beta / (input$gamma + input$mu),3),")")) + xlim(t_min,t_max)
      ggplotly(p, tooltip = "text")
    })
    
  }
  

shinyApp(ui = ui, server = server)
