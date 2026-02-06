#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(bslib)
library(tidyverse)
library(plotly)
source("threshold_phenom_solver.R")

# Define UI for application that draws a histogram
ui <- page_sidebar(
  title = "Threshold Phenomenon Simulation",
  sidebar = sidebar(
    helpText("Simulation to help understand the threshold phenomenon. Double click on any of the curves in the legend to isolate it."),
    strong("Threshold Value (Rounded):"),
    textOutput("threshold"),
    sliderInput("beta",label="Transmission Rate",min=0,max=2,value=1.5,step=0.1),
    sliderInput("gamma",label="Recovery Rate",min=0,max=1,value=0.5,step=0.1),
    sliderInput("init_sus",label="Initial Proportion of Susceptibles",min=0,max=1,value=0.5,step=0.1)
  ),
  plotlyOutput("plot_sir"),
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$plot_sir <- renderPlotly({
    
    df <- solve_sir(S0=input$init_sus, I0 = 1e-6,R0=0,b=input$beta,g=input$gamma,len=100,dt=0.1)
    df_long <- pivot_longer(df, cols = c(S, I, R), names_to = "Compartment", values_to = "Value")
    df_long$Compartment <- factor(df_long$Compartment,levels = c("S", "I", "R"),labels = c("Susceptible", "Infected", "Recovered"))
    p <- ggplot(df_long,aes(x = time, y = Value, color = Compartment,group=Compartment,text = paste0("Time: ", round(time, 2),"<br>", Compartment, ": ", round(Value,6)))) +geom_line(linewidth = 0.6) +labs(x = "Time", y = "Proportion", color = "Population")
    ggplotly(p, tooltip = "text")
  })

  output$threshold <- renderText(round(input$gamma/input$beta,3))
}

# Run the application 
shinyApp(ui = ui, server = server)
