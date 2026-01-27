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
setwd("C:/Users/anjal/Documents/euler_vs_deSolve")
source("euler_sir.R")

# Define UI for application that draws a histogram
ui <- page_sidebar(
  title = "Euler's Method vs DeSolve Simulation",
  sidebar = sidebar(
    helpText("Simulation to understand the difference between solving SIR model equations using Euler's method and the deSolve package"),
    sliderInput("beta",label="Transmission Rate",min=0,max=2,value=1.5,step=0.1),
    sliderInput("gamma",label="Recovery Rate",min=0,max=1,value=0.5,step=0.1),
    sliderInput("step",label = "For Euler's Method: Step Size",min=0,max=0.5,step=0.01,value=0.25)
  ),
  plotOutput("plot_euler"),
  plotOutput("plot_deSolve"),
  plotOutput("plot_combine")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$plot_euler <- renderPlot(
    
    euler_sir(S0=1-1e-6, I0 = 1e-6,0,beta=input$beta,gamma=input$gamma,len=50,dt=input$step)
  ) 
  
  output$plot_deSolve <- renderPlot(
    
    desolve_sir(S0=1-1e-6, I0 = 1e-6,0,b=input$beta,g=input$gamma,len=50,dt=input$step)
    
  )
  output$plot_combine <- renderPlot(
    euler_desolve_combined(S0=1-1e-6, I0 = 1e-6,0,beta=input$beta,gamma=input$gamma,len=50,dt=input$step)
    
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
