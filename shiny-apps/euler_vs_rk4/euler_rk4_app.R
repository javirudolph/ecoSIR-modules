library(shiny)
library(bslib)
source("euler_rk4_solver.R")


ui <- page_sidebar(
  title = "Euler's Method vs Runge-Kutta Method Simulation",
  sidebar = sidebar(
    helpText("Simulation to understand the difference between solving SIR model equations using Euler's method and the Fourth-Order Runge-Kutta Method using the function rk4sys in the package pracma"),
    sliderInput("beta",label="Transmission Rate",min=0,max=2,value=1.5,step=0.1),
    sliderInput("gamma",label="Recovery Rate",min=0,max=1,value=0.5,step=0.1),
    sliderInput("step",label = "Step Size",min=0,max=0.5,step=0.01,value=0.25)
  ),
  plotOutput("plot_euler"),
  plotOutput("plot_rk4"),
  plotOutput("plot_combine")
)


server <- function(input, output) {
  output$plot_euler <- renderPlot(
    
    euler_sir(S0=1-1e-6, I0 = 1e-6,0,beta=input$beta,gamma=input$gamma,len=50,dt=input$step)
  ) 
  output$plot_rk4 <- renderPlot(
    
    rk4_sir(S0=1-1e-6, I0 = 1e-6,0,b=input$beta,g=input$gamma,len=50,dt=input$step)
    
  )
  output$plot_combine <- renderPlot(
    euler_rk4_combined(S0=1-1e-6, I0 = 1e-6,0,beta=input$beta,gamma=input$gamma,len=50,dt=input$step)
    
  )
}

 
shinyApp(ui = ui, server = server)
