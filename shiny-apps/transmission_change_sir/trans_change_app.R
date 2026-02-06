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
source("trans_change_solver.R")

# Define UI for application that draws a histogram
ui <- page_sidebar(
  title = "Change in Transmission Simulation",
  sidebar = sidebar(
    helpText("Simulation to understand how the SIR curves change when the transmission rate change (Note: The transmission rate in this model is calculated directly as the product of c (the per capita contact rate per unit time t) and p (the probability of infection for each susceptible-infected contact)."),
    sliderInput("contact",label="Per Capita Contact Rate (c)",min=0,max=2,value=1.5,step=0.1),
    sliderInput("prob_inf",label="Probability of Infection for Susceptible-Infected Contacts (p)",min=0,max=1,value=0.5,step=0.1),
    sliderInput("gamma",label="Recovery Rate",min=0,max=1,value=0.5,step=0.1)
  ),
  plotOutput("plot_sir"),
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$plot_sir <- renderPlot(
    solve_sir(S0=1-1e-6, I0 = 1e-6,R0=0,c=input$contact,p=input$prob_inf,g=input$gamma,len=50,dt=0.1)
  ) 
}

# Run the application 
shinyApp(ui = ui, server = server)
