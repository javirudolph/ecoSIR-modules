library(shiny)
library(bslib)
library(tidyverse)
library(plotly)
source("equilibrium_solver.R")

# Define UI for application that draws a histogram
ui <- page_sidebar(
  title = "Effective Reproduction Number Simulation",
  sidebar = sidebar(
    helpText("Simulation to help understand disease-free and endemic equilibriums and their relationship to the basic reproduction number and initial proportion of infected individuals. Hover over each of the curves to see the effective reproduction number at that point. Double click on any of the curves in the legend to isolate it. "),
    textOutput("reprod"),
    sliderInput("beta",label="Transmission Rate",min=0,max=2,value=1.5,step=0.1),
    sliderInput("gamma",label="Recovery Rate",min=0,max=1,value=0.5,step=0.1),
    sliderInput("mu",label="Birth / Death Rate",min=0,max=1,value=0.5,step=0.1),
    sliderInput("initial_inf",label="Initial Proportion of Infected Individuals",min=0,max=1,value=0.5,step=0.01)
    ),
  plotlyOutput("plot_sir"),
)

# Define server logic required to draw a histogram
server <- function(input, output) {

     output$plot_sir <- renderPlotly({
      
      df <- solve_sir(S0=1-input$initial_inf, I0 = input$initial_inf,R0=0,b=input$beta,g=input$gamma,m=input$mu,len=500,dt=0.1)
      df_long <- pivot_longer(df, cols = c(S, I, R), names_to = "Compartment", values_to = "Value")
      df_long$Compartment <- factor(df_long$Compartment,levels = c("S", "I", "R"),labels = c("Susceptible", "Infected", "Recovered"))
      p <- ggplot(df_long,aes(x = time, y = Value, color = Compartment,group=Compartment,text = paste0("Time: ", round(time, 2),"<br>", Compartment, ": ", round(Value,6),"<br>Effective Reproduction Number: ", round(Re, 3)))) +geom_line(linewidth = 0.6) +labs(x = "Time", y = "Proportion", color = "Population")
      ggplotly(p, tooltip = "text")
    })
     output$reprod <- renderText({
       R0 <- input$beta / (input$gamma + input$mu)
       paste("Basic reproduction number:", round(R0, 3))
     })
    
  }
  

# Run the application 
shinyApp(ui = ui, server = server)
