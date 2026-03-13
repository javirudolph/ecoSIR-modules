library(shiny)
library(bslib)
library(tidyverse)
library(plotly)
source("eff_reprod_solver.R")

ui <- page_sidebar(
  title = "Effective Reproduction Number Simulation",
  sidebar = sidebar(
    helpText(HTML("Hover over the curves to see more information. Double click on any of the legend labels to isolate that curve. <br> <br> <b>Note:</b> The x-axis limits are dynamic.")),
    textOutput("reprod"),
    sliderInput("beta",label="Transmission Rate",min=0.1,max=2,value=1.5,step=0.1),
    sliderInput("gamma",label="Recovery Rate",min=0.1,max=2,value=0.5,step=0.1),
  ),
  plotlyOutput("plot_sir"),
)

server <- function(input, output) {

     output$plot_sir <- renderPlotly({
       
      df <- solve_sir(S0=1-1e-6, I0 = 1e-6,R0=0,b=input$beta,g=input$gamma,len=500,dt=0.1)
      df_long <- pivot_longer(df, cols = c(S, I, R), names_to = "Compartment", values_to = "Value")
      df_long$Compartment <- factor(df_long$Compartment,levels = c("S", "I", "R"),labels = c("Susceptible", "Infected", "Recovered"))
      
      threshold <- 1e-5
      valid_times <- df$time[df$I > threshold]
      
      if (length(valid_times)>0){
        t_min <- min(valid_times,na.rm=TRUE)
        t_max <- max(valid_times,na.rm=TRUE)
      } else{
        t_min <- 0
        t_max <- 50
      }
      
      p <- ggplot(df_long,aes(x = time, y = Value, color = Compartment,group=Compartment,text = paste0("Time: ", round(time, 2),"<br>", Compartment, ": ", round(Value,6),"<br>Effective Reproduction Number: ", round(Re, 3)))) +geom_line(linewidth = 0.6) + scale_color_manual(values=c("Susceptible"="steelblue","Infected"="firebrick","Recovered"="seagreen"))+ labs(x = "Time", y = "Proportion", color = "Compartment") + ggtitle(paste0("Epidemic SIR (R₀ = ",round(input$beta / input$gamma,3),")")) + coord_cartesian(xlim = c(t_min,t_max))
      ggplotly(p, tooltip = "text")
    })
    
  }
  

shinyApp(ui = ui, server = server)
