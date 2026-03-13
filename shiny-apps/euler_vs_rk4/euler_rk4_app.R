library(shiny)
library(bslib)
source("euler_rk4_solver.R")


ui <- page_sidebar(
  title = "Euler vs Fourth-Order Runge-Kutta Method Simulation",
  sidebar = sidebar(
    helpText(HTML("The curves solved using the Fourth-Order Runge-Kutta method are obtained using the function <code>rk4sys</code> in the <code>pracma</code> package. <br> <br> <b>Note:</b> The x-axis limits are dynamic.")),
    sliderInput("beta",label="Transmission Rate",min=0.1,max=2,value=1.5,step=0.1,),
    sliderInput("gamma",label="Recovery Rate",min=0.1,max=1,value=0.5,step=0.1),
    sliderInput("step",label = "Step Size",min=0.1,max=0.5,step=0.01,value=0.25)
  ),
  plotOutput("plot_euler"),
  plotOutput("plot_rk4"),
  plotOutput("plot_combine")
)


server <- function(input, output) {
  
  
  precomputed_euler <- reactive({
    euler_sir(S0=1-1e-6, I0 = 1e-6,0,beta=input$beta,gamma=input$gamma,len=500,dt=input$step)
  })
  
  precomputed_limits <- reactive({
    threshold <- 1e-5
    df <- precomputed_euler()
    valid_times <- df$times[df$I.vect > threshold]
    
    if (length(valid_times)>0){
      t_min <- min(valid_times,na.rm=TRUE)
      t_max <- max(valid_times,na.rm=TRUE)
    } else{
      t_min <- 0
      t_max <- 50
    }
    return(c(t_min,t_max))
  })
  
  
  output$plot_euler <- renderPlot({
    df_euler <- precomputed_euler() |> pivot_longer(cols=c(S.vect,I.vect,R.vect),names_to="n" ,values_to="v")
    
    
    limits <- precomputed_limits()
    
    ggplot(df_euler,aes(x=times,y=v,color=n)) + geom_line(linewidth=1) + labs(title="Epidemic SIR (Euler)",x="Time",y="Proportion",color="Compartment") + scale_color_manual(values = c("S.vect"="steelblue","I.vect"="firebrick","R.vect"="seagreen"),labels = c("S.vect" = "Susceptible","I.vect" = "Infected","R.vect" = "Recovered"),limits=c("S.vect","I.vect","R.vect"))+ coord_cartesian(xlim = limits)
 
     })
  
  output$plot_rk4 <- renderPlot({
    df_rk4 <- rk4_sir(S0=1-1e-6, I0 = 1e-6,0,b=input$beta,g=input$gamma,len=500,dt=input$step)
    limits <- precomputed_limits()
    
    ggplot(df_rk4,aes(x=time,y=v,color=n)) + geom_line(linewidth=1) + labs(title="Epidemic SIR (Fourth-Order Runge-Kutta)",x="Time",y="Proportion",color="Compartment") + scale_color_manual(values = c("y.S"="steelblue","y.I"="firebrick","y.R"="seagreen"),labels = c("y.S" = "Susceptible","y.I" = "Infected","y.R" = "Recovered"),limits=c("y.S","y.I","y.R"))+ coord_cartesian(xlim = limits)
    
  })
  
  output$plot_combine <- renderPlot({
    combined_df <- euler_rk4_combined(S0=1-1e-6, I0 = 1e-6,0,beta=input$beta,gamma=input$gamma,len=500,dt=input$step)
    limits <- precomputed_limits()
    
    ggplot(combined_df, aes(x = time, y = v, color = n, linetype = method, group = interaction(n, method))) + geom_line(linewidth = 1) + scale_linetype_manual(values = c("Euler" = "dotted", "Runge-Kutta" = "solid")) + scale_color_manual(values = c("S"="steelblue","I"="firebrick","R"="seagreen"), labels = c("S" = "Susceptible","I" = "Infected","R" = "Recovered"),limits=c("S","I","R")) + labs(title = "Epidemic SIR (Euler vs Fourth-Order Runge-Kutta)", x = "Time", y = "Proportion", color = "Compartment", linetype = "Solving Method") + coord_cartesian(xlim = limits)
    
  })
}

 
shinyApp(ui = ui, server = server)
