
###########################################################
# Aplicació shinny per jugar amb una binomial
# Especialment pensada per l'exercici d'estimar el %
# de carreteres pavimentades a Uganda: https://www.mapcrunch.com/
# Data 8/7/2023
##########################################################


library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Distribució Binomial:"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("num_intents", "Número d'intents:", value = 10),
      numericInput("num_encerts", "Número d'encerts:", value = 5),
      numericInput("prob_esperada", "Probabilitat esperada:", value = 0.5),
      actionButton("calcular", "Calcular")
    ),
    
    mainPanel(
      plotOutput("plot"),
      textOutput("prob_output"),
      textOutput("prob_output2"),
      textOutput("prob_output3")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$calcular, {
    num_intents <- input$num_intents
    num_encerts <- input$num_encerts
    prob_esperada <- input$prob_esperada
    
    # Càlcul de l'interval de confiança del 95% per a la probabilitat d'encert
    ci <- prop.test(num_encerts, num_intents, correct = FALSE)
    prob_estimada <- ci$estimate
    lower_bound <- ci$conf.int[1]
    upper_bound <- ci$conf.int[2]
    
    # Càlcul de la distribució binomial
    x <- 0:num_intents
    y <- dbinom(x, size = num_intents, prob = prob_esperada)
    densidad <- dnorm(x, mean = num_intents * prob_esperada, sd = sqrt(num_intents * prob_esperada * (1 - prob_esperada)))
    data <- data.frame(x, y, densidad)
    
    prob_observada<-(dbinom(num_encerts, size = num_intents, prob = prob_esperada))
    z <- num_encerts:num_intents
    prob_observada2<-sum(dbinom(z, size = num_intents, prob = prob_esperada))
    
    # Gràfic de la distribució binomial
    output$plot <- renderPlot({
      ggplot(data, aes(x)) +
        geom_bar(aes(y = y), stat = "identity", fill = "green", alpha = 0.5) +
        geom_line(aes(y = densidad), color = "darkgreen", size = 1, linetype = "dashed") +
        geom_vline(xintercept = num_encerts, color = "red", linetype = "solid", size = 2) +
        labs(x = "Número d'encerts", y = "Probabilitat / Densitat", title = "Distribució Binomial amb n = Número d'intents i  p = probabilitat esperada.") +
        theme_minimal()
    })
    
    # Output de la probabilitat estimada amb l'interval de confiança
    output$prob_output <- renderText({
      paste0("Probabilitat estimada d'encerts: ", round(prob_estimada * 100, 2), "%",
             " (IC 95%: [", round(lower_bound * 100, 2), "%, ", round(upper_bound * 100, 2), "%])")
    })
    # Output de la probabilitat observada 
    output$prob_output2 <- renderText({
      paste0("Probabilitat d'observar aquest valor: ", round(prob_observada * 100, 2))
    })
    # Output de la probabilitat observar aquest valor o major
    output$prob_output3 <- renderText({
      paste0("Probabilitat d'observar aquest valor o major: ", round(prob_observada2 * 100, 2))
    })
  })
  
}

shinyApp(ui, server)