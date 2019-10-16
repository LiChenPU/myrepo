require(shiny)
require(visNetwork)

server <- function(input, output) {
  output$network_proxy_nodes <- renderVisNetwork({
    # minimal example
    nodes <- data.frame(id = 1:3)
    edges <- data.frame(from = c(1,2), to = c(1,3))
    
    visNetwork(nodes, edges) %>%
      visOptions(manipulation = TRUE, 
                 highlightNearest = TRUE, 
                 nodesIdSelection = TRUE
                 # selectedBy = "is_metabolite"
      ) %>%
      visInteraction(navigationButtons = TRUE)
  })

  observe({
    visNetworkProxy("network_proxy_nodes") %>%
      visFocus(id = input$Focus, scale = 2)
  })
  
  observe({
    visNetworkProxy("network_proxy_nodes") %>%
      visNodes(color = input$color)
  })
}

ui <- fluidPage(
  fluidRow(
    column(
      width = 4,
      selectInput("color", "Color :",
                  c("blue", "red", "green")),
      selectInput("Focus", "Focus on node :",
                  c(1:3))
    ),
    column(
      width = 8,
      visNetworkOutput("network_proxy_nodes", height = "400px")
    )
  )
)

shinyApp(ui = ui, server = server)
