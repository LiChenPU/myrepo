library(shiny)
library(igraph)

## function ####
g <- graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = T)
## Analysis of Specific node ####
subgraph_specific_node = function(interested_node, g, step = 2)
{
  # interested_node = 1
  interested_node = interested_node
  
  g_id = g_vertex$ILP_id[g_vertex$ID==interested_node]
  
  # interested_node = as.character(g_id[1])
  g.degree <- degree(g, mode = c("all"))
  g_interest <- make_ego_graph(g, 
                              step, 
                              #1,
                              nodes = as.character(g_id[1]), 
                              mode = c("all"))[[1]]
  vertex.attributes(g_interest)$intensity[is.na(vertex.attributes(g_interest)$intensity)]=1e5
  
  # dists = distances(g_interest, g_id)
  colors <- c("black", "red", "orange", "blue", "dodgerblue", "cyan")
  # V(g_interest)$color <- colors[dists+1]
  # png(filename=paste("Subnetwork of node ", interested_node,".png",sep=""),
  #     width = 2400, height=2400,
  #     res=300)
  
  plot(g_interest,
       #vertex.color = 'white',
       vertex.label = vertex.attributes(g_interest)$formula,
       #vertex.label = vertex.attributes(g_interest)$medRt,
       vertex.label.color = "blue",
       vertex.label.cex = 1,
       #vertex.label.dist = 2,
       vertex.size = log(vertex.attributes(g_interest)$intensity)-3,
       #edge.width = edge.attributes(g_interest)$Confidence*2-2,
       # edge.color = color_palette[edge.attributes(g_interest)$color],
       #edge.label = edge.attributes(g_interest)$mz_dif,
       edge.label = edge.attributes(g_interest)$linktype,
       edge.label.color = "red",
       edge.label.cex = 1,
       edge.arrow.size = 0.5,
       edge.arrow.width = 1,
       main = paste("Subnetwork of node", interested_node),
       layout = layout_nicely(g_interest)
       #layout = layout_as_tree(g_interest)
  )
  # dev.off()
}



## -------- Shiny R --------------------####
ui <- fluidPage(
  wellPanel(
  numericInput(inputId = "Peak_id", 
              label = "Enter a peak ID",
              value = 122
  ),
  actionButton(inputId = "id_update", 
               label = "Go")
  ),
  mainPanel(
    plotOutput("graph")
  )
  
)

server <- function(input, output) {
  observeEvent(input$id_update, {
    output$graph <- renderPlot({
      subgraph_specific_node(input$Peak_id, g,1)
    })
  })
}

shinyApp(ui = ui, server = server)