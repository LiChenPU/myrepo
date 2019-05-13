library(shiny)
library(igraph)

## function ####
g <- graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = T)
## Analysis of Specific node ####
g_interest_node = function(interested_node, g, step = 2)
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
  return(g_interest)
}

Plot_g_interest = function(g_interest, interested_node)
{
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


g_show_vertice = function(g_interest)
{
  g_interest_vertice = igraph::as_data_frame(g_interest, "vertices")
  colname_vertice = c("ID", "mz", "RT", "formula", "compound_name", "is_artifact", "is_biotransform")
  g_interest_vertice = g_interest_vertice[,colname_vertice]
  return(g_interest_vertice)
}


g_show_edge = function(g_interest)
{
  g_interest_edge = igraph::as_data_frame(g_interest, "edges")
  colname_edge = c("node1", "node2", "formula1", "formula2", "linktype")
  g_interest_edge = g_interest_edge[,colname_edge]
  return(g_interest_edge)
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
    plotOutput("graph"
               , width = "100%"
               ),
    dataTableOutput("nodetable"),
    dataTableOutput("edgetable")
    
  )
)

server <- function(input, output) {
  
  g_interest <- reactive({
    g_interest_node(input$Peak_id, g,1)
    
  })

  
  output$graph <- renderPlot({
    # g_interest = g_interest_node(input$Peak_id, g,1)
    Plot_g_interest(g_interest(), input$Peak_id)
  }
  # ,height = 400, width = 800
  )
  output$nodetable <- renderDataTable({
    g_interest_vertice = g_show_vertice(g_interest())
  })
  output$edgetable <- renderDataTable({
    g_interest_edge = g_show_edge(g_interest())
  })
}

shinyApp(ui = ui, server = server)