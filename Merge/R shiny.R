library(shiny)
library(igraph)

## Read in files

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("./xi_new_neg")

g_vertex = read.csv("g_vertex.txt", stringsAsFactors = F)
g_edge = read.csv("g_edge.txt", stringsAsFactors = F)

g <- graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = T)
## function ####






{
all_two_nodes_graph = function(g, node1, node2, dist = 3)
{
  # node1 = 122
  # node2 = 389
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_id1 = g_vertex$name[g_vertex$ID==node1& g_vertex$ILP_result!=0]
  g_id2 = g_vertex$name[g_vertex$ID==node2& g_vertex$ILP_result!=0]
  g_list = list()
  i = g_id1[1]
  j = g_id2[1]
  for(i in g_id1){
    for(j in g_id2){
      if(distances(g, i, j, mode = "out") > dist){next}
      paths_connect_ij = all_simple_paths(g,i,j, mode = "out")
      ids.to.keep = sapply(paths_connect_ij, function(i) length(i)<=(dist+1))
      paths_connect_ij = paths_connect_ij[ids.to.keep]
      all_names_in_connect_graph = unique(names(unlist(paths_connect_ij)))
      if(length(all_names_in_connect_graph)==0){next}
      g_temp_vertex = g_vertex[g_vertex$name %in% all_names_in_connect_graph,]
      g_temp_edge = g_edge[g_edge$from %in% all_names_in_connect_graph & g_edge$to %in% all_names_in_connect_graph,]
      g_temp = graph_from_data_frame(d = g_temp_edge, vertices = g_temp_vertex, directed = T)
      g_list[[length(g_list)+1]] = g_temp
    }
  }
  return(g_list)
}

shortest_two_nodes_graph = function(g, node1, node2)
{
  # node1 = 122
  # node2 = 389
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_id1 = g_vertex$name[g_vertex$ID==node1& g_vertex$ILP_result!=0]
  g_id2 = g_vertex$name[g_vertex$ID==node2& g_vertex$ILP_result!=0]
  g_list = list()
  i = g_id1[1]
  j = g_id2[1]
  for(i in g_id1){
    for(j in g_id2){
      if(are.connected(g,i,j)){next}
      paths_connect_ij = shortest_paths(g,i,j, mode = "all")
      # ids.to.keep = sapply(paths_connect_ij, function(i) length(i)<=dist)
      # paths_connect_ij = paths_connect_ij[ids.to.keep]
      all_names_in_connect_graph = unique(names(unlist(paths_connect_ij$vpath)))
      if(length(all_names_in_connect_graph)==0){next}
      g_temp_vertex = g_vertex[g_vertex$name %in% all_names_in_connect_graph,]
      g_temp_edge = g_edge[g_edge$from %in% all_names_in_connect_graph & g_edge$to %in% all_names_in_connect_graph,]
      g_temp = graph_from_data_frame(d = g_temp_edge, vertices = g_temp_vertex, directed = T)
      g_list[[length(g_list)+1]] = g_temp
    }
  }
  return(g_list)
}



test = all_two_nodes_graph(g, 122,389,3)
test = shortest_two_nodes_graph(g, 121,389)
Plot_g_interest(test[[1]],1)

search_peak = function(g, mz_interest, mz_ppm, mz_lb, mz_ub, rt_lb, rt_ub)
{
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_vertex = g_vertex[g_vertex$ILP_result!=0,]
  
  g_vertex = g_vertex[g_vertex$mz<mz_interest*(1+mz_ppm/10^6) & 
                        g_vertex$mz>mz_interest*(1-mz_ppm/10^6),]
  g_vertex = g_vertex[g_vertex$mz>mz_lb & 
                        g_vertex$mz<mz_ub,]
  g_vertex = g_vertex[g_vertex$RT>rt_lb & 
                        g_vertex$RT<rt_ub,]
  
  colname_vertice = c("ID", "mz", "RT", "formula", "compound_name", "is_artifact", "is_biotransform", "is_metabolite")
  g_vertex = g_vertex[,colname_vertice]
  return(g_vertex)
}

filter_graph = function(g, intensity_low_cutoff=3, intensity_high_cutoff=10)
{
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_edge = igraph::as_data_frame(g, "edges")
  
  g_vertex = g_vertex[is.na(g_vertex$intensity) |
                        g_vertex$intensity > 10^intensity_low_cutoff &
                        g_vertex$intensity < 10^intensity_high_cutoff,]
  g_edge = g_edge[g_edge$from %in% g_vertex$name &
                    g_edge$to %in% g_vertex$name, ]
  
  g_filter = graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = T)
  
  return(g_filter)
}
interest_node_graph = function(interested_node, g, 
                           step = 1
                           )
{
  # interested_node = 122
  interested_node = interested_node

  g_vertex = igraph::as_data_frame(g, "vertices")
  g_id = g_vertex$name[g_vertex$ID==interested_node]
  
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
  colname_vertice = c("ID", "mz", "RT", "formula", "compound_name", "is_artifact", "is_biotransform", "is_metabolite")
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

}

## Shiny R --------------------####
ui <- fluidPage(
  sidebarPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Setting",
                  sliderInput(inputId = "Peak_inten_range", 
                              label = "Peak Intensity (log10)",
                              min = 0, max = 10, step = 0.01,
                              value = c(3,10)),
                  numericInput(inputId = "mz_interest", 
                               label = "Enter a mz of interest",
                               value = 180.0633
                  ),
                  numericInput(inputId = "mz_ppm", 
                               label = "ppm",
                               value = 5
                  ),
                  numericInput(inputId = "Peak_id", 
                               label = "Enter a peak ID",
                               value = 122
                  )
                  # actionButton(inputId = "id_update", 
                  #              label = "Go")
                  
                ),
                tabPanel("Advanced",
                   sliderInput(inputId = "mz_range", 
                               label = "m/z range",
                               min = 0, max = 1500, step = 1,
                               value = c(0,1500)),
                   sliderInput(inputId = "rt_range", 
                               label = "RT range",
                               min = 0, max = 20, step = .01,
                               value = c(0,20))
                )
                

    )
  ),
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Peak_table", dataTableOutput("peak_table")),
                tabPanel("Plot",  plotOutput("graph"
                                             , width = "100%"
                )),
                tabPanel("Nodes", dataTableOutput("nodetable")),
                tabPanel("Edges", dataTableOutput("edgetable"))
    )
    
    
  )
)



server <- function(input, output) {
  # g_interest <- eventReactive(
  #   input$id_update,
  
  output$peak_table <- renderDataTable({
    search_peak(g, mz_interest = input$mz_interest, mz_ppm = input$mz_ppm,
                mz_lb = input$mz_range[1], mz_ub = input$mz_range[2], 
                rt_lb = input$rt_range[1], rt_ub = input$rt_range[2])
  
  })
  
  g_filter <- reactive(
    {
    filter_graph(g,
                input$Peak_inten_range[1], input$Peak_inten_range[2])
  })
  g_interest <- reactive({
    interest_node_graph(input$Peak_id, 
                         g_filter(),1
                         )
  })
  
  output$graph <- renderPlot({
    Plot_g_interest(g_interest(), isolate(input$Peak_id))
    # tkplot(g_interest())
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



# runExample("01_hello")      # a histogram
# runExample("02_text")       # tables and data frames
# runExample("03_reactivity") # a reactive expression
# runExample("04_mpg")        # global variables
# runExample("05_sliders")    # slider bars
# runExample("06_tabsets")    # tabbed panels
# runExample("07_widgets")    # help text and submit buttons
# runExample("08_html")       # Shiny app built from HTML
# runExample("09_upload")     # file upload wizard
# runExample("10_download")   # file download wizard
# runExample("11_timer")      # an automated timer
