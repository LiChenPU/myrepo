library(shiny)
library(igraph)
library(reactlog)
library(ShinyTester)



# options(shiny.reactlog=TRUE) 
# Ctrl + F3 to view 

# Read in files ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
datapath = ("./Melanie_merge/merge3")
setwd(datapath)


g_vertex = read.csv("g_vertex.txt", stringsAsFactors = F)
g_edge = read.csv("g_edge.txt", stringsAsFactors = F)

g <- graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = F)
g_vertex = igraph::as_data_frame(g, "vertices")
g_edge = igraph::as_data_frame(g, "edges")

# g_vertex_ILP = g_vertex[g_vertex$ILP_result != 0 & !is.na(g_vertex$ILP_result),]
# write.csv(g_vertex_ILP, "g_vertex_ILP.txt", row.names = F)

# function ####
## filter_graph ####
filter_graph = function(g, 
                        intensity_lb = 3, intensity_ub= 7, 
                        mz_lb = 0, mz_ub = 1500, 
                        rt_lb = 0, rt_ub = 30,
                        is_metabolite_group)
{
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_edge = igraph::as_data_frame(g, "edges")
  
  g_vertex = g_vertex[g_vertex$ILP_result!=0 & (!is.na(g_vertex$ILP_result)),]
  g_vertex = g_vertex[is.na(g_vertex$intensity) |
                        g_vertex$intensity > 10^intensity_lb &
                        g_vertex$intensity < 10^intensity_ub,]
  g_vertex = g_vertex[g_vertex$mz>mz_lb & 
                        g_vertex$mz<mz_ub,]
  g_vertex = g_vertex[g_vertex$RT>rt_lb & 
                        g_vertex$RT<rt_ub,]
  
  is_metabolite_group = unique(g_vertex$is_metabolite)
  is_metabolite_group[is_metabolite_group==""] = NA
  g_vertex = g_vertex[g_vertex$is_metabolite %in% is_metabolite_group,]
  
  g_edge = g_edge[g_edge$from %in% g_vertex$name &
                    g_edge$to %in% g_vertex$name, ]
  
  g_filter = graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = F)
  
  return(g_filter)
}
## search_peak ####
search_peak = function(g, mz_interest, mz_ppm)
{
  
  # mz_interest = 558.51
  # mz_ppm = 5
  # mz_lb = 0
  # mz_ub = 1500
  # rt_lb = 0
  # rt_ub = 20
  
  g_vertex = igraph::as_data_frame(g, "vertices")
  # g_vertex = g_vertex[g_vertex$ILP_result!=0 & (!is.na(g_vertex$ILP_result)),]
  
  # filter
  g_vertex = g_vertex[g_vertex$mz<mz_interest*(1+mz_ppm/10^6) & 
                        g_vertex$mz>mz_interest*(1-mz_ppm/10^6),]
  
  # ranking
  g_vertex = g_vertex[with(g_vertex, order(abs(mz-mz_interest)), -ILP_result, abs(cal_mass - mz_interest)),]
  
  
  colname_vertice = c("ID", "mz", "RT", "formula", "compound_name", "is_artifact", "is_biotransform", "is_metabolite")
  g_vertex = g_vertex[,colname_vertice]
  
  return(g_vertex)
}

## search_formula ####
search_formula = function(g, peak_id){
  if(is.na(peak_id)){return(NULL)}
  g_vertex = igraph::as_data_frame(g, "vertices")
  
  # filter
  # g_vertex = g_vertex[g_vertex$ILP_result!=0 & (!is.na(g_vertex$ILP_result)),]
  g_vertex = g_vertex[g_vertex$ID==peak_id,]

  
  # ranking
  g_vertex = g_vertex[with(g_vertex, order(-ILP_result, -cplex_score)),]
  
  colname_vertice = c("ID", "mz", "RT", "formula", "compound_name", "intensity", "is_metabolite")
  g_vertex = g_vertex[,colname_vertice]
  
  return(g_vertex)
  
}

## search_partner ####
search_partner = function(g, peak_id, formula_select, step = 5 ){

  if(!is_igraph(g)) {print("return !is_igraph")
    return(NULL)}
  
  interested_node = peak_id
  
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_id = g_vertex$name[g_vertex$ID==interested_node & g_vertex$formula == formula_select]
  # print(g_id)
  if(length(g_id) == 0){
    # print("return g_id ==0")
    return(NULL)}
    
  
  # interested_node = as.character(g_id[1])
  # g.degree <- degree(g, mode = c("all"))
  g_partner <- make_ego_graph(g, 
                               step, 
                               nodes = as.character(g_id[1]), 
                               mode = c("all"))[[1]]
  
  return(g_partner)
}


## g_show_vertice_rankIntensity ####
g_show_vertice_rankIntensity = function(g_interest)
{
  
  if(!is_igraph(g_interest)) {return(NULL)}
  g_interest_vertice = igraph::as_data_frame(g_interest, "vertices")
  g_interest_vertice = g_interest_vertice[with(g_interest_vertice, order(-intensity)),]
    
  colname_vertice = c("ID", "mz", "RT", "formula", "compound_name", "intensity", "is_metabolite")
  g_interest_vertice = g_interest_vertice[,colname_vertice]
  return(g_interest_vertice)
}




## interest_node_graph ####
interest_node_graph = function(g, peak_id, formula_select, step = 1)
{
  if(!is_igraph(g)) {return(NULL)}
  # interested_node = 122
  interested_node = peak_id

  g_vertex = igraph::as_data_frame(g, "vertices")
  g_id = g_vertex$name[g_vertex$ID==interested_node & g_vertex$formula == formula_select]
  
  # interested_node = as.character(g_id[1])
  # g.degree <- degree(g, mode = c("all"))
  g_interest <- make_ego_graph(g, 
                               step, 
                               nodes = as.character(g_id[1]), 
                               mode = c("all"))[[1]]
  return(g_interest)
}

## Plot_g_interest ####
Plot_g_interest = function(g_interest, interested_node, formula_select)
{
  
  # interested_node=178
  # formula_select="C6H12O6"
  # g_interest = g_temp
  
  if(!is_igraph(g_interest)) {return(NULL)}
  vertex.attributes(g_interest)$intensity[is.na(vertex.attributes(g_interest)$intensity)]=1e5
  
  # dists = distances(g_interest, g_id)
  colors <- c("black", "red", "orange", "blue", "dodgerblue", "cyan")
  
  plot(g_interest,
       #vertex.color = 'white',
       vertex.label = vertex.attributes(g_interest)$formula,
       #vertex.label = vertex.attributes(g_interest)$medRt,
       vertex.color = "orange",
       vertex.label.color = "black",
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
       main = paste("Subnetwork of peak", interested_node, formula_select),
       layout = layout_nicely(g_interest)
       #layout = layout_as_tree(g_interest)
  )
  # dev.off()
}


## g_show_vertice & g_show_edge####
g_show_vertice = function(g_interest)
{
  if(!is_igraph(g_interest)) {return(NULL)}
  g_interest_vertice = igraph::as_data_frame(g_interest, "vertices")
  colname_vertice = c("ID", "mz", "RT", "formula", "compound_name","intensity", "is_metabolite")
  g_interest_vertice = g_interest_vertice[,colname_vertice]
  return(g_interest_vertice)
}

g_show_edge = function(g_interest)
{
  if(!is_igraph(g_interest)) {return(NULL)}
  g_interest_edge = igraph::as_data_frame(g_interest, "edges")
  colname_edge = c("node1", "node2", "formula1", "formula2", "linktype")
  g_interest_edge = g_interest_edge[,colname_edge]
  return(g_interest_edge)
}


## two_formula_neighbor_graph####

two_formula_neighbor_graph = function(g, node1, node2, formula1, formula2, dist = 3)
{
  # node1 = 178
  # node2 = 2
  # formula1 = "C6H12O6"
  # formula2 = "C16H32O2"
  # dist = 3
  if(!is_igraph(g)) {return(NULL)}
  
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_id1 = g_vertex$name[g_vertex$ID==node1 & g_vertex$formula==formula1]
  g_id2 = g_vertex$name[g_vertex$ID==node2 & g_vertex$formula==formula2]
  
  i = g_id1[1]
  j = g_id2[1]
  
  if(distances(g, i, j, mode = "all") > dist){return(NULL)}
  
  all_names_in_connect_graph = c()
  for(count in 0:dist){
    g_i = make_ego_graph(g, 
                         count, 
                         nodes = i, 
                         mode = c("all"))[[1]]
    g_j = make_ego_graph(g, 
                         dist-count, 
                         nodes = j, 
                         mode = c("all"))[[1]]
    
    g_i_vertex =  igraph::as_data_frame(g_i, "vertices")
    g_j_vertex =  igraph::as_data_frame(g_j, "vertices")
    
    intersect_names = intersect(g_i_vertex$name, g_j_vertex$name)
    all_names_in_connect_graph = c(all_names_in_connect_graph, intersect_names)
  }

  all_names_in_connect_graph = unique(all_names_in_connect_graph)

  if(length(all_names_in_connect_graph)==0){return(NULL)}
  g_temp_vertex = g_vertex[g_vertex$name %in% all_names_in_connect_graph,]
  g_temp_edge = g_edge[g_edge$from %in% all_names_in_connect_graph & g_edge$to %in% all_names_in_connect_graph,]
  g_temp = graph_from_data_frame(d = g_temp_edge, vertices = g_temp_vertex, directed = T)
  
  return(g_temp)
}

## shortest_two_nodes_shortest####
two_formula_shortest_path_graph = function(g, node1, node2, formula1, formula2)
{
  # node1 = 178
  # node2 = 2
  # formula1 = "C6H12O6"
  # formula2 = "C16H32O2"

  if(!is_igraph(g)) {return(NULL)}
  
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_id1 = g_vertex$name[g_vertex$ID==node1 & g_vertex$formula==formula1]
  g_id2 = g_vertex$name[g_vertex$ID==node2 & g_vertex$formula==formula2]
  
  i = g_id1[1]
  j = g_id2[1]
  
  paths_connect_ij = shortest_paths(g,i,j, mode = "all")
  # ids.to.keep = sapply(paths_connect_ij, function(i) length(i)<=dist)
  # paths_connect_ij = paths_connect_ij[ids.to.keep]
  all_names_in_connect_graph = unique(names(unlist(paths_connect_ij$vpath)))
  
  if(length(all_names_in_connect_graph)==0){return(NULL)}
  g_temp_vertex = g_vertex[g_vertex$name %in% all_names_in_connect_graph,]
  g_temp_edge = g_edge[g_edge$from %in% all_names_in_connect_graph & g_edge$to %in% all_names_in_connect_graph,]
  g_temp = graph_from_data_frame(d = g_temp_edge, vertices = g_temp_vertex, directed = T)

  return(g_temp)
}




# Shiny R --------------------####
## ui ####
ui <- fluidPage(
  sidebarPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Setting",
                  numericInput(inputId = "mz_interest", 
                               label = "Enter a mz of interest",
                               value = 180.0633
                  ),
                  selectInput(inputId = "ion_form", 
                               label = "Select ionization",
                              choices = c("M+H","M-H","M"),
                              selected = "M"
                  ),
                  numericInput(inputId = "mz_ppm", 
                               label = "ppm",
                               value = 5
                  ),
                  numericInput(inputId = "Peak_id", 
                               label = "Peak ID",
                               value = numeric(0)
                  ),
                  selectInput(inputId = "formula_select",
                              label = "Formula",
                              choices = character(0)
                              ),
                  
                  wellPanel(
                    numericInput(inputId = "Partner_level", 
                                 label = "Level of partner",
                                 value = 2
                    ),
                    numericInput(inputId = "Partner_id", 
                                 label = "Partner ID",
                                 value = numeric(0)
                    ),
                    selectInput(inputId = "Partner_formula",
                                label = "Partner formula",
                                choices = character(0)
                    )
                  )


                  
                ),
                tabPanel("Advanced",
                   sliderInput(inputId = "Peak_inten_range", 
                               label = "Peak Intensity (log10)",
                               min = 0, max = 10, step = 0.01,
                               value = c(3,10)),
                   sliderInput(inputId = "mz_range", 
                               label = "m/z range",
                               min = 0, max = 1500, step = 1,
                               value = c(0,1500)),
                   sliderInput(inputId = "rt_range", 
                               label = "RT range",
                               min = 0, max = 20, step = .01,
                               value = c(0,20)),
                   checkboxGroupInput("is_metabolite", "Select whether to include:",
                                      choices = unique(g_vertex$is_metabolite), 
                                      selected = unique(g_vertex$is_metabolite)
                   )
                )
                
    ),
    actionButton("one_node_graph", "one_node_graph"),
    actionButton("two_nodes_shortest", "two_nodes_shortest"),
    actionButton("two_nodes_all_graph", "two_nodes_all_graph")
    
  ),
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Peak_table", dataTableOutput("peak_table")),
                tabPanel("Partner_table", dataTableOutput("Partner_table")),
                tabPanel("Plot",  plotOutput("graph"
                                             , width = "100%"
                )),
                tabPanel("Nodes", dataTableOutput("nodetable")),
                tabPanel("Edges", dataTableOutput("edgetable"))
    )
  )
)


## server ####
server <- function(input, output, session) {
  ## Filter graph based on intensity, mz range and rt range
  g_filter <- reactive({
    print("enter g_filter")
    filter_graph(g,
                 input$Peak_inten_range[1], input$Peak_inten_range[2],
                 mz_lb = input$mz_range[1], mz_ub = input$mz_range[2], 
                 rt_lb = input$rt_range[1], rt_ub = input$rt_range[2],
                 input$is_metabolite)
  })
  
  ## adjust mz based on ionization and ppm
  mz_interest <- reactive({
    
    print("enter mz_interest")
    ion_form = input$ion_form
    if(ion_form == "M"){mz_adjust = input$mz_interest}
    if(ion_form == "M+H"){mz_adjust = input$mz_interest - 1.007276}
    if(ion_form == "M-H"){mz_adjust = input$mz_interest + 1.007276}
    mz_adjust
  })
  
  ## show peak table that meet the mz requirement
  peak_table <- reactive({
    print("enter peak_table")
    search_peak(g_filter(), mz_interest = mz_interest(), mz_ppm = input$mz_ppm)
  })
  
  output$peak_table <- renderDataTable({
    print("output peak_table")
    peak_table()
  })
  
  ## update Peak_id
  observe({
    print("enter update peak_id")
    x<-peak_table()
    if (nrow(x)==0){
      value <- numeric(0)
    } else{
      value <- x$ID[1]
    }
    updateNumericInput(session, "Peak_id",
                      value = value
    )
  })
  
  ## update formula_select
  observe({
    print("enter update formula_select")
    x = search_formula(isolate(g_filter()), input$Peak_id)
    if (is.null(x)){
      x <- " "
    } else{
      x <- x$formula
    }
    updateSelectInput(session, "formula_select",
                      # label = paste("Formula", length(x)),
                      choices = x,
                      selected = head(x, 1)
    )
  })
  
  
  ## Partner graph
  g_partner <- reactive({
    print("enter g_partner")
    # g_partner <-   eventReactive(input$mz_interest, {
    # req(input$formula_select, cancelOutput = F)
    search_partner(isolate(g_filter()), isolate(input$Peak_id), input$formula_select, 5)
  })
  
  ## Partner table
  Partner_table <- reactive({
    print("enter Partner_table")
    g_show_vertice_rankIntensity(g_partner())
  })
  
  
  ## Output partner table
  output$Partner_table <- renderDataTable({
    print("output Partner_table")
    Partner_table()
    })
  
  observe({
    print("enter update Partner_id")
    # observeEvent(Partner_table(), {
    x = Partner_table()
    if (is.null(x)){
      x <- numeric(0)
    } else{
      x <- x$ID[1]
    }
    
    updateNumericInput(session, "Partner_id",
                       value = x
    )
  })
  
  observe({
    print("enter update Partner_formula")
    x = search_formula(isolate(g_filter()), input$Partner_id)
    if (is.null(x)){
      x <- character(0)
    } else{
      x <- x$formula
    }
    updateSelectInput(session, "Partner_formula",
                      # label = paste("Formula", length(x)),
                      choices = x,
                      selected = head(x, 1)
    )
  })
  
  
  ## select graph of interest 
  
  g_interest <- reactiveValues(data = NULL)
  observeEvent(input$one_node_graph, {
    print("enter g_interest interest_node_graph")
    ### basic graph based on selected peak and formula
    g_interest$data = interest_node_graph(g_partner(), 
                                          isolate(input$Peak_id), 
                                          isolate(input$formula_select), 
                                          step = 1)
  })

  observeEvent(input$two_nodes_shortest, {
    print("enter g_interest two_formula_shortest_path_graph")
    g_interest$data = two_formula_shortest_path_graph(isolate(g_partner()), 
                                                      isolate(input$Peak_id), 
                                                      isolate(input$Partner_id), 
                                                      isolate(input$formula_select), 
                                                      input$Partner_formula)
  })
  
  observeEvent(input$two_nodes_all_graph, {
    print("enter g_interest two_formula_neighbor_graph")
    g_interest$data = two_formula_neighbor_graph(isolate(g_partner()), 
                                                 isolate(input$Peak_id), 
                                                 isolate(input$Partner_id), 
                                                 isolate(input$formula_select), 
                                                 input$Partner_formula, 
                                                 dist = input$Partner_level)
  })
  
  
  # g_interest <- reactive({
  #   print("enter g_interest two_formula_shortest_path_graph")
  #   two_formula_shortest_path_graph(isolate(g_partner()), isolate(input$Peak_id), isolate(input$Partner_id), isolate(input$formula_select), input$Partner_formula)
  # })


  
  ## output graph, related node & edge table
  output$graph <- renderPlot({
    print("enter Plot_g_interest")
    Plot_g_interest(g_interest$data, isolate(input$Peak_id), isolate(input$formula_select))
    # tkplot(g_interest())
  }
  # ,height = 400, width = 800
  )
  
  output$nodetable <- renderDataTable({
    print("enter show g_interest_vertice")
    g_interest_vertice = g_show_vertice(g_interest$data)
  })
  output$edgetable <- renderDataTable({
    print("enter show g_interest_edge")
    g_interest_edge = g_show_edge(g_interest$data)
  })
  
  
}

## Run shiny ####
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

# Debug ####
# {
#   output = list()
#   input = list()
# 
#   input$mz_interest = 180.0633
#   input$ion_form = "M"
#     if(input$ion_form == "M"){mz_interest = input$mz_interest}
#     if(input$ion_form == "M+H"){mz_interest = input$mz_interest - 1.007276}
#     if(input$ion_form == "M-H"){mz_interest = input$mz_interest + 1.007276}
# 
# 
# 
#   input$mz_ppm = 5
# 
# 
# 
#   peak_table <- search_peak(g, mz_interest = input$mz_interest, mz_ppm = input$mz_ppm)
#   peak_table <- search_peak(g, mz_interest = mz_interest, mz_ppm = input$mz_ppm)
#   output$peak_table <- peak_table
# 
#   input$Peak_id = peak_table$ID[1]
# 
#   x = search_formula(g, input$Peak_id)
#   if (is.null(x)){
#     x <- character(0)
#   } else{
#     x <- x$formula
#   }
# 
#   input$formula_select = list(choices = x,
#                               selected = head(x, 1))
# 
#   search_partner = function(g, peak_id, formula_select, step = 3 ){
#     # interested_node = 122
#     formula_select = input$formula_select$selected
#     interested_node = input$Peak_id
# 
#     g_vertex = igraph::as_data_frame(g, "vertices")
#     g_id = g_vertex$name[g_vertex$ID==interested_node & g_vertex$formula == formula_select]
# 
#     # g.degree <- degree(g, mode = c("all"))
#     g_interest <- make_ego_graph(g_filter,
#                                  order = step,
#                                  # nodes = as.character(g_id[1]),
#                                  nodes = g_id[1],
#                                  mode = c("in"))[[1]]
# 
#     g_vertex_test = igraph::as_data_frame(g_interest, "vertices")
#     return(g_interest)
#   }
# 
# 
# 
# 
#   input$Peak_inten_range = c(3,10)
#   input$mz_range = c(0,1500)
#   input$rt_range = c(0,20)
# 
#   g_filter <- filter_graph(g,
#                  input$Peak_inten_range[1], input$Peak_inten_range[2],
#                  mz_lb = input$mz_range[1], mz_ub = input$mz_range[2],
#                  rt_lb = input$rt_range[1], rt_ub = input$rt_range[2])
# 
#   g_filter_vertex = igraph::as_data_frame(g_filter, "vertice")
# 
#   g_interest <- interest_node_graph(g_filter,input$Peak_id,input$formula_select$selected, step = 1)
# 
# 
# 
#   output$graph <- Plot_g_interest(g_interest, isolate(input$Peak_id), isolate(input$formula_select$selected))
#   output$nodetable <- g_show_vertice(g_interest)
#   output$edgetable <- g_show_edge(g_interest)
# 
#   output$partner_table = 0
# 
#   test = all_two_nodes_shortest(g, 178,2,3)
#   test = shortest_two_nodes_shortest(g, 178,2)
#   Plot_g_interest(test[[1]],input$Peak_id,input$formula_select$selected)
# }
# 
# 
# Unused function ####

## all_two_nodes_shortest 
# 
# all_two_nodes_shortest = function(g, node1, node2, dist = 3)
# {
#   # node1 = 122
#   # node2 = 389
#   g_vertex = igraph::as_data_frame(g, "vertices")
#   g_id1 = g_vertex$name[g_vertex$ID==node1& g_vertex$ILP_result!=0]
#   g_id2 = g_vertex$name[g_vertex$ID==node2& g_vertex$ILP_result!=0]
#   g_list = list()
#   i = g_id1[1]
#   j = g_id2[1]
#   for(i in g_id1){
#     for(j in g_id2){
#       if(distances(g, i, j, mode = "out") > dist){next}
#       paths_connect_ij = all_simple_paths(g,i,j, mode = "out")
#       ids.to.keep = sapply(paths_connect_ij, function(i) length(i)<=(dist+1))
#       paths_connect_ij = paths_connect_ij[ids.to.keep]
#       all_names_in_connect_graph = unique(names(unlist(paths_connect_ij)))
#       if(length(all_names_in_connect_graph)==0){next}
#       g_temp_vertex = g_vertex[g_vertex$name %in% all_names_in_connect_graph,]
#       g_temp_edge = g_edge[g_edge$from %in% all_names_in_connect_graph & g_edge$to %in% all_names_in_connect_graph,]
#       g_temp = graph_from_data_frame(d = g_temp_edge, vertices = g_temp_vertex, directed = T)
#       g_list[[length(g_list)+1]] = g_temp
#     }
#   }
#   return(g_list)
# }
# 
# ## shortest_two_nodes_shortest
# shortest_two_nodes_shortest = function(g, node1, node2)
# {
#   # node1 = 178
#   # node2 = 2
#   g_vertex = igraph::as_data_frame(g, "vertices")
#   g_id1 = g_vertex$name[g_vertex$ID==node1& g_vertex$ILP_result!=0]
#   g_id2 = g_vertex$name[g_vertex$ID==node2& g_vertex$ILP_result!=0]
#   g_list = list()
#   i = g_id1[1]
#   j = g_id2[1]
#   for(i in g_id1){
#     for(j in g_id2){
#       if(are.connected(g,i,j)){next}
#       paths_connect_ij = shortest_paths(g,i,j, mode = "all")
#       # ids.to.keep = sapply(paths_connect_ij, function(i) length(i)<=dist)
#       # paths_connect_ij = paths_connect_ij[ids.to.keep]
#       all_names_in_connect_graph = unique(names(unlist(paths_connect_ij$vpath)))
#       if(length(all_names_in_connect_graph)==0){next}
#       g_temp_vertex = g_vertex[g_vertex$name %in% all_names_in_connect_graph,]
#       g_temp_edge = g_edge[g_edge$from %in% all_names_in_connect_graph & g_edge$to %in% all_names_in_connect_graph,]
#       g_temp = graph_from_data_frame(d = g_temp_edge, vertices = g_temp_vertex, directed = T)
#       g_list[[length(g_list)+1]] = g_temp
#     }
#   }
#   return(g_list)
# }