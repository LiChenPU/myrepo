## server ##
server <- function(input, output, session) {
  
  
  # navbarPage1 - Peak list ####
  {
    peak_list <- reactive({
      print("enter show_peak_list")
      show_peak_list(ilp_nodes,
                     mz_interest = input$mz_interest, 
                     ion_form = input$ion_form, 
                     mz_ppm = input$mz_ppm,
                     input$Peak_inten_range[1], input$Peak_inten_range[2],
                     mz_lb = input$mz_range[1], mz_ub = input$mz_range[2],
                     rt_lb = input$rt_range[1], rt_ub = input$rt_range[2])
      
    })
    
    output$peak_list <- DT::renderDataTable({
      peak_list()
    }, options = list(pageLength = 10),
    filter = 'bottom',
    rownames = FALSE
    )
  }
  
  
  # navbarPage2 - Network visualization ####
  {
    ## selected peak ####
    node_selected <- reactive({
      print("enter node_selected")
      node_selected = ilp_nodes %>%
        filter(node_id == input$peak_id) %>%
        arrange(-ilp_result, -cplex_score)
      if(input$optimized_only){
        node_selected = node_selected %>%
          slice(1)
      }
      node_selected
    })
    
    query_ilp_id = reactive({
      print("enter query_ilp_id")
      node_selected() %>%
        filter(formula == input$formula) %>%
        filter(class == input$class) %>%
        pull(ilp_node_id)
    })
    
    observe({
      x = node_selected()
      peak_formula = x$formula
      updateSelectInput(session, "formula",
                        choices = peak_formula,
                        selected = head(peak_formula, 1)
      )
    })
    
    observe({
      x = input$formula
      peak_class = node_selected() %>% filter(formula == x) %>% pull(class)
      updateSelectInput(session, "class",
                        choices = peak_class,
                        selected = head(peak_class, 1)
      )
    })
    
    output$peak_mz <- renderText({
      paste0("m/z ", node_selected()$medMz[1])
    })
    
    output$peak_rt <- renderText({
      paste0("RT ", node_selected()$medRt[1])
    })
    
    output$peak_inten <- renderText({
      paste0("log10_inten ", node_selected()$log10_inten[1])
    })
    
    
    ## Network graph ####
    {
      g_parent <- reactive({
        print("enter g_parent")
        g_interest = NULL
        if(input$class == "Metabolite"){
          g_interest = network_annotation_met(query_ilp_id = query_ilp_id(),
                                              g_annotation = g_met,
                                              core_ilp_node = core_met,
                                              optimized_only = input$optimized_only)
        }
        if(input$class == "Artifact"){
          g_interest = network_annotation_nonmet(query_ilp_id = query_ilp_id(), 
                                                 g_annotation = g_nonmet, 
                                                 core_ilp_node = core_nonmet, 
                                                 weight_tol = 1,
                                                 optimized_only = input$optimized_only)
        }
        g_interest
      })
      
      g_child <- reactive({
        print("enter g_child")
        query_ilp_id = node_selected() %>%
          filter(formula == input$formula) %>%
          filter(class == input$class) %>%
          pull(ilp_node_id)
        g_child_artifact = network_child_nonmet(query_ilp_id = query_ilp_id(), 
                                                g_annotation = g_nonmet,
                                                connect_degree = 1,
                                                optimized_only = input$optimized_only)
      })
      
      g_interest <- eventReactive(input$plot_network, {
        print("enter g_interest")
        merge_nodes = NULL
        merge_edges = NULL
        if(input$parent_graph){
          g1_nodes = igraph::as_data_frame(g_parent(), "vertices")
          g1_edges = igraph::as_data_frame(g_parent(), "edges")
          merge_nodes = bind_rows(merge_nodes, g1_nodes)
          merge_edges = bind_rows(merge_edges, g1_edges)
        }
        
        if(input$child_graph){
          g2_nodes = igraph::as_data_frame(g_child(), "vertices")
          g2_edges = igraph::as_data_frame(g_child(), "edges")
          merge_nodes = bind_rows(merge_nodes, g2_nodes)
          merge_edges = bind_rows(merge_edges, g2_edges)
        }
        
        if(is.null(merge_edges)){
          return(NULL)
        }
        
        merge_nodes = merge_nodes %>%
          distinct()
        
        g_interest = graph_from_data_frame(merge_edges,
                                           directed = T,
                                           merge_nodes)
        
      })
      
      output$Network_plot <- renderVisNetwork({
        print("enter output$Network_plot renderVisNetwork")
        Plot_g_interest(g_interest(), 
                        query_ilp_node = isolate(query_ilp_id()), 
                        show_node_labels = input$node_labels, 
                        show_edge_labels = input$edge_labels
        )
      })
    }
    
    
    ## show structure_table ####
    {
      structure_table_trigger <- reactiveVal()
      observeEvent(input$plot_network, {
        structure_table_trigger("plot_network")
      })
      observeEvent(input$click, {
        structure_table_trigger("click")
      })
      
      # Use reactiveVal + observeEvant to update internal value
      structure_table <- reactive({
        print("enter structure_table")
        req(structure_table_trigger())
        if(structure_table_trigger() == "plot_network"){
          print(query_ilp_id())
          core_rank %>%
            filter(ilp_node_id == query_ilp_id())
        } else if(structure_table_trigger() == "click"){
          core_rank %>%
            filter(ilp_node_id == input$click)
        }
      })
      
      output$structure_list <- DT::renderDataTable({
        print("enter output$structure_list")
        structure_table() %>%
          dplyr::select(class, core_annotate, origin, note)
      }, options = list(pageLength = 5),
      filter = 'bottom',
      autoHideNavigation = T,
      rownames = FALSE
      )
    }
    
    
    ## show structure_plot ####
    {
      structure_plot_counter = reactiveVal(1)
      
      # goto next smiles when click
      observeEvent(input$struct_plot_click, {
        structure_plot_counter(structure_plot_counter()+1)
      })
      
      # goto previous smiles when dblclick
      observeEvent(input$struct_plot_dblclick, {
        req(structure_plot_counter()>1)
        structure_plot_counter(structure_plot_counter()-1)
      })
      
      # reset when new structure table is trigger
      observeEvent(structure_table(), {
        structure_plot_counter(1)
      })
      
      output$structure <- renderPlot({
        print("enter output$structure")
        smiles = structure_table() %>%
          pull(SMILES)
        my_SMILES2structure(smiles[structure_plot_counter()])
      })
    }
    
    
   
    
  }
  
  
  
  
}





## old ####
# server <- function(input, output, session) {
#   ## Filter graph based on intensity, mz range and rt range ####
#   g_filter <- reactive({
#     print("enter g_filter")
#     filter_graph(g,
#                  input$Peak_inten_range[1], input$Peak_inten_range[2],
#                  mz_lb = input$mz_range[1], mz_ub = input$mz_range[2], 
#                  rt_lb = input$rt_range[1], rt_ub = input$rt_range[2])
#   })
#   
#   ## adjust mz based on ionization and ppm
#   mz_interest <- reactive({
#     print("enter mz_interest")
#     ion_form = input$ion_form
#     if(ion_form == "M"){mz_adjust = input$mz_interest}
#     if(ion_form == "M+H"){mz_adjust = input$mz_interest - 1.007276}
#     if(ion_form == "M-H"){mz_adjust = input$mz_interest + 1.007276}
#     mz_adjust
#   })
#   
#   ## show peak table that meet the mz requirement
#   peak_table <- reactive({
#     print("enter peak_table")
#     search_peak(g_filter(), mz_interest = mz_interest(), mz_ppm = input$mz_ppm)
#   })
#   
#   output$peak_table <- DT::renderDataTable({
#     print("output peak_table")
#     peak_table()
#   })
#   
#   ## update Peak_id
#   observe({
#     print("enter update peak_id")
#     x<-peak_table()
#     if (nrow(x)==0){
#       value <- numeric(0)
#     } else{
#       value <- x$ID[1]
#     }
#     updateNumericInput(session, "Peak_id",
#                       value = value
#     )
#   })
#   
#   ## update formula_select
#   observe({
#     print("enter update formula_select")
#     x = search_formula(isolate(g_filter()), input$Peak_id)
#     if (is.null(x)){
#       x <- " "
#     } else{
#       x <- x$formula
#     }
#     updateSelectInput(session, "formula_select",
#                       # label = paste("Formula", length(x)),
#                       choices = x,
#                       selected = head(x, 1)
#     )
#   })
#   
#   ## Partner graph ####
#   g_partner <- reactive({
#     print("enter g_partner")
#     search_partner(isolate(g_filter()), isolate(input$Peak_id), input$formula_select, 5)
#   })
#   
#   ## Partner table
#   Partner_table <- reactive({
#     print("enter Partner_table")
#     g_show_vertice_rankIntensity(g_partner())
#   })
#   
#   
#   ## Output partner table
#   output$Partner_table <- DT::renderDataTable({
#     print("output Partner_table")
#     Partner_table()
#     })
#   
#   observe({
#     print("enter update Partner_id")
#     # observeEvent(Partner_table(), {
#     x = Partner_table()
#     if (is.null(x)){
#       x <- numeric(0)
#     } else{
#       x <- x$ID[1]
#     }
#     
#     updateNumericInput(session, "Partner_id",
#                        value = x
#     )
#   })
#   
#   observe({
#     print("enter update Partner_formula")
#     x = search_formula(isolate(g_filter()), input$Partner_id)
#     if (is.null(x)){
#       x <- character(0)
#     } else{
#       x <- x$formula
#     }
#     updateSelectInput(session, "Partner_formula",
#                       # label = paste("Formula", length(x)),
#                       choices = x,
#                       selected = head(x, 1)
#     )
#   })
#   
#   
#   ## select graph of interest ####
#   g_interest <- reactiveValues(data = NULL)
#   observeEvent(input$one_node_graph, {
#     print("enter g_interest search_partner")
#     ### basic graph based on selected peak and formula
#     g_interest$data = search_partner(g_partner(), 
#                                      isolate(input$Peak_id), 
#                                      isolate(input$formula_select), 
#                                      isolate(input$depth))
#   })
# 
#   observeEvent(input$two_nodes_shortest, {
#     print("enter g_interest two_formula_shortest_path_graph")
#     g_interest$data = two_formula_shortest_path_graph(isolate(g_partner()), 
#                                                       isolate(input$Peak_id), 
#                                                       isolate(input$Partner_id), 
#                                                       isolate(input$formula_select), 
#                                                       input$Partner_formula)
#   })
#   
#   observeEvent(input$two_nodes_all_graph, {
#     print("enter g_interest two_formula_neighbor_graph")
#     g_interest$data = two_formula_neighbor_graph(isolate(g_partner()), 
#                                                  isolate(input$Peak_id), 
#                                                  isolate(input$Partner_id), 
#                                                  isolate(input$formula_select), 
#                                                  input$Partner_formula, 
#                                                  dist = input$Partner_level)
#   })
#   
#   g_aesthetic <- reactive({
#     print("enter g_interest aesthetic_filter")
#     aesthetic_filter(g_interest$data,
#                      isolate(input$Peak_id), 
#                      isolate(input$formula_select), 
#                      show_metabolite_group = input$is_metabolite, 
#                      show_artifact_edges = input$show_artifact_edges,
#                      show_biotransform_edges = input$show_biotransform_edges,
#                      show_duplicated_formulas = input$show_duplicated_formulas,
#                      show_library_nodes = input$show_library_nodes)
#   })
#   
#   
#   ## output graph, related node & edge table ####
#   output$network_proxy_nodes <- renderVisNetwork({
#     print("enter Plot_g_interest")
#     Plot_g_interest(g_aesthetic(),
#                     show_metabolite_labels = input$show_metabolite_labels, 
#                     show_artifact_labels = input$show_artifact_labels,
#                     show_biotransform_edge_labels = input$show_biotransform_edge_labels, 
#                     show_artifact_edge_labels = input$show_artifact_edge_labels)
#   })
#   
#   output$nodetable <- DT::renderDataTable({
#     print("enter show g_interest_vertice")
#     g_interest_vertice = g_show_vertice(g_interest$data)
#   })
#   output$edgetable <- DT::renderDataTable({
#     print("enter show g_interest_edge")
#     g_interest_edge = g_show_edge(g_interest$data)
#   })
# }
