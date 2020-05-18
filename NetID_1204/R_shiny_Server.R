## server ##
server <- function(input, output, session) {
  
  
  # top part - Peak list ####
  {
    peak_list <- reactive({
      print("enter show_peak_list")
      show_peak_list(ilp_nodes,
                     input_interest = input$input_interest, 
                     ion_form = input$ion_form, 
                     mz_ppm = input$mz_ppm)
      
    })
    
    output$peak_list <- DT::renderDataTable({
      peak_list()
    }, options = list(pageLength = 5,
                      aoColumnDefs = list(list(sClass="alignLeft"))
                      ),
    
    filter = 'bottom',
    rownames = FALSE
    )
  }
  
  
  # bottom part - Network visualization ####
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
    
    observeEvent(peak_list(), {
      x = peak_list() %>% slice(1)
      updateNumericInput(session, "peak_id",
                         label = "Peak ID",
                         value = x$peak_id
      )
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
      fillContainer = F,
      rownames = TRUE
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





