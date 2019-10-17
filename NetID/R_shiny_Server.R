## server ####
server <- function(input, output, session) {
  ## Filter graph based on intensity, mz range and rt range ####
  g_filter <- reactive({
    print("enter g_filter")
    filter_graph(g,
                 input$Peak_inten_range[1], input$Peak_inten_range[2],
                 mz_lb = input$mz_range[1], mz_ub = input$mz_range[2], 
                 rt_lb = input$rt_range[1], rt_ub = input$rt_range[2])
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
  
  ## Partner graph ####
  g_partner <- reactive({
    print("enter g_partner")
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
  
  
  ## select graph of interest ####
  g_interest <- reactiveValues(data = NULL)
  observeEvent(input$one_node_graph, {
    print("enter g_interest search_partner")
    ### basic graph based on selected peak and formula
    g_interest$data = search_partner(g_partner(), 
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
  
  g_aesthetic <- reactive({
    print("enter g_interest aesthetic_filter")
    aesthetic_filter(g_interest$data,
                     show_metabolite_group = input$is_metabolite, 
                     show_artifact_edge = input$show_artifact_edges,
                     show_library_node = input$show_library_nodes)
  })
  
  
  ## output graph, related node & edge table ####
  output$network_proxy_nodes <- renderVisNetwork({
    print("enter Plot_g_interest")
    Plot_g_interest(g_aesthetic(), isolate(input$Peak_id), isolate(input$formula_select))
  })
  
  output$nodetable <- renderDataTable({
    print("enter show g_interest_vertice")
    g_interest_vertice = g_show_vertice(g_interest$data)
  })
  output$edgetable <- renderDataTable({
    print("enter show g_interest_edge")
    g_interest_edge = g_show_edge(g_interest$data)
  })
}
