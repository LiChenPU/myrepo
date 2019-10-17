
## ui ####
ui <- shinyUI({
  fluidPage(
  ## sidebarPanel ####
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
                                      choices = c("Yes", "No", "Maybe", NA), 
                                      selected = c("Yes", "No", "Maybe", NA)
                   ),
                   checkboxInput("show_artifact_edges", "show_artifact_edges",
                                 value = T
                   ),
                   checkboxInput("show_library_nodes", "show_library_nodes",
                                 value = T
                   )
                )
                
    ),
    actionButton("one_node_graph", "one_node_graph"),
    actionButton("two_nodes_shortest", "two_nodes_shortest"),
    actionButton("two_nodes_all_graph", "two_nodes_all_graph")
    
  ),
  
  ## mainPanel ####
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Peak_table", dataTableOutput("peak_table")),
                tabPanel("Partner_table", dataTableOutput("Partner_table")),
                tabPanel("Network_plot",  
                         visNetworkOutput("network_proxy_nodes", height = "400px")
                ),
                tabPanel("Nodes", dataTableOutput("nodetable")),
                tabPanel("Edges", dataTableOutput("edgetable"))
    )
  )
)
})

