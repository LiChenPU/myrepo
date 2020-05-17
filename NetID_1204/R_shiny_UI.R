
## ui ####
ui <- tagList(
  # shinythemes::themeSelector(),
  navbarPage(
    theme = shinytheme("cerulean"), # other theme can be viewed from themeSelector()
    "NetID",
    ## Peak list tab ####
    tabPanel("Peak list",
      ## sidebarPanel ####
      sidebarPanel(
        width = 3,
        sidebarLayout(
          wellPanel(
            numericInput(inputId = "mz_interest", 
                         label = "Enter a mz of interest",
                         value = 0
            ),
            selectInput(inputId = "ion_form", 
                        label = "Select ionization",
                        choices = c("M+H","M-H","M"),
                        selected = "M"
            ),
            numericInput(inputId = "mz_ppm", 
                         label = "ppm",
                         value = 3
            )
          ),
          
          wellPanel(
            sliderInput(inputId = "Peak_inten_range", 
                        label = "Peak Intensity (log10)",
                        min = 0, max = 10, step = 0.01,
                        value = c(0,10)),
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
    
      ## mainPanel ####
      mainPanel(
        tabPanel("peak_table", dataTableOutput("peak_list"))
      )
    ),
      
    ## Network visualization tab ####
    tabPanel("Network visualization",
             ## First row ####
             fluidRow(
               column(3,
                      numericInput(inputId = "peak_id",
                                   label = "Peak ID",
                                   value = 0
                      )
               ),
               column(3,
                      verbatimTextOutput("peak_mz")
               ),
               column(3,
                      verbatimTextOutput("peak_rt")
               ),
               column(3,
                      verbatimTextOutput("peak_inten")
               )
             ),
             
             ## Second row ####
             fluidRow(
               column(3, 
                      selectInput(inputId = "formula",
                                  label = "Formula",
                                  choices = character(0)
                      )
               ),
               column(3,
                      selectInput(inputId = "class",
                                  label = "Class",
                                  choices = character(0)
                      )
               ),
               column(2,
                      checkboxInput("optimized_only", "optimized_only",
                                    value = T)
               )
             ),
             hr(),
             ## Network + structure ####
             fluidRow(
               ## Network ####
               column(7,
                      visNetworkOutput("Network_plot", height = "100%", width = "100%"),
                      column(2,
                             checkboxInput("node_labels", "Node labels",
                                           value = T)
                      ),
                      column(2,
                             checkboxInput("edge_labels", "Edge labels",
                                           value = T)
                      ),
                      column(2,
                             checkboxInput("parent_graph", "Parent graph",
                                           value = T)
                      ),
                      column(2,
                             checkboxInput("child_graph", "Child graph",
                                           value = F)
                      ),
                      column(2, offset = 1,
                             actionButton("plot_network", "Plot network")
                      )
                      
               ),
               ## Structure ####
               column(5,
                      plotOutput("structure"),
                      dataTableOutput("structure_list")
                      
                      ## Structure here ##
                      )
             )
    )
  )
)




### old ######
# ui <- shinyUI({
#   fluidPage(
#   ## sidebarPanel ####
#   sidebarPanel(
#     tabsetPanel(type = "tabs",
#                 tabPanel("Setting",
#                   numericInput(inputId = "mz_interest",
#                                label = "Enter a mz of interest",
#                                value = 180.0633
#                   ),
#                   selectInput(inputId = "ion_form",
#                                label = "Select ionization",
#                               choices = c("M+H","M-H","M"),
#                               selected = "M"
#                   ),
#                   numericInput(inputId = "mz_ppm",
#                                label = "ppm",
#                                value = 5
#                   ),
#                   numericInput(inputId = "Peak_id",
#                                label = "Peak ID",
#                                value = numeric(0)
#                   ),
#                   selectInput(inputId = "formula_select",
#                               label = "Formula",
#                               choices = character(0)
#                               ),
# 
#                   wellPanel(
#                     numericInput(inputId = "Partner_level",
#                                  label = "Level of partner",
#                                  value = 2
#                     ),
#                     numericInput(inputId = "Partner_id",
#                                  label = "Partner ID",
#                                  value = numeric(0)
#                     ),
#                     selectInput(inputId = "Partner_formula",
#                                 label = "Partner formula",
#                                 choices = character(0)
#                     )
#                   )
# 
# 
# 
#                 ),
#                 tabPanel("Advanced",
# 
#                    numericInput(inputId = "depth",
#                                 label = "Graph depth",
#                                 value = 1
#                    ),
#                    checkboxGroupInput("is_metabolite", "Select whether to include:",
#                                       choices = c("Yes", "No", "Maybe", NA),
#                                       selected = c("Yes","Maybe")
#                    ),
#                    checkboxInput("show_biotransform_edges", "show_biotransform_edges",
#                                  value = T
#                    ),
#                    checkboxInput("show_artifact_edges", "show_artifact_edges",
#                                  value = F
#                    ),
#                    checkboxInput("show_duplicated_formulas", "show_duplicated_formulas",
#                                  value = F
#                    ),
#                    checkboxInput("show_library_nodes", "show_library_nodes",
#                                  value = F
#                    ),
#                    checkboxInput("show_metabolite_labels", "show_metabolite_labels",
#                                  value = T
#                    ),
#                    checkboxInput("show_artifact_labels", "show_artifact_labels",
#                                  value = F
#                    ),
#                    checkboxInput("show_biotransform_edge_labels", "show_biotransform_edge_labels",
#                                  value = T
#                    ),
#                    checkboxInput("show_artifact_edge_labels", "show_artifact_edge_labels",
#                                  value = F
#                    )
# 
# 
#                 )
# 
#     ),
#     actionButton("one_node_graph", "one_node_graph"),
#     actionButton("two_nodes_shortest", "two_nodes_shortest"),
#     actionButton("two_nodes_all_graph", "two_nodes_all_graph")
# 
#   ),
# 
#   ## mainPanel ####
#   mainPanel(
#     tabsetPanel(type = "tabs",
#                 tabPanel("Peak_table", dataTableOutput("peak_table")),
#                 tabPanel("Partner_table", dataTableOutput("Partner_table")),
#                 tabPanel("Network_plot",
#                          visNetworkOutput("network_proxy_nodes", height = "800px")
#                 ),
#                 tabPanel("Nodes", dataTableOutput("nodetable")),
#                 tabPanel("Edges", dataTableOutput("edgetable"))
#     )
#   )
# )
# })
# 
