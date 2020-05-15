
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('R_shiny_functions.R')
source("NetID_function.R")

# options(shiny.reactlog=TRUE) # Enable reactlog here
# options(digits = 8)

# Read in files ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
datapath = ("./Xi_new_neg")
setwd(datapath)

# Prepare 
{
  dt = readRDS("output.rds")
  
  g_met = initiate_g_met(dt$ilp_nodes, dt$ilp_edges)
  
  core_met = dt$ilp_nodes %>%
    filter(steps == 0) %>%
    filter(class == "Metabolite")
  
  core_nonmet = dt$ilp_nodes %>%
    filter(steps %% 1 == 0) %>%
    filter(class != "Unknown")
  
  g_nonmet = initiate_g_nonmet(dt$ilp_nodes, dt$ilp_edges, dt$heterodimer_ilp_edges)
  
  ilp_edges_annotate_met = igraph::as_data_frame(g_met, "edges")
  ilp_edges_annotate_nonmet = igraph::as_data_frame(g_nonmet, "edges")
  
  test = dt$ilp_edges %>%
    filter(category == "Library_MS2_fragment")
}


# function
{
  # given a metabolite ilp
  # show all shortest connections of from core
  # for each core, show all possible formulas
  g_met_interest = network_annotation_met(query_ilp_id = 34, 
                                          g_annotation = g_met, 
                                          core_ilp_node = core_met,
                                          optimized_only = T)
  Plot_g_interest(g_met_interest, query_ilp_node = 34)
  
  g_artifact_interest = network_annotation_nonmet(query_ilp_id = 523,  #2863
                                                  g_annotation = g_nonmet, 
                                                  core_ilp_node = core_nonmet, 
                                                  weight_tol = 3,
                                                  optimized_only = T)
  Plot_g_interest(g_artifact_interest, query_ilp_node = 523)
  
  g_child_artifact = network_child_nonmet(query_ilp_id = 2, 
                                          g_annotation = g_nonmet,
                                          connect_degree = 1,
                                          optimized_only = T)
  Plot_g_interest(g_child_artifact, query_ilp_node = 2)
  
  
  core_rank = core_annotate(dt$ilp_nodes, dt$FormulaSet_df, dt$LibrarySet)
  
  

}





## Run shiny ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('R_shiny_UI.R', local = TRUE)
source('R_shiny_Server.R')
shinyApp(ui = ui, server = server)


