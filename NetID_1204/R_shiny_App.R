
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('R_shiny_functions.R')
source("NetID_function.R")

# options(shiny.reactlog=TRUE) # Enable reactlog here
# shiny::showReactLog() # Run after closing the app
# options(digits = 8)

# Read in files ####
datapath = ("./Xi_new_neg")
setwd(datapath)

# Prepare 
{
  dt = readRDS("20200518002230_output.rds")
  
  ilp_nodes = dt$ilp_nodes %>%
    mutate(medMz = signif(medMz, 7),
           medRt = round(medRt, 2),
           log10_inten = round(log10_inten, 2),
           ppm_error = round(ppm_error, 2))
  
  g_met = initiate_g_met(ilp_nodes, dt$ilp_edges)
  
  core_met = ilp_nodes %>%
    filter(steps == 0) %>%
    filter(class == "Metabolite")
  
  core_nonmet = ilp_nodes %>%
    filter(steps %% 1 == 0) %>%
    filter(class != "Unknown")
  
  g_nonmet = initiate_g_nonmet(ilp_nodes, dt$ilp_edges, dt$heterodimer_ilp_edges)
  
  # ilp_edges_annotate_met = igraph::as_data_frame(g_met, "edges")
  # ilp_edges_annotate_nonmet = igraph::as_data_frame(g_nonmet, "edges")
  core_rank = core_annotate(ilp_nodes, dt$FormulaSet_df, dt$LibrarySet)
  
  data(isotopes)

}

## Run shiny ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  source("NetID_function.R")
  source('R_shiny_functions.R')
  source('R_shiny_UI.R', local = TRUE)
  source('R_shiny_Server.R')
  app = shinyApp(ui = ui, server = server)
  runApp(app)
}



#----####

# function ####
{
  # given a metabolite ilp
  # show all shortest connections of from core
  # for each core, show all possible formulas
  g_met_interest = network_annotation_met(query_ilp_id = 2, 
                                          g_annotation = g_met, 
                                          core_ilp_node = core_met,
                                          optimized_only = T)
  Plot_g_interest(g_met_interest, query_ilp_node = 2)
  
  g_artifact_interest = network_annotation_nonmet(query_ilp_id = 6030, 
                                                  g_annotation = g_nonmet, 
                                                  core_ilp_node = core_nonmet, 
                                                  weight_tol = 1,
                                                  optimized_only = T)
  Plot_g_interest(g_artifact_interest, query_ilp_node = 6030)
  
  g_child_artifact = network_child_nonmet(query_ilp_id = 2, 
                                          g_annotation = g_nonmet,
                                          connect_degree = 1,
                                          optimized_only = T)
  Plot_g_interest(g_child_artifact, query_ilp_node = 2)
  
  # test = layout_with_fr(g_child_artifact, grid = "grid")
  
  
  
  
  my_SMILES2structure(core_rank$SMILES[8])
  
  
  
  
  
  
}






