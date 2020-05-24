
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('R_shiny_functions.R')
# options(shiny.reactlog=TRUE) # Enable reactlog here
# options(digits = 8)

# Read in files ####
datapath = "./WL_liver_neg"
setwd(datapath)
filename = "20200519181142_output.rds"


# Global parameter in the backgroun ####
{
  dt = readRDS(filename)
  
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
  reactlogReset()
  source('R_shiny_functions.R')
  source('R_shiny_UI.R', local = TRUE)
  source('R_shiny_Server.R')
  app = shinyApp(ui = ui, server = server)
  runApp(app)
  
  shiny::reactlogShow() # Run after closing the app
}

test = dt$ilp_nodes %>%
  filter(node_id == 1802)
  filter(category == "Library_MS2_fragment")
