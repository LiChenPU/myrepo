
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('R_shiny_functions.R')

# options(digits = 8)

# Read in files ####
datapath = "./WL_liver_pos/200711"
# datapath = "./WL_liver_neg/200711"
datapath = "./Unknown in IO RT/Sc_neg"
setwd(datapath)
# filename = "20200704005808_output.rds" # Sc_pos
filename = "20200709234041_output.rds" # liver_neg
filename = "20200711005313_output.rds" # liver_pos
filename = "20200704085443_output.rds" # Sc_neg

testing = F
# Global parameter in the background ####
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
  # reactlogReset()
  # options(shiny.reactlog=TRUE) # Enable reactlog here
  
  source('R_shiny_functions.R')
  source('R_shiny_UI.R', local = TRUE)
  source('R_shiny_Server.R')
  app = shinyApp(ui = ui, server = server)
  runApp(app)
  
  # shiny::reactlogShow() # Run after closing the app
}


