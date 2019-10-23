# Load library ####
{
  library(shiny)
  library(igraph)
  library(reactlog)
  library(ShinyTester)
  # install.packages("visNetwork")
  library(visNetwork)
  library(dplyr)
  library(RColorBrewer)
  library(stringr)
}

# options(digits = 8)
# 
# test_edge = g_edge %>% filter(ILP_result !=0, grepl("\\[37", linktype))
# test_node = g_vertex %>% filter(ID == 1156)
# {
#   select_node_id = 177
#   select_formula_ls = formula_list %>%
#     filter(ID == select_node_id)
# 
# }
# 
# {
#   select_ilp_id = 14051
#   select_ilp_id = 359
#   select_edge_info_sum = g_edge %>%
#     filter(from == select_ilp_id | to == select_ilp_id)
# 
# }
# 
# {
#   pred_formula_ls = CPLEXset$formula$pred_formula_ls
#   pred_formula_ls[[1797]]
#   pred_formula_ls[[5109]]
# }


# options(shiny.reactlog=TRUE) # Enable reactlog here
# Read in files ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
datapath = ("./Xi_new_neg")
setwd(datapath)

g_vertex = read.csv("g_vertex2.txt", stringsAsFactors = F) %>%
  mutate(mz = round(mz, 4)) %>%
  mutate(RT = round(RT, 2)) %>%
  mutate(intensity = signif(intensity, 6))
g_edge = read.csv("g_edge2.txt", stringsAsFactors = F)

g <- graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = T)
g_vertex = igraph::as_data_frame(g, "vertices")
g_edge = igraph::as_data_frame(g, "edges")

# test_g = search_partner(g, peak_id = 122, formula_select = "C4H7K3O8S1", step=1)
# test_g = search_partner(g, peak_id = 127, formula_select = "C20H32N6O12S2", step=1)
# test_g = search_partner(g, peak_id = 178, formula_select = "C6H12O6", step=1)


## Run shiny ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('R_shiny_functions.R')
source('R_shiny_UI.R', local = TRUE)
source('R_shiny_Server.R')
shinyApp(ui = ui, server = server)


