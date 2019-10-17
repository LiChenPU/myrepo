library(shiny)
library(igraph)
library(reactlog)
library(ShinyTester)
# install.packages("visNetwork")
library(visNetwork)
library(dplyr)
library(RColorBrewer)
library(stringr)


# options(shiny.reactlog=TRUE) # Enable reactlog here
# Read in files ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
datapath = ("./Xi_new_neg")
setwd(datapath)

g_vertex = read.csv("g_vertex.txt", stringsAsFactors = F) %>%
  mutate(mz = round(mz, 4)) %>%
  mutate(RT = round(RT, 2)) %>% 
  mutate(intensity = signif(intensity, 6))
g_edge = read.csv("g_edge.txt", stringsAsFactors = F)

g <- graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = T)
g_vertex = igraph::as_data_frame(g, "vertices")
g_edge = igraph::as_data_frame(g, "edges")

# test_g = search_partner(g, peak_id = 122, formula_select = "C4H7K3O8S1", step=1)
# test_g = search_partner(g, peak_id = 127, formula_select = "C20H32N6O12S2", step=1)
# test_g = search_partner(g, peak_id = 178, formula_select = "C6H12O6", step=1)


# Shiny R --------------------####
## Run shiny ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('R_shiny_functions.R')
source('R_shiny_UI.R', local = TRUE)
source('R_shiny_Server.R')
shinyApp(ui = ui, server = server)


