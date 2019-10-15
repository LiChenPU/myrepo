library(shiny)
library(igraph)
library(reactlog)
library(ShinyTester)
# install.packages("visNetwork")
library(visNetwork)
library(dplyr)
library(RColorBrewer)
library(stringr)


# options(shiny.reactlog=TRUE) 
# Read in files ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
datapath = ("./Xi_new_neg")
setwd(datapath)


g_vertex = read.csv("g_vertex.txt", stringsAsFactors = F) %>%
  mutate(mz = round(mz, 4)) %>%
  mutate(RT = round(RT, 2)) %>% 
  mutate(intensity = signif(intensity, 6))
g_edge = read.csv("g_edge.txt", stringsAsFactors = F)

g <- graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = F)
g_vertex = igraph::as_data_frame(g, "vertices")
g_edge = igraph::as_data_frame(g, "edges")

# test_g = interest_node_graph(g, peak_id = 122, formula_select = "C4H7K3O8S1", step=1)
test_g = interest_node_graph(g, peak_id = 127, formula_select = "C20H32N6O12S2", step=1)

display.brewer.all()
my_palette = brewer.pal(4, "Set3")

# test_g = interest_node_graph(g, peak_id = 178, formula_select = "C6H12O6", step=1)
nodes = igraph::as_data_frame(test_g, "vertices") %>%
  # dplyr::select(name) %>%
  mutate(id = name) %>%
  mutate(label = ifelse(is_metabolite == "Yes", formula, "")) %>%
  # mutate(group = is_metabolite) %>%
  mutate(size = ifelse(is.na(intensity), 5, log10(intensity)) * 2) %>%
  mutate(color = case_when(
    is_metabolite == "Yes" ~ my_palette[1],
    is_metabolite == "Maybe" ~ my_palette[2],
    is_metabolite == "No" ~ "AAAAAA",
    is.na(is_metabolite) ~ "AAAAAA"
  )) %>%
  # [:digit:] means 0-9, \\. means ".", + means one or more, \\1 means the content in (), <sub> is HTML language
  mutate(formula_sub = str_replace_all(formula,"([[:digit:]|\\.|-]+)","<sub>\\1</sub>")) %>%
  mutate(title = paste0("Formula:", formula_sub, "<br>",
                        "ID:", ID, "<br>",
                        "mz:", mz, "<br>",
                        "RT:", RT, "<br>",
                        "TIC:", intensity)
  )

edges = igraph::as_data_frame(test_g, "edges") %>%
  # dplyr::select(from, to) %>%
  mutate(arrows = ifelse(direction==-1, "from", "to")) %>%
  # mutate(length = 100) %>%
  mutate(label = ifelse(category == "biotransform", linktype, "")) %>%
  mutate(color = case_when(
    category == "biotransform" ~ my_palette[1],
    category != "biotransform" ~ "AAAAAA"
  )) %>%
  mutate(title = paste0(category, "<br>",
                        ifelse(direction==-1, paste0("-",linktype), linktype), "<br>"
  )
  )

visNetwork(nodes, edges, height = "100%", width = "100%") %>% 
  visLegend() %>%
  visOptions(manipulation = TRUE, 
             highlightNearest = TRUE, 
             nodesIdSelection = TRUE,
             selectedBy = "is_metabolite"
  ) %>%
  visGroups(groupname = "Yes", color = my_palette[1]) %>%
  visGroups(groupname = "Maybe", color = my_palette[2]) %>%
  visGroups(groupname = "No", color = my_palette[3]) %>%
  visInteraction(navigationButtons = TRUE)

