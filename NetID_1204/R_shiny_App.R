
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('R_shiny_functions.R')

# options(digits = 8)

# Read in files ####
# datapath = "./WL_liver_pos/200711"
datapath = "./Unknown in IO RT/Sc_pos"
setwd(datapath)
filename = "20200704005808_output.rds"
# filename = "20200709234041_output.rds"
# filename = "20200704085443_output.rds"

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

## Testing ####
if(testing){
  
  RColorBrewer::display.brewer.all()
  my_palette = brewer.pal(5, "Set3")
  
  g_met2_node = igraph::as_data_frame(g_met, "vertices") %>%
    filter(ilp_result > 0.01) %>%
    mutate(color = ifelse(class == "Metabolite", my_palette[1], my_palette[2]))
  
  g_met2_edges = igraph::as_data_frame(g_met, "edges") %>%
    filter(from %in% g_met2_node$name & to %in% g_met2_node$name) %>%
    mutate(color = my_palette[1])
  
  g_met2 = graph_from_data_frame(g_met2_edges,
                                 vertices = g_met2_node,
                                 directed = F)
  
  g_nonmet2_node = igraph::as_data_frame(g_nonmet, "vertices") %>%
    filter(ilp_result > 0.01) %>%
    mutate(color = case_when(
      class == "Metabolite" ~ my_palette[1],
      class == "Putative Metabolite" ~ my_palette[2],
      class == "Artifact" ~ my_palette[3],
      class == "Unknown" ~ "#666666"
    ))
  g_nonmet2_edges = igraph::as_data_frame(g_nonmet, "edges") %>%
    filter(from %in% g_nonmet2_node$name & to %in% g_nonmet2_node$name,
           ilp_result > 0.01) %>%
    mutate(color = my_palette[3])
  g_nonmet2 = graph_from_data_frame(g_nonmet2_edges,
                                 vertices = g_nonmet2_node,
                                 directed = F)
  

  
  g_all_nodes = bind_rows(g_met2_node, g_nonmet2_node) %>%
    distinct(name, .keep_all = T) %>%
    filter(class != "Unknown")
  g_all_edges = bind_rows(g_met2_edges, g_nonmet2_edges) %>%
    distinct() %>%
    filter(from %in% g_all_nodes$name & to %in% g_all_nodes$name)
  
  g_all = graph_from_data_frame(g_all_edges,
                                directed = F,
                                vertices = g_all_nodes)
  
  
  
  test3 = g_all %>% # g_all
    # set_vertex_attr("color", value = my_palette[1]) %>%
    set_vertex_attr("size", value = 2) %>%
    set_vertex_attr("label", value = NA) %>%
    # set_edge_attr("color", value = my_palette[1]) %>%
    set_edge_attr("width", value = 1)
  
  
  g_sub=test3
  clu=components(g_sub)
  #subnetwork criteria 
  subnetwork = igraph::groups(clu)[table(clu$membership)>= 30]
  test4 = make_ego_graph(graph=g_sub, 
                         order=diameter(g_sub), 
                         nodes = subnetwork[[1]][1],
                         mode="all")[[1]]

  
  pdf("network_plot2.pdf",
      width = 10,
      height = 10)
  plot.igraph(test4,
              layout = layout_with_graphopt
              # layout = layout_with_fr
              )
  dev.off()
  
  

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


