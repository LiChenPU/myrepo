# function ####
## filter_graph ####
filter_graph = function(g, 
                        intensity_lb = 3, intensity_ub= 7, 
                        mz_lb = 0, mz_ub = 1500, 
                        rt_lb = 0, rt_ub = 30)
{
  # intensity_lb = 3
  # intensity_ub= 7
  # mz_lb = 0
  # mz_ub = 1500
  # rt_lb = 0
  # rt_ub = 30
  
  # Select which peaks to retain, and include all library nodes
  g_vertex = igraph::as_data_frame(g, "vertices") %>%
    filter((ILP_result != 0 & !is.na(ILP_result)) | RT == -1) %>%
    filter(is.na(intensity) | (intensity > 10^intensity_lb & intensity < 10^intensity_ub)) %>%
    filter(mz > mz_lb, mz < mz_ub) %>%
    filter((RT > rt_lb & RT < rt_ub) | RT == -1) 
  
  metabolite_id = g_vertex %>%
    filter(Biotransform) %>%
    pull(name)
  
  g_edge = igraph::as_data_frame(g, "edges") %>%
    filter(from %in% g_vertex$name, 
           to %in% g_vertex$name) %>%
    # Filter out biotransform connect between artifacts
    filter(category != "biotransform" | (from %in% metabolite_id & to %in% metabolite_id))
  
  g_filter = graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = T)
  
  return(g_filter)
}


## search_peak ####
search_peak = function(g, mz_interest, mz_ppm)
{
  # mz_interest = 180.0631
  # mz_ppm = 5
  
  g_vertex = igraph::as_data_frame(g, "vertices") %>%
    filter(mz < mz_interest*(1+mz_ppm/10^6),
           mz > mz_interest*(1-mz_ppm/10^6)) %>%
    arrange(-intensity, abs(mz-mz_interest), -ILP_result, abs(cal_mass - mz_interest)) %>%
    dplyr::select("ID", "mz", "RT", "formula", "compound_name", "intensity", "is_metabolite")
  
  return(g_vertex)
}

## search_formula ####
search_formula = function(g, peak_id){
  if(is.na(peak_id)){return(NULL)}
  g_vertex = igraph::as_data_frame(g, "vertices") %>%
    filter(ID==peak_id) %>%
    arrange(-ILP_result, -cplex_score) %>%
    dplyr::select("ID", "mz", "RT", "formula", "compound_name", "intensity", "is_metabolite")
  
  return(g_vertex)
  
}

## g_show_vertice_rankIntensity ####
g_show_vertice_rankIntensity = function(g_interest)
{
  
  if(!is_igraph(g_interest)) {return(NULL)}
  g_interest_vertice = igraph::as_data_frame(g_interest, "vertices")
  g_interest_vertice = g_interest_vertice[with(g_interest_vertice, order(-intensity)),]
  
  colname_vertice = c("ID", "mz", "RT", "formula", "compound_name", "intensity", "is_metabolite")
  g_interest_vertice = g_interest_vertice[,colname_vertice]
  return(g_interest_vertice)
}




## search_partner ####
search_partner = function(g_filter, peak_id, formula_select, step = 1)
{
  # peak_id = 178
  # formula_select = "C6H12O6"
  # step = 1
  # g_filter = filter_graph(g)
  
  if(!is_igraph(g_filter)) {return(NULL)}
  
  g_vertex = igraph::as_data_frame(g_filter, "vertices") 
  g_edge = igraph::as_data_frame(g_filter, "edges") 
  g_id = g_vertex %>%
    filter(ID==peak_id, formula == formula_select) %>%
    pull(name)
  
  if(length(g_id) == 0){
    # print("return g_id ==0")
    return(NULL)}
  
  g_partner <- make_ego_graph(g_filter, 
                              step, 
                              nodes = as.character(g_id[1]), 
                              mode = c("all"))[[1]]
  
  g_partner_vertex = igraph::as_data_frame(g_partner, "vertices") 
  
  g_partner_edge = igraph::as_data_frame(g_filter, "edges") %>%
    filter(from %in% g_partner_vertex$name, 
           to %in% g_partner_vertex$name) 
  
  g_partner = graph_from_data_frame(g_partner_edge, vertices = g_partner_vertex, directed = T)
  return(g_partner)
}

## aesthetic_filter ####
aesthetic_filter = function(g_partner,
                            peak_id,
                            formula_select,
                            show_metabolite_group, 
                            show_artifact_edges = T,
                            show_biotransform_edges = T,
                            show_library_nodes = T,
                            show_duplicated_formulas = T)
{
  
  # show_metabolite_group = c("No", NA)
  # show_artifact_edges = F
  # show_library_nodes = F
  # show_biotransform_edges = T
  # show_duplicated_formulas = F
  if(!is.igraph(g_partner)){return (NULL)}
  
  g_interest_vertex = igraph::as_data_frame(g_partner, "vertices") 
  
  center_vertex = g_interest_vertex %>%
    filter(ID == peak_id, formula == formula_select)
  
  show_metabolite_group[show_metabolite_group == ""] = NA
  g_interest_vertex = g_interest_vertex %>%
    filter(is_metabolite %in% show_metabolite_group) 
  
  if(!show_library_nodes){
    g_interest_vertex = g_interest_vertex %>%
      filter(RT != -1)
  }
  
  if(!show_duplicated_formulas){
    g_interest_vertex_library = g_interest_vertex %>%
      filter(RT == -1)
    g_interest_vertex = g_interest_vertex %>%
      arrange(-intensity) %>%
      distinct(formula, .keep_all=T) %>%
      rbind(g_interest_vertex_library) %>%
      distinct()
  }
  
  g_interest_vertex = g_interest_vertex %>%
    rbind(center_vertex) %>%
    distinct()
  
  g_interest_edges = igraph::as_data_frame(g_partner, "edges") %>%
    filter(from %in% g_interest_vertex$name, 
           to %in% g_interest_vertex$name) 
  
  if(!show_artifact_edges){
    g_interest_edges = g_interest_edges %>%
      filter(category == "biotransform")
  }
  
  if(!show_biotransform_edges){
    g_interest_edges = g_interest_edges %>%
      filter(category != "biotransform")
  }
  
  g_interest_vertex = g_interest_vertex %>%
    filter(name %in% g_interest_edges$from | name %in% g_interest_edges$to)

  
  
  g_aesthetic_filter = graph_from_data_frame(g_interest_edges, vertice = g_interest_vertex, directed = T)
  
  return(g_aesthetic_filter)
  
}

## Plot_g_interest ####
Plot_g_interest = function(g_interest, show_metabolite_labels = T, show_artifact_labels = F,
                           show_biotransform_edge_labels = T, show_artifact_edge_labels = F)
{
  # browser()
  # interested_node=178
  # formula_select="C6H12O6"
  # g_interest = g_temp
  # display.brewer.all()
  # display.brewer.pal(4, "Set3")
  # g_interest = search_partner(g, peak_id = 178, formula_select = "C6H12O6", step=1)
  
  # g_interest = g_temp
  if(!is.igraph(g_interest)){return(NULL)}
  
  my_palette = brewer.pal(4, "Set3")
  nodes = igraph::as_data_frame(g_interest, "vertices") %>%
    mutate(id = name) %>%
    mutate(size = ifelse(is.na(intensity), 5, log10(intensity)) * 2) %>%
    mutate(label = "") %>%
    mutate(color = case_when(
      # assigned to the first color, not overwitten by later assignment
      is_metabolite == "No" ~ my_palette[3],
      RT == -1 ~ my_palette[4],
      is_metabolite == "Yes" ~ my_palette[1],
      is_metabolite == "Maybe" ~ my_palette[2],
      is.na(is_metabolite) ~ my_palette[3]
      )) %>%
    # [:digit:] means 0-9, \\. means ".", + means one or more, \\1 means the content in (), <sub> is HTML language
    mutate(formula_sub = str_replace_all(formula,"([[:digit:]|\\.|-]+)","<sub>\\1</sub>")) %>%
    mutate(title = paste0("Formula:", formula_sub, "<br>",
                          "ID:", ID, "<br>",
                          "mz:", mz, "<br>",
                          "RT:", RT, "<br>",
                          "TIC:", intensity)
    )
  
  if(show_metabolite_labels){
    nodes = nodes %>% 
      mutate(label = ifelse(is_metabolite %in% c("Yes", "Maybe"), formula, label))
  }
  if(show_artifact_labels){
    nodes = nodes %>% 
      mutate(label = ifelse(is_metabolite %in% c("Maybe", "No", "NA"), formula, label))
  }
  
  edges = igraph::as_data_frame(g_interest, "edges") %>%
    mutate(arrows = ifelse(direction==-1, "from", "to")) %>%
    mutate(label = "") %>%
    mutate(color = case_when(
      category == "biotransform" ~ my_palette[1],
      category != "biotransform" ~ my_palette[3]
    )) %>%
    mutate(title = paste0(category, "<br>",
                          ifelse(direction==-1, paste0("-",linktype), linktype), "<br>"
    )
    )
  
  if(show_biotransform_edge_labels){
    edges = edges %>% 
      mutate(label = ifelse(category == "biotransform", linktype, label))
  }
  if(show_artifact_edge_labels){
    edges = edges %>% 
      mutate(label = ifelse(category != "biotransform", linktype, label))
  }
  
  
  visNetwork(nodes, edges, height = "200%", width = "200%") %>% 
    visLegend() %>%
    visOptions(manipulation = TRUE, 
               highlightNearest = TRUE, 
               nodesIdSelection = TRUE,
               selectedBy = "is_metabolite"
    ) %>%
    visGroups(groupname = "Yes", color = my_palette[1]) %>%
    visGroups(groupname = "Maybe", color = my_palette[2]) %>%
    visGroups(groupname = "No", color = my_palette[3]) %>%
    visInteraction(navigationButtons = TRUE) %>%
    visPhysics(stabilization = FALSE) %>%
    visLayout(randomSeed = 123)
    # visLayout(hierarchical = TRUE) 
    # visIgraphLayout(layout = "layout_in_circle")
  
}


## g_show_vertice & g_show_edge####
g_show_vertice = function(g_interest)
{
  if(!is_igraph(g_interest)) {return(NULL)}
  g_interest_vertice = igraph::as_data_frame(g_interest, "vertices")
  colname_vertice = c("ID", "mz", "RT", "formula", "compound_name","intensity", "is_metabolite")
  g_interest_vertice = g_interest_vertice[,colname_vertice]
  return(g_interest_vertice)
}

g_show_edge = function(g_interest)
{
  if(!is_igraph(g_interest)) {return(NULL)}
  g_interest_edge = igraph::as_data_frame(g_interest, "edges")
  colname_edge = c("node1", "node2", "formula1", "formula2", "linktype")
  g_interest_edge = g_interest_edge[,colname_edge]
  return(g_interest_edge)
}


## two_formula_neighbor_graph####
two_formula_neighbor_graph = function(g, node1, node2, formula1, formula2, dist)
{
  # node1 = 178
  # node2 = 2
  # formula1 = "C6H12O6"
  # formula2 = "C16H32O2"
  # dist = 3
  if(!is_igraph(g)) {return(NULL)}
  
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_id1 = g_vertex$name[g_vertex$ID==node1 & g_vertex$formula==formula1]
  g_id2 = g_vertex$name[g_vertex$ID==node2 & g_vertex$formula==formula2]
  
  i = g_id1[1]
  j = g_id2[1]
  
  if(distances(g, i, j, mode = "all") > dist){return(NULL)}
  
  all_names_in_connect_graph = c()
  for(count in 0:dist){
    g_i = make_ego_graph(g, 
                         count, 
                         nodes = i, 
                         mode = c("all"))[[1]]
    g_j = make_ego_graph(g, 
                         dist-count, 
                         nodes = j, 
                         mode = c("all"))[[1]]
    
    g_i_vertex =  igraph::as_data_frame(g_i, "vertices")
    g_j_vertex =  igraph::as_data_frame(g_j, "vertices")
    
    intersect_names = intersect(g_i_vertex$name, g_j_vertex$name)
    all_names_in_connect_graph = c(all_names_in_connect_graph, intersect_names)
  }
  
  all_names_in_connect_graph = unique(all_names_in_connect_graph)
  
  if(length(all_names_in_connect_graph)==0){return(NULL)}
  g_temp_vertex = g_vertex %>% 
    filter(name %in% all_names_in_connect_graph)
  g_temp_edge = g_edge %>%
    filter(from %in% all_names_in_connect_graph & to %in% all_names_in_connect_graph)
  g_temp = graph_from_data_frame(d = g_temp_edge, vertices = g_temp_vertex, directed = T)
  
  return(g_temp)
}

## shortest_two_nodes_shortest####
two_formula_shortest_path_graph = function(g, node1, node2, formula1, formula2)
{
  # node1 = 178
  # node2 = 2
  # formula1 = "C6H12O6"
  # formula2 = "C16H32O2"
  
  if(!is_igraph(g)) {return(NULL)}
  
  g_vertex = igraph::as_data_frame(g, "vertices")
  g_id1 = g_vertex$name[g_vertex$ID==node1 & g_vertex$formula==formula1]
  g_id2 = g_vertex$name[g_vertex$ID==node2 & g_vertex$formula==formula2]
  
  i = g_id1[1]
  j = g_id2[1]
  
  paths_connect_ij = shortest_paths(g,i,j, mode = "all")
  # ids.to.keep = sapply(paths_connect_ij, function(i) length(i)<=dist)
  # paths_connect_ij = paths_connect_ij[ids.to.keep]
  all_names_in_connect_graph = unique(names(unlist(paths_connect_ij$vpath)))
  
  if(length(all_names_in_connect_graph)==0){return(NULL)}
  g_temp_vertex = g_vertex[g_vertex$name %in% all_names_in_connect_graph,]
  g_temp_edge = g_edge[g_edge$from %in% all_names_in_connect_graph & g_edge$to %in% all_names_in_connect_graph,]
  g_temp = graph_from_data_frame(d = g_temp_edge, vertices = g_temp_vertex, directed = T)
  
  return(g_temp)
}