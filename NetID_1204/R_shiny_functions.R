# Shiny library ####
{
  library(shiny)
  library(shinyjs)
  library(igraph)
  library(reactlog)
  library(ShinyTester)
  library(shinythemes)
  library(visNetwork)
  library(dplyr)
  library(RColorBrewer)
  library(stringr)
  library(ChemmineR)
  library(ChemmineOB)
  library(enviPat)
  library(lc8)
}

# function ####

## show_peak_list ####
show_peak_list = function(ilp_nodes, input_interest, ion_form, mz_ppm){
  
  ilp_nodes_filter = ilp_nodes %>%
    filter(ilp_result > 1e-6) %>%
    dplyr::select(node_id, medMz, medRt, log10_inten, class, formula, ppm_error) %>%
    dplyr::rename(peak_id = node_id)
  
  mz_adjust = 0
  if(!is.na(as.numeric(input_interest))){
    mz_interest = as.numeric(input_interest)
    if(ion_form == "M"){mz_adjust = mz_interest}
    if(ion_form == "M+H"){mz_adjust = mz_interest - 1.007276}
    if(ion_form == "M-H"){mz_adjust = mz_interest + 1.007276}
  } else {
    input_interest = gsub(" ", "", input_interest)
    formula_check = enviPat::check_chemform(isotopes, input_interest)
    if(!formula_check$warning){
      mz_adjust = formula_mz(formula_check$new_formula)
    } 
  }
    
  if(mz_adjust != 0){
    ilp_nodes_filter = ilp_nodes_filter %>%
      filter(abs(medMz - mz_adjust) < medMz * mz_ppm * 1e-6)
  }
  
  return(ilp_nodes_filter)
}

## initiate_g_met ####
initiate_g_met = function(ilp_nodes, ilp_edges){
  ilp_nodes_met = ilp_nodes %>%
    filter(grepl("Metabolite", class)) %>%
    arrange(-ilp_result)
  ilp_edges_met = ilp_edges %>%
    filter(category == "Biotransform") %>%
    mutate(direction = ifelse(ilp_nodes1 < ilp_nodes2, 1, -1)) %>%
    dplyr::select(ilp_nodes1, ilp_nodes2, everything())
  
  g_met = graph_from_data_frame(ilp_edges_met, 
                                directed = T, 
                                vertices = ilp_nodes_met)
  
}

## initiate_g_nonmet ####
initiate_g_nonmet = function(ilp_nodes, ilp_edges, heterodimer_ilp_edges){
  
  ilp_edges_merge = merge(ilp_edges, heterodimer_ilp_edges, all = T) 
  
  ilp_edges_nonmet = ilp_edges_merge %>%
    filter(category != "Biotransform") %>%
    dplyr::select(ilp_nodes1, ilp_nodes2, everything())
  
  ilp_edges_nonmet_reorder1 = ilp_edges_nonmet %>%
    filter(direction == 1) %>%
    dplyr::rename(from = ilp_nodes1, to = ilp_nodes2)
  
  ilp_edges_nonmet_reorder2 = ilp_edges_nonmet %>%
    filter(direction == -1) %>%
    dplyr::rename(from = ilp_nodes2, to = ilp_nodes1)
  
  ilp_edges_nonmet_directed = bind_rows(ilp_edges_nonmet_reorder1, ilp_edges_nonmet_reorder2)
  
  ilp_edges_weights = ilp_edges_nonmet_directed %>%
    mutate(edge_weight = 2 - ilp_result - 0.1*cplex_score) %>%
    mutate(edge_weight = ifelse(category == "Heterodimer", edge_weight + 1.5, edge_weight)) %>%
    arrange(edge_weight)
  
  g_nonmet = graph_from_data_frame(ilp_edges_weights, 
                                   directed = T, 
                                   vertices = ilp_nodes)
  
  return(g_nonmet)
}

## core_annotate ####
core_annotate = function(ilp_nodes, FormulaSet_df, LibrarySet){
  # Make annotation to core peaks (where steps == 0)
  
  core_rank = FormulaSet_df %>%
    # filter(category == "Metabolite") %>%
    filter(steps == 0) %>%
    merge(LibrarySet %>%
            dplyr::select(library_id, name, note, origin, SMILES, status) %>% 
            dplyr::rename(parent_id = library_id)) %>%
    mutate(score_source = case_when(
      origin == "Manual_library" ~ 0.5,
      origin == "Library_known" ~ 0.5,
      status == "quantified" ~ 0.6,
      status == "detected" ~ 0.3,
      status == "expected" ~ 0,
      status == "predicted" ~ -0.5
    )) %>% 
    mutate(score_source = ifelse(transform == "", score_source, score_source-0.7)) %>%
    mutate(rank_score = score_source + prior_score - database_prior) %>%
    arrange(-rank_score) %>% 
    mutate(core_annotate = case_when(
      transform == "" ~ paste(name, parent_formula),
      direction == 1 ~ paste(name, parent_formula, "+", transform),
      direction == -1 ~ paste(name, parent_formula, "-", transform)
    )) %>%
    mutate(class = case_when(
      category == "Unknown" ~ "Unknown",
      category == "Metabolite" ~ "Metabolite", 
      T ~ "Artifact"
    )) %>%
    inner_join(ilp_nodes %>%
                 dplyr::select(node_id, formula, category, ilp_node_id))
  
  return(core_rank)
}
## network_annotation_met ####  
network_annotation_met = function(query_ilp_id, 
                                  g_annotation = g_met, 
                                  core_ilp_node = core_met,
                                  optimized_only = T){
  
  # g_annotation = g_met
  # query_ilp_id = 2
  # core_ilp_node = core_met
  
  if(optimized_only){
    g_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(ilp_result != 0)
    g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(from %in% g_nodes$name, to %in% g_nodes$name)
    
    g_annotation = graph_from_data_frame(g_edges,
                                         directed = T,
                                         vertices = g_nodes)
    
    core_ilp_node = core_ilp_node %>%
      filter(ilp_result != 0)
  }
  
  core_ilp_id = as.character(core_ilp_node$ilp_node_id)
  query_ilp_id = as.character(query_ilp_id)
  
  valid_ilp_ids = igraph::as_data_frame(g_annotation, "vertices") %>% pull(name)
  
  if(!any(query_ilp_id == valid_ilp_ids)){
    return(NULL)
  }
  
  query_distMatrix <- shortest.paths(g_annotation,
                                     v=core_ilp_id,
                                     to=query_ilp_id,
                                     mode = "all")
  
  query_distMatrix_min = min(query_distMatrix[query_distMatrix!=0])
  
  if(is.infinite(query_distMatrix_min)){
    sub_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(name == query_ilp_id)
    
    sub_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(FALSE)
    
    g_met_interest = graph_from_data_frame(sub_edges,
                                           directed = T,
                                           vertices = sub_nodes)
    return(g_met_interest)
  }
  
  shortest_ilp_nodes = which(query_distMatrix == query_distMatrix_min)
  parent_selected = core_ilp_id[shortest_ilp_nodes]
  paths_connect_ij_nodes = lapply(parent_selected, function(x){
    shortest_paths(g_annotation, 
                   from = x, 
                   to = query_ilp_id, 
                   mode = "all", 
                   output = "vpath")
  })
  
  ilp_node_path = lapply(paths_connect_ij_nodes, function(x){
    x[[1]] %>% unlist() %>% names() %>% as.numeric()
  })
  
  ilp_node_interest = unlist(ilp_node_path)
  
  sub_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
    filter(name %in% as.character(ilp_node_interest))
  
  sub_edges = igraph::as_data_frame(g_annotation, "edges") %>%
    filter(from %in% as.character(ilp_node_interest) & 
             to %in% as.character(ilp_node_interest))
  
  g_met_interest = graph_from_data_frame(sub_edges,
                                         directed = F,
                                         vertices = sub_nodes)
  # plot.igraph(g_met_interest)
  
  return(g_met_interest)
}


## network_annotation_nonmet ####  
network_annotation_nonmet = function(query_ilp_id, 
                                     g_annotation = g_nonmet, 
                                     core_ilp_node = core_nonmet, 
                                     weight_tol = 1,
                                     optimized_only = T){
  
  # query_ilp_id = 6030
  # g_annotation = g_nonmet
  # core_ilp_node = core_nonmet
  # weight_tol = 1
  # optimized_only = T
  
  if(optimized_only){
    g_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(ilp_result != 0)
    # g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
    #   filter(ilp_result != 0)
    g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(from %in% g_nodes$name, to %in% g_nodes$name) %>%
      arrange(-ilp_result) %>%
      distinct(from, to, .keep_all=T) %>%
      filter(ilp_result>1e-6)
    
    g_annotation = graph_from_data_frame(g_edges, 
                                         directed = T, 
                                         vertices = g_nodes)
    
    core_ilp_node = core_ilp_node %>%
      filter(ilp_result != 0)
  }
  
  core_ilp_id = as.character(core_ilp_node$ilp_node_id)
  query_ilp_id = as.character(query_ilp_id)
  
  valid_ilp_ids = igraph::as_data_frame(g_annotation, "vertices") %>% pull(name)
  
  if(!any(query_ilp_id == valid_ilp_ids)){
    return(NULL)
  }
  
  
  g_edge = igraph::as_data_frame(g_annotation, "edges")
  
  query_distMatrix <- shortest.paths(g_annotation,
                                     v=core_ilp_id,
                                     to=query_ilp_id,
                                     mode = "out",
                                     weights = g_edge$edge_weight)
  
  query_distMatrix_min = min(query_distMatrix[query_distMatrix!=0])
  
  if(is.infinite(query_distMatrix_min)){
    sub_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(name == query_ilp_id)
    
    sub_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(FALSE)
    
    g_interest = graph_from_data_frame(sub_edges,
                                       directed = F,
                                       vertices = sub_nodes)
    return(g_interest)
  }
  
  shortest_ilp_nodes = which(query_distMatrix < query_distMatrix_min+weight_tol)
  # shortest_ilp_nodes = which(query_distMatrix == query_distMatrix_min)
  parent_selected = core_ilp_id[shortest_ilp_nodes]
  paths_connect_ij_nodes = lapply(parent_selected, function(x){
    all_shortest_paths(g_annotation, 
                       from = x, 
                       to = query_ilp_id, 
                       mode = "out")
  })
  
  ilp_node_path = lapply(paths_connect_ij_nodes, function(x){
    x[[1]] %>% unlist() %>% names() %>% as.numeric()
  })
  
  ilp_node_interest = unlist(ilp_node_path)
  
  sub_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
    filter(name %in% as.character(ilp_node_interest))
  
  sub_edges = igraph::as_data_frame(g_annotation, "edges") %>%
    filter(from %in% as.character(ilp_node_interest) & 
             to %in% as.character(ilp_node_interest))
  
  g_interest = graph_from_data_frame(sub_edges,
                                     directed = T,
                                     vertices = sub_nodes)
  # plot.igraph(g_interest)
  return(g_interest)
}

## network_child_nonmet ####  
network_child_nonmet = function(query_ilp_id, 
                                g_annotation = g_nonmet,
                                connect_degree = 1, 
                                optimized_only = T){
  
  # query_ilp_id = 3
  # g_annotation = g_nonmet
  # core_ilp_node = core_nonmet
  
  query_ilp_id = as.character(query_ilp_id)
  
  if(optimized_only){
    g_nodes = igraph::as_data_frame(g_annotation, "vertices") %>%
      filter(ilp_result != 0)
    # g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
    #   filter(ilp_result != 0)
    g_edges = igraph::as_data_frame(g_annotation, "edges") %>%
      filter(from %in% g_nodes$name, to %in% g_nodes$name) %>%
      arrange(-ilp_result) %>%
      distinct(from, to, .keep_all=T) %>%
      filter(ilp_result > 0.01)
    
    g_annotation = graph_from_data_frame(g_edges, 
                                         directed = T, 
                                         vertices = g_nodes)
  }
  
  g_partner = make_ego_graph(g_annotation, 
                             connect_degree,
                             nodes = query_ilp_id, 
                             mode = c("out"))[[1]]
  
  # plot.igraph(g_partner)
  return(g_partner)
}

## Plot_g_interest ####
Plot_g_interest = function(g_interest, query_ilp_node, show_node_labels = T, show_edge_labels = T)
{
  
  # query_ilp_node = 871
  # show_node_labels = T
  # show_edge_labels = T
  
  if(!is.igraph(g_interest)){return(NULL)}
  
  my_palette = brewer.pal(4, "Set3")
  nodes = igraph::as_data_frame(g_interest, "vertices") %>%
    dplyr::rename(id = name) %>%
    mutate(size = log10_inten * 2) %>%
    mutate(label = "") %>%
    mutate(color = case_when(
      # assigned to the first color, not overwitten by later assignment
      id == as.character(query_ilp_node) ~ my_palette[4],
      class == "Metabolite" ~ my_palette[1],
      class == "Artifact" ~ my_palette[3],
      class == "Unknown" ~ "#666666"
    )) %>%
    # [:digit:] means 0-9, \\. means ".", + means one or more, \\1 means the content in (), <sub> is HTML language
    mutate(formula_sub = str_replace_all(formula,"([[:digit:]|\\.|-]+)","<sub>\\1</sub>")) %>%
    mutate(title = paste0("Formula:", formula_sub, "<br>",
                          "ID:", node_id, "<br>",
                          "mz:", round(medMz,4), "<br>",
                          "RT:", round(medRt,2), "<br>",
                          "TIC:", format(signif(10^log10_inten,5), scientific = T))
    )
  
  
  
  if(show_node_labels){
    nodes = nodes %>%
      mutate(label = formula)
  }
  
  edges = igraph::as_data_frame(g_interest, "edges") %>%
    mutate(arrows = "to") %>%
    mutate(label = "") %>%
    mutate(color = "#666666") %>%
    filter(T)
  
  if(show_edge_labels & nrow(edges)!=0){
    edges = edges %>%
      mutate(label = case_when(
        category == "Multicharge" ~ paste0(linktype, "-charge"),
        category == "Oligomer" ~ paste0(linktype, "-mer"),
        category == "Heterodimer" ~ paste0("Peak ", linktype),
        category == "Library_MS2_fragment" ~ linktype,
        direction == -1 ~ paste0("-", linktype),
        T ~ linktype
      ))
  }
  
  visNetwork(nodes, edges) %>% 
    visLegend() %>%
    visOptions(manipulation = TRUE, 
               highlightNearest = TRUE
               # height = "100%", width = "100%"
    ) %>%
    visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
              ;}") %>%
    visInteraction(navigationButtons = TRUE) %>%
    visIgraphLayout(layout = "layout_with_fr", 
                    randomSeed = 123)
  
}
## my_SMILES2structure ####
my_SMILES2structure = function(SMILES){
  SDF = try(smiles2sdf(SMILES), silent=T)
  if(inherits(SDF, "try-error")){
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    # text(x=0.5, y=0.5, paste("Invalid SMILES"))
    return(0)
  }
  tryError = try(ChemmineR::plotStruc(SDF[[1]]), silent=T)
  if(inherits(tryError, "try-error")){
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    # text(x=0.5, y=0.5, paste("Invalid SDF"))
    return(0)
  }
  return(1)
  
}