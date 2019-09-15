source("~/myrepo/mz_calculator/mz_calculator_function.R")

library(readr)
library(igraph)
library(ggplot2)
library(ggrepel)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Read HMDB_files & library_data
{
 
  HMDB_clean = read_tsv('hmdb_structure_sdf_unique_neutral_formula.tsv')
  HMDB_clean2 = HMDB_clean %>%
    dplyr::mutate(mz = Exact_Mass, ID = 1:nrow(.))
  
  {
    data("isotopes")
    HMDB_xml = read_csv("hmdb_xml_checkchemform.csv") 
    
    HMDB_quantified_detected_unique = HMDB_xml %>% 
      filter(status %in% c("detected", "quantified")) %>%
      distinct(chemical_formula, .keep_all=T)
    
    HMDB_quantified_detected_formula = HMDB_xml %>%
      filter(chemical_formula %in% HMDB_quantified_detected_unique$chemical_formula)
    
    HMDB_clean_quantified_detected = HMDB_clean %>%
      filter(HMDB_ID %in% HMDB_quantified_detected_formula$accession)
    
    HMDB_clean_quantified_detected2 = HMDB_clean_quantified_detected %>%
      dplyr::mutate(mz = Exact_Mass, ID = 1:nrow(.))
    
    
    HMDB_clean2 = HMDB_clean_quantified_detected2
  }

  
  library_data = expand_formula_to_library("C6H12O6")
  
}

{
  Mset = list()
  Mset[["Raw_data"]] = HMDB_clean2
  Mset[["Library"]] = library_data
  Mset[["Connect_rules"]]=Read_rule_table(rule_table_file = 'biotransform2.csv', extend_rule = F)
  Mset[["Global_parameter"]]=  list(mode = 0)
  Mset[["Data"]] = Peak_cleanup(Mset)
}


{
  Mset[["NodeSet"]]=Form_node_list(Mset)
  EdgeSet = list()
  
  EdgeSet[["Connect_rules"]] = Edge_Connect_rules(Mset, 
                                                  mass_abs = 0.0002, 
                                                  mass_ppm = .5/10^6)
  
  EdgeSet[["Connect_rules"]] = Edge_score(EdgeSet$Connect_rules, plot_graph = T)
  EdgeSet[["Merge"]] = Merge_edgeset(EdgeSet)
  
  Mset[["NodeSet_network"]] = Network_prediction(Mset, 
                                                 edge_biotransform = EdgeSet$Connect_rules,
                                                 biotransform_step = 100,
                                                 propagation_score_threshold = 0.9,
                                                 top_n = 50,
                                                 print_step = T)
  
  CPLEXset = Prepare_CPLEX(Mset, EdgeSet)
}


### evaluate proposed formula and 
{
  all_formula =bind_rows(CPLEXset$data$pred_formula_ls)
  
  merge_HMDB_propagation = HMDB_clean2 %>%  
    mutate(ID = 1:nrow(.)) %>%
    merge(all_formula, by.x = "ID", by.y = "id", all=T )
  
  all_formula_unique_num = table(table(all_formula$id))
  all_formula_unique = all_formula %>% distinct(id, .keep_all=T)
  
  all_formula_highscore = all_formula %>%
    filter(score>0.9)
  
  merge_HMDB_propagation_correct = merge_HMDB_propagation %>% filter(MF == formula)
  merge_HMDB_propagation_wrong = merge_HMDB_propagation %>% filter(!ID %in% merge_HMDB_propagation_correct$ID, !is.na(formula))
  
  basic_linktype = c("H2", "C1H2", "H2O1","H1N1", "H3N1", "C16H30O1", "H1O3P1", "C5H8", "H2S1", 
                     "C1O2", "C2H4", "C6H10O5", "O1", "C2H2O1", "C6H4O1")
}


##Network analysis


# identify if there is large cluster ignored, 
# 2 > 20
{
  merge_node_list = Mset$NodeSet %>% 
    merge(HMDB_clean2, by = "ID", all=T)
  
  merge_edge_list = EdgeSet$Connect_rules %>% filter(edge_massdif_score > 0.9) 
  g <- graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list, directed = FALSE)
  clu = components(g)
  edge_size=table(clu$csize)
  mainnetwork = igraph::groups(clu)[table(clu$membership)>2000]
  subnetwork = igraph::groups(clu)[table(clu$membership)>50&table(clu$membership)<300]
  
  sub_nodelist = HMDB_clean2[subnetwork[[1]],]
  g_subnetwork = make_ego_graph(subnetwork[[1]][1], graph = g, order = diameter(g), mode = "all")[[1]]
  
  merge_node_list["degree"]=degree(g, mode = c("all"))
  g <- graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list, directed = FALSE)
  # pdf("full_all_connect_graph.pdf")
  # plot.igraph(g,
  #      vertex.color = 'cyan',
  #      vertex.label = "",
  #      frame.color = "grey",
  #      # vertex.label.color = "black",
  #      # vertex.label.cex = 1,
  #      edge.color = 'grey',
  #      # edge.label = edge.attributes(g_subnetwork)$linktype,
  #      vertex.size = sqrt(vertex.attributes(g)$degree)+0.5,
  #      edge.arrow.size = .05,
  #      layout = layout_with_graphopt(g)
  #      # main = paste("Subnetwork of node", interested_node)
  # )
  # dev.off()
}

# evaluate from basic 10 biotransform, adding which one has most gain
{
  network_size_linktype = list()
  i = unique(EdgeSet$Connect_rules$linktype)[1]
  for(i in unique(EdgeSet$Connect_rules$linktype)){
    selected_linktype = c(basic_linktype, i)
    merge_node_list = merge_HMDB_propagation_correct
    merge_edge_list = EdgeSet$Connect_rules %>%
      filter(node1 %in% merge_node_list$ID, node2 %in% merge_node_list$ID) %>%
      filter(linktype %in% selected_linktype) %>%
      filter(edge_massdif_score > 0.9)
    merge_node_list = merge_node_list %>% 
      filter(ID %in% merge_edge_list$node1|ID %in% merge_edge_list$node2)
    # table(merge_edge_list$linktype)
    #Generate network graph using node_list and edge_list above
    
    g <- graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list, directed = FALSE)
    clu = components(g)
    edge_size=max(as.numeric(names(table(clu$csize))))
    network_size_linktype[[length(network_size_linktype)+1]] = data.frame(linktype=i, max_size = edge_size, stringsAsFactors = F)
  }
  
  summary_network_gain = bind_rows(network_size_linktype) %>%
    mutate(gain = max_size - min(max_size))
  
}


# evaluate if there is biotransform that is not importnat
{
  merge_node_list = merge_HMDB_propagation_correct
  merge_edge_list = EdgeSet$Connect_rules %>%
    filter(node1 %in% merge_node_list$ID, node2 %in% merge_node_list$ID) %>%
    filter(linktype %in% basic_linktype) %>%
    filter(edge_massdif_score > 0.9)
  merge_node_list = merge_node_list %>% 
    filter(ID %in% merge_edge_list$node1|ID %in% merge_edge_list$node2)
  
  #Generate network graph using node_list and edge_list above
  g <- graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list, directed = FALSE)
  clu = components(g)
  edge_size=table(clu$csize)
  mainnetwork = igraph::groups(clu)[table(clu$membership) == max(table(clu$membership))]
  g_main = make_ego_graph(g, order=diameter(g), nodes = mainnetwork[[1]][1], mode = c("all"))[[1]]
  
  main_node_list = igraph::as_data_frame(g_main, "vertices")
  main_edge_list = igraph::as_data_frame(g_main, "edges")
  
  edge_size= edge_membership = list()
  i = "H2O1"
  for (i in unique(main_edge_list$linktype)){
    temp_edge_list=main_edge_list[main_edge_list$linktype!=i,]
    temp_merge_node_list = main_node_list %>% filter(name %in% temp_edge_list$from|name %in% temp_edge_list$to)
    temp_g <- graph_from_data_frame(d = temp_edge_list, vertices = temp_merge_node_list, directed = FALSE)
    temp_clu = components(temp_g)
    edge_size[[i]]=table(temp_clu$csize)
    edge_membership[[i]]=temp_clu$membership
  }

  drop_out_size=sapply(edge_size, function(x) max(as.numeric(names(x))))
  
  
  number_link = data.frame(table(main_edge_list$linktype), stringsAsFactors = F)
  drop_out_size2 = data.frame(Var1 = names(drop_out_size), drop_out_size=drop_out_size, stringsAsFactors = F) %>%
    mutate(drop_out_num = nrow(main_node_list) - drop_out_size)
  
  summary_linktype = merge(number_link,drop_out_size2)
}

# Plot step network
{
  merge_node_list = merge_HMDB_propagation_correct
  merge_edge_list = EdgeSet$Connect_rules %>%
    filter(node1 %in% merge_node_list$ID, node2 %in% merge_node_list$ID) %>%
    filter(linktype %in% basic_linktype) %>%
    filter(edge_massdif_score > 0.9)
  
  merge_node_list = merge_node_list %>% 
    filter(ID %in% merge_edge_list$node1|ID %in% merge_edge_list$node2)
  g <- graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list, directed = FALSE)
  
  clu = components(g)
  edge_size=table(clu$csize)
  mainnetwork = igraph::groups(clu)[table(clu$membership) == max(table(clu$membership))]
  g_main = make_ego_graph(g, order=diameter(g), nodes = mainnetwork[[1]][1], mode = c("all"))[[1]]
  
  main_node_list = igraph::as_data_frame(g_main, "vertices")
  main_edge_list = igraph::as_data_frame(g_main, "edges")
  
  step_info = as.data.frame(table(main_node_list$steps), stringsAsFactors = F) %>%
    transmute(Step = Var1, No.metabolites = Freq) %>%
    mutate(Step = as.numeric(Step)) %>%
    mutate(Total.metabolites = cumsum(No.metabolites))  
  step_info2 = step_info %>%
    filter(!Step>9) %>%
    add_row(Step = ">9", No.metabolites=sum(step_info$No.metabolites[step_info$Step>9]), Total.metabolites = max(step_info$Total.metabolites))
    
            
  
  library(RColorBrewer)
  display.brewer.all(n=30, exact.n=FALSE)
  # color_set1 = c(brewer.pal(9, "Set1"),rep("#999999", 50))
  # color_set3_3color = rep(brewer.pal(3, "Set3"), c(3,3,50))
  # color_set1_3color = rep(c(brewer.pal(3, "Set1"),"#999999"), c(1,1,1,50))
  # color_set2_3color = rep(c(brewer.pal(3, "Set2"),"#999999"), c(1,2,3,50))
  
  color_spectral = c(brewer.pal(9, "Spectral"), rep("#666666", 50))
  
  vertex_color = color_spectral[main_node_list$steps]
  
  # pdf("full_connect_graph_graphopt.pdf")
  # for(i in 1:1){
  #   plot.igraph(g_main,
  #               vertex.color = vertex_color,
  #               vertex.label = "",
  #               vertex.frame.color = "black",
  #               # vertex.label.color = "black",
  #               # vertex.label.cex = 1,
  #               edge.color = 'grey',
  #               # edge.label = edge.attributes(g_subnetwork)$linktype,
  #               vertex.size = sqrt(degree(g_main, mode = c("all")))+0.5,
  #               # vertex.size = 2.5,
  #               # edge.arrow.size = .05,
  #               layout = layout_with_graphopt(g_main)
  #               # main = paste("Subnetwork of node", interested_node)
  #   )
  # }
  # dev.off()
  
  step_info2$Step = as.factor(step_info2$Step)
  step_info2$Step  = factor(step_info2$Step, levels=c(levels(step_info2$Step)[-1],levels(step_info2$Step)[1]))
  ggplot(step_info2, aes(x = Step)) +
    geom_bar(aes(y = No.metabolites), stat = 'identity', fill = color_spectral[1:nrow(step_info2)]) + 
    # scale_x_discrete(limits = temp_position) +
    geom_point(aes(y = Total.metabolites/4)) + 
    geom_line(aes(y = Total.metabolites/4), group=1) +
    geom_text_repel(aes(y = Total.metabolites/4, label = Total.metabolites), vjust = -.5) +
    scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Total.metabolites"), limits = c(0,2500), expand = c(0,0)) +
    
    # theme(legend.title = element_blank()) + 
    theme(legend.position = "right") +
    theme_classic()
   
}
# 
# {
#   #subnetwork criteria 
#   clu = components(g)
#   main_network = igraph::groups(clu)[table(clu$membership)>1000]
#   subnetwork = igraph::groups(clu)[table(clu$membership)>30&table(clu$membership)<1000]
#   g_subnetwork_list = lapply(subnetwork, make_ego_graph, graph=g, order=diameter(g), mode="all")
#   sub_nodelist = merge_library[merge_library$ID %in% subnetwork[[2]],]
#   #write.csv(sub_nodelist, "maybe wrong entry in HMDB.csv")
#   
#   g_sub = graph_from_data_frame(d = edge_list_sub, vertices = merge_node_list[merge_node_list$ID %in% c(edge_list_sub[,1], edge_list_sub[,2]),], directed = FALSE)
#   colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")
#   V(g_sub)$color = colors[vertex.attributes(g_sub)$category+1]
#   E(g_sub)$color = colors[edge_list_sub$category+1]
#   
#   #Basic graph characteristics, distance, degree, betweeness
#   {
#     #distance
#     farthest_vertices(g) 
#     #Degree
#     g.degree <- degree(g, mode = c("all"))
#     table(g.degree)
#     hist(g.degree)
#     # which.max(g.degree)
#     node_list_0degree = node_list[g.degree==0,]
#     # node_list$Exact_Mass[g.degree==0]
#     # node_list$Exact_Mass[71]
#     # node_list$MF[71]
#     #Betweenness
#     g.b <- betweenness(g, directed = T)
#     colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")
#     V(g)$color <- colors[merge_node_list$category+2]
#     plot(g, 
#          vertex.label = NA,
#          edge.color = 'black',
#          vertex.size = sqrt(degree(g, mode = c("all")))+1,
#          edge.arrow.size = 0.05,
#          layout = layout_nicely(g))
#   }
#   
#   
#   #Analyze the network/subgraph of specific node
#   {
#     interested_node = "1941"
#     g_intrest <- make_ego_graph(g,2, nodes = interested_node, mode = c("all"))[[1]]
#     dists = distances(g_intrest, interested_node)
#     colors <- c("black", "red", "orange", "blue", "dodgerblue", "cyan")
#     
#     #V(g_intrest)$color <- colors[dists+1]
#     plot(g_intrest,
#          #vertex.color = 'white',
#          #vertex.label = vertex.attributes(g_intrest)$Predict_formula,
#          #vertex.label = vertex.attributes(g_intrest)$MF,
#          vertex.label = vertex.attributes(g_intrest)$MF,
#          #vertex.label = vertex.attributes(g_intrest)$ID,
#          vertex.label.color = "red",
#          vertex.label.cex = 1,
#          edge.color = 'black',
#          edge.label = fun_group_1$fun_group[edge.attributes(g_intrest)$linktype],
#          vertex.size = 10,
#          edge.arrow.size = .05,
#          main = paste("Subnetwork of node", interested_node)
#     )
#   }
#   
#   #merge_node_list[which(merge_node_list$MF=="C6H13O4PS2"),]
#   
#   #Analysis of Subnetwork 
#   {
#     clu=components(g)
#     #table(clu$csize)
#     #subnetwork criteria 
#     main_network = igraph::groups(clu)[table(clu$membership)>1000]
#     subnetwork = igraph::groups(clu)[table(clu$membership)>30&table(clu$membership)<1000]
#     g_subnetwork_list = lapply(subnetwork, make_ego_graph, graph=g, order=diameter(g), mode="all")
#     sub_nodelist = merge_library[merge_library$ID %in% subnetwork[[1]],]
#     #write.csv(sub_nodelist, "maybe wrong entry in HMDB.csv")
#     
#     
#     for (i in 1:length(subnetwork)){
#       plot(g_subnetwork_list[[i]][[1]],
#            #vertex.color = 'white',
#            vertex.label = vertex.attributes(g_subnetwork_list[[i]][[1]])$MF,
#            #vertex.label = vertex.attributes(g_subnetwork_list[[i]][[1]])$Exact_Mass,
#            vertex.label.color = "red",
#            vertex.label.cex = 1,
#            edge.color = 'black',
#            edge.label = fun_group_1$fun_group[edge.attributes(g_subnetwork_list[[i]][[1]])$linktype],
#            vertex.size = 10,
#            edge.arrow.size = .05,
#            main = paste("Subnetwork",names(subnetwork)[[i]])
#       )
#     }
#     merge_node_list[subnetwork[["52"]],]
#   }
#   
# }





mz_calculator = function(raw_data, 
                         library_data,
                         connect_depth = 5,
                         ion_mode = 1,
                         rule_table_file = "~/myrepo/mz_calculator/dependent/connect_rules.csv"
                         
)
{
  # raw_data = test_rawdata
  # library_data = expand_formula_to_library("C5H7N3O")
  {
    Mset = list()
    Mset[["Raw_data"]] = raw_data
    Mset[["Library"]] = library_data
    Mset[["Connect_rules"]]=Read_rule_table(rule_table_file = rule_table_file)
    Mset[["Global_parameter"]]=  list(mode = ion_mode)
    Mset[["Data"]] = Peak_cleanup(Mset)
  }
  
  {
    Mset[["NodeSet"]]=Form_node_list(Mset)
    EdgeSet = list()
    
    EdgeSet[["Connect_rules"]] = Edge_Connect_rules(Mset, 
                                                    mass_abs = 0.002, 
                                                    mass_ppm = 25/10^6)
    
    EdgeSet[["Connect_rules"]] = Edge_score(EdgeSet$Connect_rules)
    EdgeSet[["Merge"]] = Merge_edgeset(EdgeSet)
    
    Mset[["NodeSet_network"]] = Network_prediction(Mset, 
                                                   edge_biotransform = EdgeSet$Connect_rules,
                                                   biotransform_step = connect_depth,
                                                   propagation_score_threshold = 0.25,
                                                   top_n = 50)
    
    CPLEXset = Prepare_CPLEX(Mset, EdgeSet)
  }
  
  # Run CPLEX ####
  {
    CPLEXset$data$unknown_formula = Score_formula(Mset, CPLEXset,
                                                  rdbe=T, step_score=F)
    edge_info_sum = Score_edge_cplex(CPLEXset, edge_bonus = 0.1)
    obj_cplex = c(CPLEXset$data$unknown_formula$cplex_score, edge_info_sum$edge_score)
    
    CPLEXset[["Init_solution"]] = list(Run_CPLEX(CPLEXset, obj_cplex))
    # CPLEXset[["Screen_solution"]] = CPLEX_screen_edge(CPLEXset, edge_bonus_range = seq(-.6, -0.9, by=-0.1))
    CPLEXset[["Pmt_solution"]] = CPLEX_permutation(CPLEXset, n_pmt = 20, sd_rel_max = 0.2)
  }
  
  # Read CPLEX result ####
  {
    CPLEX_all_x = Read_CPLEX_result(CPLEXset$Init_solution)
    CPLEX_all_x = Read_CPLEX_result(CPLEXset$Pmt_solution)
    
    CPLEX_x = rowMeans(CPLEX_all_x,na.rm=T)
    #CPLEX_x = result_solution$x
  }
  
  {
    unknown_nodes = CPLEXset$data$unknown_nodes[,1:3]
    unknown_formula = CPLEXset$data$unknown_formula
    
    unknown_formula["ILP_result"] = CPLEX_x[1:nrow(unknown_formula)]
    unknown_formula_CPLEX = unknown_formula[unknown_formula$ILP_result !=0,]
    print(paste("pred formula num =", nrow(unknown_formula_CPLEX)))
    
    unknown_node_CPLEX = merge(unknown_nodes,unknown_formula_CPLEX,by.x = "ID", by.y = "id",all=T)
  }
  
  return(unknown_node_CPLEX)
  
} 

