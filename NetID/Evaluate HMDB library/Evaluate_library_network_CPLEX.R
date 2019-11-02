source("~/myrepo/NetID/NetID_function.R")
library(readr)
library(igraph)
library(ggplot2)
library(ggrepel)

# display.brewer.pal(4, "Set3")
my_palette = c(brewer.pal(4, "Set3"), rep("#666666", 50))
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
result = list()
for(ppm_error in c(0.5, 1, 2, 5)){
  sigma = ppm_error
# Read HMDB_files & library_data ####
{
  HMDB_clean = read_tsv('hmdb_structure_sdf_unique_neutral_formula.tsv')
  ## Run these if use detected/quantified HMDB
  {
    HMDB_xml = read_csv("hmdb_xml_checkchemform.csv")

    HMDB_quantified_detected_unique = HMDB_xml %>%
      filter(status %in% c("detected", "quantified")) %>%
      distinct(chemical_formula, .keep_all=T)

    HMDB_quantified_detected_formula = HMDB_xml %>%
      filter(chemical_formula %in% HMDB_quantified_detected_unique$chemical_formula)

    HMDB_clean_quantified_detected = HMDB_clean %>%
      filter(HMDB_ID %in% HMDB_quantified_detected_formula$accession)


    HMDB_clean = HMDB_clean_quantified_detected
  }
  
  HMDB_clean2 = HMDB_clean %>%
    mutate(medMz = Exact_Mass, ID = 1:nrow(.)) %>%
    rename(Compound_name = Name) %>%
    dplyr::select(ID, HMDB_ID, Compound_name, Exact_Mass, MF, medMz) %>%
    mutate(medRt = 0, pseudo_inten = 10)
  
  library_formula_num = 1
  
  set.seed(random_seed)
  library_formula = sample(HMDB_clean2$MF, library_formula_num)
  if(library_formula_num==1){library_formula = "C6H12O6"} 
  library_data = expand_formula_to_library(library_formula)
  
}

## adding variation to absolute mz
{
  ppm_error = ppm_error 
  
  set.seed(random_seed )
  HMDB_clean2 = HMDB_clean2 %>%
    mutate(medMz = medMz * (1+ rnorm(nrow(.), 0, ppm_error)/10^6)) # add ppm error
  # mutate(mz = mz * (1+ rnorm(nrow(.), 0, 2)/10^3)) # add abs error
  
  HMDB_clean_test = HMDB_clean2 %>%
    mutate(ppm = (medMz-Exact_Mass)/Exact_Mass*10^6) %>%
    mutate(abs = medMz-Exact_Mass) 
  
  # HMDB_clean_test%>%
  #   arrange(abs) %>% dplyr::select(abs) %>% tail()
  
}


## NetID run scripts ####
{
  Mset = list()
  Mset[["Library"]] = library_data
  
  Mset[["Biotransform"]]=Read_rule_table(rule_table_file = biotransform_file)
  Mset[["Artifacts"]]=Read_rule_table(rule_table_file = artifact_file)
  
  Mset[["Raw_data"]] = HMDB_clean2
}

## Initialise ####
{
  Mset[["Global_parameter"]]=  list(mode = ion_mode,
                                    normalized_to_col_median = F)
  Mset[["Cohort"]]=Cohort_Info(Mset, first_sample_col_num = first_sample_col_num)
  print(Mset$Cohort)
  
  #Clean-up duplicate peaks 
  Mset[["Data"]] = Peak_cleanup(Mset,
                                ms_dif_ppm=0/10^6, 
                                rt_dif_min=0,
                                detection_limit=0,
                                remove_high_blank_ratio = 2,
                                first_sample_col_num = first_sample_col_num)
  print(c(nrow(Mset$Raw_data), nrow(Mset$Data)))
  
  Mset[["ID"]] = Mset$Data$ID
}

# Network ####
{
  EdgeSet = list()
  
  Mset[["NodeSet"]]=Form_node_list(Mset)
  
  EdgeSet[["Biotransform"]] = Edge_biotransform(Mset, 
                                                mass_abs = 0,
                                                mass_ppm = 10)
  
  mass_dist_sigma = sigma
  EdgeSet[["Biotransform"]] = Edge_score(EdgeSet$Biotransform, mass_dist_sigma = mass_dist_sigma)
  
  EdgeSet[["Merge"]] = Merge_edgeset(EdgeSet)
  
  Mset[["NodeSet_network"]] = Network_prediction(Mset, 
                                                 EdgeSet,
                                                 biotransform_step = 16,
                                                 artifact_step = 0,
                                                 propagation_score_threshold = 0.2,
                                                 propagation_intensity_threshold = 0,
                                                 max_formula_num = 1e6,
                                                 top_n = 50)
}

## CPLEX optimiazation ####
{
  CPLEXset = list()
  CPLEXset[["formula"]] = Prepare_CPLEX_formula(Mset, mass_dist_sigma = mass_dist_sigma,
                                                rdbe=F, step_score=F, iso_penalty_score=F)
  CPLEXset[["edge"]] = Prepare_CPLEX_edge(EdgeSet, 
                                          CPLEXset,
                                          edge_bonus = edge_bonus, 
                                          isotope_bonus = isotope_bonus,
                                          artifact_bonus = artifact_bonus)
  CPLEXset[["para"]] = Prepare_CPLEX_para(Mset, EdgeSet, CPLEXset)
}
# Run CPLEX
{
  CPLEXset[["Init_solution"]] = list(Run_CPLEX(CPLEXset, obj_cplex = CPLEXset$para$obj))
  # CPLEXset[["Pmt_solution"]] = CPLEX_permutation(CPLEXset, n_pmt = 10, sd_rel_max = 0.2)
  
  # Test_para_CPLEX(CPLEXset, obj_cplex = CPLEXset$para$obj, test_para =c(1,3))
}

# Read CPLEX result 
{
  CPLEX_all_x = Read_CPLEX_result(CPLEXset$Init_solution)
  # CPLEX_all_x = Read_CPLEX_result(CPLEXset$Pmt_solution)
  
  CPLEX_x = rowMeans(CPLEX_all_x,na.rm=T)
  
  {
    unknown_nodes = CPLEXset$formula$unknown_nodes[,1:3]
    unknown_formula = CPLEXset$formula$unknown_formula %>% mutate(ILP_result = CPLEX_x[1:nrow(.)])
    unknown_formula_CPLEX = unknown_formula %>% filter(ILP_result!=0 )
    edge_info_sum = CPLEXset$edge %>%
      mutate(ILP_result = 0)
    edge_info_sum = edge_info_sum %>%
      filter(CPLEX_score >=0) %>% 
      mutate(ILP_result = CPLEX_x[(nrow(unknown_formula)+1):length(CPLEX_x)]) %>%
      rbind(edge_info_sum) %>%
      distinct(edge_ilp_id, .keep_all = T) %>%
      arrange(edge_ilp_id)
    print(paste("pred formula num =", nrow(unknown_formula_CPLEX)))
  }
}



# Network analysis ####
# Figure 1B - formula network #### 
## how many metabolites can be connected?
## # of metabolites that are Unconnected, subnetwork, main_network?
{
  ## Identify all correct_connects
  {
    all_connects = EdgeSet$Merge %>% filter(node1 %in% 1:nrow(Mset$Data), node2 %in% 1:nrow(Mset$Data))
    
    node1_formula = HMDB_clean2$MF[all_connects$node1]
    node2_formula = HMDB_clean2$MF[all_connects$node2]
    linktype_formula = all_connects$linktype
    correct_log = rep(F, nrow(all_connects))
    for(i in 1:nrow(all_connects)){
      if(linktype_formula[i] == ""){next}
      correct_log[i] = my_calculate_formula(node1_formula[i], linktype_formula[i]) == node2_formula[i]
    }
    correct_connects = all_connects %>% filter(correct_log)
  }
  
  ## All network with correct linkage
  {
    merge_node_list = Mset$NodeSet %>% 
      merge(HMDB_clean2, by = "ID", all=T) %>%
      filter(category == 1)
    merge_edge_list = correct_connects %>%
      filter(node1 %in% merge_node_list$ID, node2 %in% merge_node_list$ID)
    
    g_all <- graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list, directed = FALSE)
    
    merge_node_list = igraph::as_data_frame(g_all, "vertices")
    merge_edge_list = igraph::as_data_frame(g_all, "edges")
    
    g_clu = components(g_all)
    g_edge_size=table(g_clu$csize)
    
    basic_linktype = c("H2", "C1H2", "H2O1","H1N1", "H3N1", "C16H30O1", "H1O3P1", "C5H8", "H2S1", 
                       "C1O2", "C2H4", "C6H10O5", "O1", "C2H2O1", "C6H4O1")
    merge_node_list_basic = merge_node_list
    merge_edge_list_basic = merge_edge_list %>%
      filter(linktype %in% basic_linktype) 
    # merge_node_list_basic = merge_node_list_basic %>% 
      # filter(name %in% merge_edge_list_basic$to|name %in% merge_edge_list_basic$from)
    
    g_basic = graph_from_data_frame(d = merge_edge_list_basic, vertices = merge_node_list_basic, directed = FALSE)
    g_basic_clu = components(g_basic)
    g_basic_edge_size=table(g_basic_clu$csize)
  }
  
  ## identify if there is large cluster ignored ####
  # {
  #   # Use all nodes and all correct connects
  #   g = g_basic
  #   
  #   clu = components(g)
  #   edge_size=table(clu$csize)
  #   
  #   
  #   mainnetwork = igraph::groups(clu)[table(clu$membership)>1000]
  #   subnetwork = igraph::groups(clu)[table(clu$membership)>10&table(clu$membership)<300]
  #   
  #   sub_nodelist = HMDB_clean2[subnetwork[[1]],]
  #   g_subnetwork = make_ego_graph(subnetwork[[1]][1], graph = g, order = diameter(g), mode = "all")[[1]]
  #   
  #   
  #   
  #   temp_vertex_color = igraph::as_data_frame(g, "vertices") %>%
  #     mutate(membership = clu$membership) %>%
  #     mutate(color = case_when(
  #       clu$csize[membership] == max(clu$csize) ~ my_palette[1],
  #       clu$csize[membership] != 1 ~ my_palette[2],
  #       clu$csize[membership] == 1 ~ my_palette[5]
  #     )) %>%
  #     pull(color)
  #   
  #   # pdf("full_all_connect_graph2.pdf")
  #   # plot.igraph(g,
  #   #      vertex.color = temp_vertex_color,
  #   #      vertex.label = "",
  #   #      vertex.frame.color = "#666666",
  #   #      # vertex.label.color = "black",
  #   #      # vertex.label.cex = 1,
  #   #      edge.color = 'grey',
  #   #      # edge.label = edge.attributes(g_subnetwork)$linktype,
  #   #      vertex.size = sqrt(degree(g, mode = c("all")))+0.5,
  #   #      # vertex.size = 2,
  #   #      # edge.arrow.size = .05,
  #   #      layout = layout_with_graphopt(g)
  #   #      # main = paste("Subnetwork of node", interested_node)
  #   # )
  #   # dev.off()
  # }
  
 
  ## Get results of
  # {
  #   # full rules
  #   g_edge_size
  #   # basic 15 rules
  #   g_basic_edge_size
  #   ## DQ 1910 metabolites
  #   ### biotransform_all 77 57 1776
  #   ### biotransform_basic15 137  75 1698
  #   ## Full 10304 metabolites 
  #   ### biotransform_all 318  324 9662
  #   ### biotransform_basic15 615  934 8755
  #   temp_size = g_edge_size
  #   unconnected = as.numeric(temp_size["1"])
  #   main_network = as.numeric(names(temp_size)[length(temp_size)])
  #   rest_network = nrow(merge_node_list) -unconnected-main_network
  #   print(c(unconnected, rest_network, main_network))
  # }
  
 
  ## ggplots
  connection_summary = data.frame(`Detected` = c(1698, 75, 137), 
                                  `All` = c(8755, 934, 615), 
                                  category = c("Main network", "Subnetworks", "Unconnected")) %>%
    gather(key = "metabolites", value = "number", -category) %>%
    mutate(metabolites = gsub("\\.", " ", metabolites))
    
  
  # pdf("bar_HMDB_connectivity.pdf", 
  #     width = 4,
  #     height = 3)
  # ggplot(connection_summary, aes(y = number, x = reorder(metabolites, -number), fill = forcats::fct_rev(category), label = number)) + 
  #   geom_bar(stat = "identity", 
  #            position = "fill" # make percentage graph
  #            ) +
  #   geom_text(size = 3,position = position_fill(vjust = 0.5)) +
  #   labs(# title = "Connectivity of HMDB metabolite formulas",
  #        x = NULL,
  #        y = "Fraction") + 
  #   guides(fill = guide_legend(
  #     title = "Metabolite status",
  #     reverse = F
  #     )) + 
  #   scale_y_continuous(expand = c(0,0),
  #                      labels = scales::percent,
  #                      breaks = scales::pretty_breaks(n = 8)
  #                      ) +
  #   scale_x_discrete(limits = c("Detected", "All")) +
  #   scale_fill_manual(values = c("Main network" = my_palette[1],
  #                     "Subnetworks" = my_palette[2],
  #                     "Unconnected" = my_palette[5])) +
  #   theme_classic(base_size = 12 # edit font size for all non-data text
  #                 ) +
  #   theme(plot.title = element_text(hjust = 0.5))
  # dev.off()
}

# Figure 2 ####
## Assume sigma = 1ppm
{
  ## evaluate proposed formula 
  {
    all_formula =bind_rows(Mset$NodeSet_network)
    
    merge_HMDB_propagation = HMDB_clean2 %>%  
      # mutate(ID = 1:nrow(.)) %>%
      merge(all_formula, by = "ID", all=T )
    
    merge_HMDB_propagation_correct = merge_HMDB_propagation %>% filter(MF == formula) %>%
      filter(ID %in% correct_connects$node1|ID %in% correct_connects$node2)
    merge_HMDB_propagation_wrong = merge_HMDB_propagation %>% filter(!ID %in% merge_HMDB_propagation_correct$ID, !is.na(formula))
  }
  
  ## evaluate CPLEX results
  {
    unknown_formula_CPLEX_correct = unknown_formula_CPLEX %>%
      merge(merge_HMDB_propagation_correct) 
    
    
    correct_connects_cplex = correct_connects %>%
      filter(node1 %in% unknown_nodes$ID  & node2 %in% unknown_nodes$ID)
    
    edge_info_sum_distinct = edge_info_sum %>%
      distinct(edge_id, .keep_all=T)
    
    edge_info_sum_distinct_correct = edge_info_sum_distinct %>%
      merge(correct_connects) 
    
    edge_info_sum_cplex = edge_info_sum %>%
      filter(ILP_result == 1)
    
    edge_info_sum_cplex_correct = edge_info_sum_cplex %>%
      merge(correct_connects) 
  }
  
  ## Main network with correct formula and correct linkage 
  {
    ## extended 51 connections
    mainnetwork = igraph::groups(g_clu)[table(g_clu$membership) == max(table(g_clu$membership))]
    g_main = make_ego_graph(g_all, order=diameter(g_all), nodes = mainnetwork[[1]][1], mode = c("all"))[[1]]
    # table(components(g_main)$csize)
    
    main_node_list = igraph::as_data_frame(g_main, "vertices") %>%
      rename(ID = name) %>%
      merge(merge_HMDB_propagation_correct, by = "ID",  all = T)
    
    main_edge_list = igraph::as_data_frame(g_main, "edges")
    g_main = graph_from_data_frame(d = main_edge_list, vertices = main_node_list, directed = FALSE)
    
    all_main_connects = all_connects %>%
      filter(node1 %in% main_node_list$name & node2 %in% main_node_list$name )
    ## Basic 15 
    
    mainnetwork = igraph::groups(g_basic_clu)[table(g_basic_clu$membership) == max(table(g_basic_clu$membership))]
    g_main_basic = make_ego_graph(g_basic, order=diameter(g_basic), nodes = mainnetwork[[1]][1], mode = c("all"))[[1]]
    # table(components(g_main_basic)$csize)
    
    main_node_list_basic = igraph::as_data_frame(g_main_basic, "vertices") %>%
      rename(ID = name) %>%
      merge(merge_HMDB_propagation_correct, by = "ID",  all = T)
    main_edge_list_basic = igraph::as_data_frame(g_main_basic, "edges")
    
    all_main_basic_connects = all_connects %>%
      filter(node1 %in% main_node_list_basic$name & node2 %in% main_node_list_basic$name )
  }
  
  
  ## Plot HMDB connection and formulas bar graph summary ####
  # {
  #   formula_summary = data.frame(`Optimized` = c(nrow(unknown_formula_CPLEX_correct), nrow(unknown_formula_CPLEX) - nrow(unknown_formula_CPLEX_correct)), 
  #                                `Pre-optimized` =  c(nrow(merge_HMDB_propagation_correct), nrow(unknown_formula) - nrow(merge_HMDB_propagation_correct)), 
  #                                category = c("Correct", "Incorrect")) %>%
  #     gather(key = "cohorts", value = "number", -category) %>%
  #     mutate(cohorts = gsub("\\.", "-", cohorts))
  #   connection_summary = data.frame(`Optimized` = c(nrow(edge_info_sum_cplex_correct), nrow(edge_info_sum_cplex)-nrow(edge_info_sum_cplex_correct)),
  #                                `Pre-optimized` =  c(nrow(edge_info_sum_distinct_correct), nrow(edge_info_sum_distinct)-nrow(edge_info_sum_distinct_correct)),
  #                                category = c("Correct", "Incorrect")) %>%
  #     gather(key = "cohorts", value = "number", -category) %>%
  #     mutate(cohorts = gsub("\\.", "-", cohorts))
  #   
  #   pdf("bar_HMDB_connections.pdf",
  #       width = 4,
  #       height = 4)
  #   # dev.new(width = 4, height = 3, unit = "in")
  #   ggplot(connection_summary, aes(y = number, x = reorder(cohorts, -number), fill = forcats::fct_rev(category), label = number)) +
  #     geom_bar(stat = "identity",
  #              position = "fill" # make percentage graph
  #              ) +
  #     geom_text_repel(size = 3,position = position_fill(vjust = 0.5), max.iter=1,arrow=T) +
  #     labs(x = NULL,
  #          title = "Connection assignment",
  #          y = "Fraction") +
  #     guides(fill = guide_legend(
  #       title = NULL,
  #       reverse = F
  #       )) +
  #     scale_y_continuous(expand = c(0,0),
  #                        labels = scales::percent,
  #                        breaks = scales::pretty_breaks(n = 8)
  #                        ) +
  #     expand_limits(y = 1.05) +
  #     # scale_x_discrete(limits = c("Detected", "All")) +
  #     scale_fill_manual(values = c("Correct" = my_palette[1],
  #                                  "Incorrect" = my_palette[5])) +
  #     theme_classic(base_size = 11 # edit font size for all non-data text
  #                   ) +
  #     theme(plot.title = element_text(hjust = 0.5),
  #           plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
  #           axis.text.x = element_text(angle = 45, hjust = 1))
  #   dev.off()
  # }
  
  # evaluate brute force formula results ####
  {
    Brute_force_formula = HMDB_clean2 %>%
      dplyr::select(ID, Exact_Mass, medMz, MF) %>%
      filter(ID %in% unknown_formula_CPLEX$ID)

    # Max number of element
    element_range = c()
    for(i in c("C","H","N","O","P","S")){
      element_range = c(element_range, max(sapply(Brute_force_formula$MF, elem_num_query, i)))
    }

    # rdbes = sapply(Brute_force_formula$MF, formula_rdbe)
    # max(rdbes)


    Elem_ratio_rule_formula_ls = list()
    for(i in 1:nrow(Brute_force_formula)){
      Elem_ratio_rule_formula_ls[[Brute_force_formula$ID[i]]] = mz_formula(Brute_force_formula$medMz[i],
                                                                           C_range = 0:element_range[1],
                                                                           H_range = 0:element_range[2],
                                                                           N_range = 0:element_range[3],
                                                                           O_range = 0:element_range[4],
                                                                           P_range = 0:element_range[5],
                                                                           S_range = 0:element_range[6],
                                                                           db_max = 99,
                                                                           charge=0,
                                                                           Elem_ratio_rule = T,
                                                                           ppm = 5)
    }

    Brute_force_formula["elem_ratio_rule_position"] = apply(Brute_force_formula, 1, function(x){
      temp_ID = as.numeric(x["ID"])
      temp_formula = x["MF"]
      temp_df = Elem_ratio_rule_formula_ls[[temp_ID]]
      if(!is.data.frame(temp_df)){return(-2)} # means no matched formulas
      if(!temp_formula %in% temp_df$formula){return(-1)} # means no matched formula
      which(temp_df$formula==temp_formula)
    })

    test = Brute_force_formula %>%
      filter(elem_ratio_rule_position < 0)
    table(Brute_force_formula["elem_ratio_rule_position"] )

    brute_force_summary = as.data.frame(table(Brute_force_formula$elem_ratio_rule_position), stringsAsFactors = F) %>%
      mutate(Var1 = as.numeric(Var1))

    brute_force_top1 = brute_force_summary %>%
      filter(Var1==1) %>%
      dplyr::select(Freq) %>% sum()

    brute_force_top3 = brute_force_summary %>%
      filter(Var1<=3, Var1>=1) %>%
      dplyr::select(Freq) %>% sum()
    
    brute_force_top10 = brute_force_summary %>%
      filter(Var1<=10, Var1>=1) %>%
      dplyr::select(Freq) %>% sum()
  }
}

result[[length(result)+1]] = list(sigma = rep(sigma,4),
                                  category = c("optimized_correct",
                                               "optimized_all",
                                               "propagated_correct",
                                               "propagated_all"),
                                  formula = c(nrow(unknown_formula_CPLEX_correct), 
                                              nrow(unknown_formula_CPLEX), 
                                              nrow(merge_HMDB_propagation_correct), 
                                              nrow(unknown_formula)),
                                  connection = c(nrow(edge_info_sum_cplex_correct), 
                                                 nrow(edge_info_sum_cplex), 
                                                 nrow(edge_info_sum_distinct_correct), 
                                                 nrow(edge_info_sum_distinct)),
                                  brute_force_category = c("top1",
                                                           "top3",
                                                           "top10",
                                                           "all"),
                                  brute_force = c(brute_force_top1,
                                                  brute_force_top3,
                                                  brute_force_top10,
                                                  nrow(unknown_formula_CPLEX))
                                  )
}
                                  
result_summary = lapply(result, bind_cols)


printtime = Sys.time()
timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
sink(paste(timestamp,"log.txt"))
print(paste0("random_seed",random_seed))
print(bind_cols(result))
sink()


## Assignment accuracy ####
{
  ppm_errors = sapply(result_summary, function(x) return(x$sigma[1]) )
  optimized_correct = sapply(result_summary, function(x) return(x$formula[1]) )
  optimized_all = sapply(result_summary, function(x) return(x$formula[2]) )
  top1 = sapply(result_summary, function(x) return(x$brute_force[1]) )
  top3 = sapply(result_summary, function(x) return(x$brute_force[2]) )
  
  
  colfunc = colorRampPalette(c(my_palette[3], "white"))
  colfunc(5)[1]
  
  HMDB_ppm_summary = data.frame(ppm_errors = ppm_errors,
                                NetID = optimized_correct/optimized_all,
                                `Heuristic top1` = top1 / optimized_all,
                                `Heuristic top3` = top3 / optimized_all) %>%
    gather(key = "cohorts", value = "number", -ppm_errors) %>%
    mutate(cohorts = gsub("\\.", " ", cohorts)) %>%
    mutate(color = case_when(
      cohorts == "NetID" ~ "red",
      cohorts == "Heuristic top1" ~ "black",
      cohorts == "Heuristic top3" ~ "grey"
    ))
    
  num_x = 4
  
  pdf("bar_error_summary.pdf",
      width = 4,
      height = 4)
  # dev.new(width = 1, height = 1, unit = "in")
  ggplot(HMDB_ppm_summary, aes(y = number, x = factor(ppm_errors), group = forcats::fct_rev(cohorts), color = cohorts)) +
    geom_line(stat = "identity",
              size = 1.5
              # linetype = rep(c("solid", rep("dashed", 2)),2),
              # color = rep(c("red", "black", "grey"), num_x)
              ) + 
    geom_point(shape = 16,
               size = 4) + 
    labs(x = "Gaussian noise level (ppm)",
         title = "Formula assignment accuracy",
         y = "Fraction") +
    guides(color = guide_legend(
      title = NULL,
      reverse = T
      )) +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,1.05),
                       breaks = scales::pretty_breaks(n = 5)
                       ) +
    # scale_x_discrete(limits = c("0.5", "1")) + 
    # scale_colour_manual(values=c("red", "black", "grey")) +
    scale_colour_manual(values=c("NetID" = my_palette[1],"Heuristic top1" = my_palette[5], "Heuristic top3" = my_palette[3])) +
    theme_classic(base_size = 11 # edit font size for all non-data text
                  ) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))
  
  dev.off()
  
}


## Plot step network ####
{
  # Plot step network
  {
    g = g_main
    
    node_list = main_node_list
    # node_list = unknown_formula
    edge_list = main_edge_list
    
    step_info = as.data.frame(table(node_list$steps), stringsAsFactors = F) %>%
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
    
    # pdf("step_graph_equalsize.pdf")
    # # for(i in 1:1){
    #   plot.igraph(g_main,
    #               vertex.color = vertex_color,
    #               vertex.label = "",
    #               vertex.frame.color = "black",
    #               # vertex.label.color = "black",
    #               # vertex.label.cex = 1,
    #               edge.color = 'grey',
    #               # edge.label = edge.attributes(g_subnetwork)$linktype,
    #               # vertex.size = sqrt(degree(g_main, mode = c("all")))+0.5,
    #               vertex.size = 2,
    #               # edge.arrow.size = .05,
    #               edge.size = 0.001,
    #               layout = layout_with_graphopt(g_main)
    #               # main = paste("Subnetwork of node", interested_node)
    #   )
    # # }
    # dev.off()
    pdf("step_bar_graph.pdf", 
        width = 4,
        height = 3)
    step_info2$Step = as.factor(step_info2$Step)
    step_info2$Step  = factor(step_info2$Step, levels=c(levels(step_info2$Step)[-1],levels(step_info2$Step)[1]))
    ggplot(step_info2, aes(x = Step)) +
      geom_bar(aes(y = No.metabolites), stat = 'identity', fill = color_spectral[1:nrow(step_info2)]) + 
      # scale_x_discrete(limits = temp_position) +
      geom_point(aes(y = Total.metabolites/4)) + 
      geom_line(aes(y = Total.metabolites/4), group=1) +
      geom_text(aes(y = Total.metabolites/4, label = Total.metabolites), vjust = -.5) +
      scale_y_continuous(# sec.axis = sec_axis(~.*4, name = "Total.metabolites"), 
                         limits = c(0,500), expand = c(0,0)) +
      # theme(legend.title = element_blank()) + 
      theme(legend.position = "right") +
      theme_classic()
    dev.off()
  }
}


## Evaluate choice of biotransform rules ####
{ 
  # evaluate from basic 1 biotransform, adding which one has most gain
  {
    network_size_linktype = list()
    i = unique(EdgeSet$Merge$linktype)[1]
    for(i in unique(EdgeSet$Merge$linktype)){
      selected_linktype = c(basic_linktype, i)
      node_list = main_node_list
      edge_list = main_edge_list %>%
        filter(linktype %in% selected_linktype) 
      node_list = node_list %>% 
        filter(name %in% edge_list$to|name %in% edge_list$from)
      # table(merge_edge_list$linktype)
      #Generate network graph using node_list and edge_list above
      
      g <- graph_from_data_frame(d = edge_list, vertices = node_list, directed = FALSE)
      clu = components(g)
      edge_size=max(as.numeric(names(table(clu$csize))))
      network_size_linktype[[length(network_size_linktype)+1]] = data.frame(linktype=i, max_size = edge_size, stringsAsFactors = F)
    }
    
    summary_network_gain = bind_rows(network_size_linktype) %>%
      mutate(gain = max_size - min(max_size))
    
  }

  # evaluate if there is biotransform that is not importnat
  {
    node_list = main_node_list_basic
    edge_list = main_edge_list_basic
    
    edge_size= edge_membership = list()
    for (i in unique(edge_list$linktype)){
      temp_edge_list=edge_list[edge_list$linktype!=i,]
      temp_merge_node_list = main_node_list %>% filter(name %in% temp_edge_list$from|name %in% temp_edge_list$to)
      temp_g <- graph_from_data_frame(d = temp_edge_list, vertices = temp_merge_node_list, directed = FALSE)
      temp_clu = components(temp_g)
      edge_size[[i]]=table(temp_clu$csize)
      edge_membership[[i]]=temp_clu$membership
    }
    
    drop_out_size=sapply(edge_size, function(x) max(as.numeric(names(x))))
    
    number_link = data.frame(table(edge_list$linktype), stringsAsFactors = F)
    drop_out_size2 = data.frame(Var1 = names(drop_out_size), drop_out_size=drop_out_size, stringsAsFactors = F) %>%
      mutate(drop_out_num = nrow(node_list) - drop_out_size)
    
    summary_linktype = merge(number_link,drop_out_size2)
  }
}




## Evaluate different scoring ####
results = list()
for(edge_bonus_test in seq(1,2.5,.5)){
 

# evaluate CPLEX results
{
  # cplex formulas
  unknown_nodes = CPLEXset$formula$unknown_nodes[,1:3]
  unknown_formula = CPLEXset$formula$unknown_formula %>%
    mutate(ILP_result = CPLEX_x[1:nrow(.)])
  unknown_formula_CPLEX = unknown_formula %>%
    filter(ILP_result!=0)
  unknown_node_CPLEX = merge(unknown_nodes,unknown_formula_CPLEX,by.x = "ID", by.y = "ID",all=T)
  
  cplex_formula = unknown_formula %>% filter(ILP_result !=0)
  cplex_formula_correct = cplex_formula %>% 
    inner_join(main_node_list[,c("ID", "MF")], by= c("ID"="ID", "formula"="MF"))
  # table(cplex_formula_correct$ILP_result)
  cplex_formula_wrong = cplex_formula %>% 
    anti_join(main_node_list[,c("ID", "MF")], by= c("ID"="ID", "formula"="MF"))
  
  # cplex connects
  edge_info_sum = edge_info_sum %>%
    mutate(ILP_result = CPLEX_x[(nrow(unknown_formula)+1):length(CPLEX_x)]) 
  cplex_connects = edge_info_sum %>%
    filter(ILP_id1 <= nrow(unknown_formula), ILP_id2 <= nrow(unknown_formula)) %>%
    mutate(ID1 = unknown_formula$ID[ILP_id1], ID2 = unknown_formula$ID[ILP_id2]) %>%
    filter(ILP_result==1)
  
  cplex_connects_correct = cplex_connects %>% filter(edge_id %in% main_edge_list$edge_id)
  
  cplex_connects_wrong = cplex_connects %>% filter(!edge_id %in% main_edge_list$edge_id)
}

  
  results[[length(results)+1 ]] = list(cplex_formula_correct, cplex_formula_wrong, all_formula)
} 
results2[[length(results2)+1]] = results
# }
lapply(results2, function(results){
  sapply(results, function(x){c(nrow(x[[1]]),nrow(x[[2]]), nrow(x[[3]]))})
})



