# Main Codes ####
## Read files ####

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
printtime = Sys.time()
timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
source("NetID_function.R")

work_dir = "Unknown in IO RT"
work_dir = "matt"
# filename_wl = "pks_liver_pos_buff_2020-02-21.xlsx"
# filename_wl = "pks_Io_pos__2020-06-12.xlsx"
# filename_wl = "pks_Io_neg__2020-06-11.xlsx"
# filename_wl = "pks_Rt_pos__2020-06-12.xlsx"
# filename_wl = "Lu-Table-S4-final_Sc_neg.xlsx"
# filename_wl = "pks_Rt_pos__2020-06-12.xlsx"
# work_dir = "WL_liver_pos"
# filename_wl = "pks_liver_pos_buff_2020-02-21.xlsx"
ion_mode = 1
MS2_library_file = "./dependent/HMDB_pred_MS2_pos.rds"
MS2_folder = "MS2_pos_200524"
LC_method = "lipids" # "Hilic_Rutgers_QEPlus" "Hilic_25min_QE", lipids is empty

MS2_library = readRDS(MS2_library_file)
mass_dist_gamma_rate = 1 # smaller means more penalty on mass error, similar to sd
Empirical_rules_file = "./dependent/Empirical_rules.csv"
propagation_rule = read.csv("./dependent/propagation_rule.csv", row.names = 1)

# sink(paste(timestamp,"log.txt"))
# sink()

{
  Mset = list()
  # Mset[["Library"]] = read.csv("./dependent/HMDB_CHNOPS_clean.csv", stringsAsFactors = F)
  Mset[["Library_HMDB"]] = read.csv("./dependent/hmdb_library.csv", stringsAsFactors = F)
  
  Mset[["Library_known"]] = read.csv("./dependent/known_library.csv", stringsAsFactors = F) %>%
    filter(!is.na(.[,eval(LC_method)]))
  
  Mset[["Empirical_rules"]]=Read_rule_table(rule_table_file = Empirical_rules_file)
  
  Mset[["Global_parameter"]]=list(mode = ion_mode,
                                  LC_method = LC_method,
                                  propagation_rule = propagation_rule)
}

## WL data format ####
{
  setwd(work_dir)

  filename = "raw_data.csv"
  Mset[["Raw_data"]] <- read_raw_data(filename)

  # WL = readxl::read_xlsx(filename_wl) %>%
  # # WL = readxl::read_xlsx(filename_wl,
  # #                        sheet = "Sheet1",
  # #                        guess_max = 1e6
  # #                        ) %>%
  #   dplyr::rename(medRt = rt,
  #                 medMz = mz) %>%
  #   dplyr::rename(Formula = Formula...32,
  #                 Feature = Feature...33,
  #                 Background = Background...25) %>%
  #   filter(T)
  # 
  # 
  # raw_data_WL = Mset$Raw_data %>%
  #   merge(WL, all = T) %>%
  #   # mutate(id = 1:nrow(.),
  #   #        liver = 10^sig,
  #   #        liver2 = 10^sig) %>%
  #   mutate(id = 1:nrow(.),
  #          yeast = 10^sig,
  #          yeast2 = 10^sig) %>%
  #   filter(is.na(Background)) %>%
  #   dplyr::select(colnames(Mset$Raw_data))
  # Mset[["Raw_data"]] = raw_data_WL
}

## Initialise ####
{
  manual_library_file = "manual_library.csv"
  Mset[["Manual_library"]] = read_manual_library(manual_library_file)

  Mset[["Cohort"]]=Cohort_Info(Mset, first_sample_col_num = 15)
  print(Mset$Cohort)
  
  #Clean-up duplicate peaks 
  Mset[["Data"]] = Peak_cleanup(Mset,
                                mz_tol=2/10^6, 
                                rt_tol=0.1,
                                inten_cutoff=0e4,
                                high_blank_cutoff = 0,
                                first_sample_col_num = 15)
  
  print(c(nrow(Mset$Raw_data), nrow(Mset$Data)))
  
  manual_library_file = "manual_library.csv"
  Mset[["Manual_library"]] = read_manual_library(manual_library_file)
}



## NodeSet and FormulaSet ####
{
  NodeSet = Initiate_nodeset(Mset)
  NodeSet = Add_MS2_nodeset(MS2_folder = MS2_folder, 
                            NodeSet)
  
  test = sapply(NodeSet,"[[", "MS2")
  sum(sapply(test, is.null))
  
  LibrarySet = Initiate_libraryset(Mset)
  
  FormulaSet = Initilize_empty_formulaset(NodeSet)
  
  FormulaSet = Match_library_formulaset(FormulaSet, 
                                        Mset, NodeSet, 
                                        LibrarySet,
                                        expand = F,
                                        ppm_tol = 5e-6)
  
  Sys_msr_error = Check_sys_error(NodeSet, FormulaSet, LibrarySet, 
                                  RT_match = F)
  
  mass_adjustment = abs(Sys_msr_error$ppm_adjust + Sys_msr_error$abs_adjust/250*1e6) > 0.2
  # mass_adjustment = Sys_msr_error$fitdistData$estimate["mean"] > 10*Sys_msr_error$fitdistData$sd["mean"]
  
  if(mass_adjustment){
    NodeSet = lapply(NodeSet, function(Node){
      Node$mz = Node$mz * (Sys_msr_error$ppm_adjust * 1e-6 + 1) + Sys_msr_error$abs_adjust
      return(Node)
    })
    FormulaSet = Initilize_empty_formulaset(NodeSet)
    FormulaSet = Match_library_formulaset(FormulaSet, 
                                          Mset, NodeSet, 
                                          LibrarySet,
                                          expand = F,
                                          ppm_tol = 5e-6)
  }
  
}

## EdgeSet ####
{
  EdgeSet = Initiate_edgeset(Mset, NodeSet, 
                             mz_tol_abs = 0, mz_tol_ppm = 10, 
                             rt_tol_bio = Inf, rt_tol_nonbio = 0.2)
  
  # Extension edgeset
  peak_group = Peak_grouping(NodeSet, RT_cutoff = 0.2, inten_cutoff = 1e4)
  EdgeSet_ring_artifact = Ring_artifact_connection(peak_group, ppm_range_lb = 50, ppm_range_ub = 1000, ring_fold = 50, inten_threshold = 1e6)
  EdgeSet_oligomer_multicharge = Oligomer_multicharge_connection(peak_group, ppm_tol = 10)
  EdgeSet_multicharge_isotope = Multicharge_isotope_connection(EdgeSet, EdgeSet_oligomer_multicharge)
  EdgeSet_heterodimer = Heterodimer_connection(peak_group, NodeSet, ppm_tol = 10, inten_threshold = 1e6)
  EdgeSet_experiment_MS2_fragment = Experiment_MS2_fragment_connection(peak_group, NodeSet, ppm_tol = 10, inten_threshold = 1e5)
  EdgeSet_library_MS2_fragment_df = Library_MS2_fragment_connection(peak_group, FormulaSet, MS2_library,
                                                                    inten_threshold = 1e5,
                                                                    ppm_tol = 10, abs_tol = 1e-4)
  
  EdgeSet_all_df = merge_edgeset(EdgeSet, 
                                 EdgeSet_ring_artifact, 
                                 EdgeSet_oligomer_multicharge, 
                                 EdgeSet_multicharge_isotope,
                                 EdgeSet_heterodimer,
                                 EdgeSet_experiment_MS2_fragment,
                                 EdgeSet_library_MS2_fragment_df)
  
  
}

## Candidate formula pool ####
{
  FormulaSet = Propagate_formulaset(Mset, 
                                    NodeSet,
                                    FormulaSet,
                                    EdgeSet_all_df,
                                    biotransform_step = 2,
                                    artifact_step = 3,
                                    propagation_ppm_threshold = 5e-6,
                                    propagation_abs_threshold = 2e-4,
                                    record_RT_tol = 0.15,
                                    record_ppm_tol = 5e-6)
  
  all_formulas = bind_rows(FormulaSet) %>%
    distinct(node_id, formula, category, .keep_all = T) %>%
    filter(T)
  
}

## CplexSet & Scoring ####
{
  CplexSet = list()
  FormulaSet_df = Score_formulaset(FormulaSet,
                                   database_match = 0.5,
                                   MS2_match = 1,
                                   MS2_match_cutoff = 0.8,
                                   MS2_similarity = 0.5,
                                   MS2_similarity_cutoff = 0.5,
                                   manual_match = 0.5,
                                   rt_match = 1, 
                                   known_rt_tol = 0.5,
                                   bio_decay = -1,
                                   artifact_decay = -0.5)
  # Initialize
  CplexSet[["ilp_nodes"]] = initiate_ilp_nodes(FormulaSet_df, artifact_decay = 1) 
  CplexSet[["ilp_edges"]] = initiate_ilp_edges(EdgeSet_all_df, CplexSet, 
                                               Exclude = "")
  CplexSet[["heterodimer_ilp_edges"]] = initiate_heterodimer_ilp_edges(EdgeSet_all_df, CplexSet, NodeSet)
  
  # Score
  CplexSet[["ilp_nodes"]] = score_ilp_nodes(CplexSet, mass_dist_gamma_rate = mass_dist_gamma_rate, 
                                            formula_score = 0, metabolite_score = 0, 
                                            isotope_Cl_penalty = -1, unassigned_penalty = -0.5)
  CplexSet[["ilp_edges"]] = score_ilp_edges(CplexSet, NodeSet,
                                            rule_score_biotransform = 0, rule_score_adduct = 0.5, 
                                            rule_score_oligomer = 0.5, rule_score_natural_abundance = 1,
                                            rule_score_fragment = 0.3, rule_score_ring_artifact = 2,
                                            rule_score_radical = 0.2, rule_score_multicharge = 0.5, 
                                            rule_score_multicharge_isotope = 1, 
                                            rule_score_experiment_MS2_fragment = 1, rule_score_library_MS2_fragment = 0.3,
                                            rt_penalty_artifact_cutoff = 0.05, rt_penalty_artifact_ratio = 5,
                                            inten_score_isotope = 1, 
                                            MS2_score_similarity = 1, MS2_similarity_cutoff = 0.3,
                                            MS2_score_experiment_fragment = 0.5)
  
  CplexSet[["heterodimer_ilp_edges"]] = score_heterodimer_ilp_edges(CplexSet, rule_score_heterodimer = 0,
                                                                    MS2_score_experiment_fragment = 0.5)
  
  CplexSet[["para"]] = Initiate_cplexset(CplexSet)
  CplexSet[["para_lp"]] = Initiate_cplexset_lp(CplexSet)
}
  
save.image(paste0(timestamp,".RData"))
print(Sys.time()-printtime)
print(paste("Complexity is", CplexSet$para$nc, "variables and", CplexSet$para$nr, "constraints."))
print(paste("Complexity is", CplexSet$para_lp$nc, "variables and", CplexSet$para_lp$nr, "constraints."))
sapply(CplexSet$para_lp, length)

## Run_cplex ####
{
  CplexSet[["init_solution"]] = list(Run_cplex(CplexSet, obj_cplex = CplexSet$para$obj,
                                               relative_gap = 1e-3, total_run_time = 2500))
  # CplexSet[["lp_solution"]] = list(Run_cplex_lp(CplexSet, obj_cplex = CplexSet$para_lp$obj,
  #                                               total_run_time = 5000))
  # Test_para_CPLEX(CplexSet, obj_cplex = CplexSet$para$obj, 
  #                 para = c(0), para2 = NA, 
  #                 relative_gap = 1e-1, total_run_time = 3000)
}

## Read out CPLEX results ####
{
  CPLEX_all_x = Read_cplex_result(CplexSet$init_solution)
  # CPLEX_all_x = Read_cplex_result(CplexSet$lp_solution)
  # CPLEX_all_x = Read_CPLEX_result(CPLEXset$Pmt_solution)
  CPLEX_x = rowMeans(CPLEX_all_x,na.rm=T)
  
  ilp_nodes = CplexSet$ilp_nodes %>%
    mutate(ilp_result = CPLEX_x[1:nrow(.)]) %>%
    inner_join(Mset$Data %>% dplyr::select(id, Input_id, log10_inten, medMz, medRt), by = c("node_id"="id")) %>%
    filter(T)
  
  ilp_edges = CplexSet$ilp_edges %>%
    mutate(ilp_result = CPLEX_x[1:nrow(.) + nrow(ilp_nodes)]) %>%
    # filter(ilp_result == 1) %>%
    filter(T)
  
  heterodimer_ilp_edges = CplexSet$heterodimer_ilp_edges %>%
    mutate(ilp_result = CPLEX_x[1:nrow(.) + nrow(ilp_nodes) + nrow(ilp_edges)])
  
  ilp_nodes_ilp = ilp_nodes %>%
    filter(ilp_result > 1e-6) 

}

# # Add comparison annotation
# {
#   
#   CPLEX_all_x = Read_cplex_result(CplexSet$lp_solution)
#   # CPLEX_all_x = Read_CPLEX_result(CPLEXset$Pmt_solution)
#   CPLEX_x = rowMeans(CPLEX_all_x,na.rm=T)
#   
#   ilp_nodes = ilp_nodes %>%
#     mutate(lp_result = CPLEX_x[1:nrow(.)])
#   
#   ilp_edges = ilp_edges %>%
#     mutate(lp_result = CPLEX_x[1:nrow(.) + nrow(ilp_nodes)]) %>%
#     # filter(ilp_result == 1) %>%
#     filter(T)
#   
#   heterodimer_ilp_edges = heterodimer_ilp_edges %>%
#     mutate(lp_result = CPLEX_x[1:nrow(.) + nrow(ilp_nodes) + nrow(ilp_edges)])
#   
#   test = ilp_nodes %>%
#     arrange(-log10_inten, -lp_result) %>%
#     # filter(lp_result > 1e-6)
#     filter(ilp_result != 0 | lp_result != 0)
#     # filter(ilp_result != lp_result)
#   # hist(ilp_nodes$lp_result)
#   
#   
#   
#   # tabyl(test, lp_result)
# }

# Path annotation ####
{
  core_annotation = core_annotate(ilp_nodes, FormulaSet_df, LibrarySet)
  core_annotation_unique = core_annotate_unique(core_annotation)
  g_met = initiate_g_met(ilp_nodes, ilp_edges)
  
  met_dist_mat = initiate_met_dist_mat(g_met, ilp_nodes, core_annotation_unique)
  
  g_nonmet = initiate_g_nonmet(ilp_nodes, ilp_edges, heterodimer_ilp_edges)
  
  nonmet_dist_mat = initiate_nonmet_dist_mat(g_nonmet, ilp_nodes, only_ilp_result = T)
  
  core_annotation_nonmet = igraph::as_data_frame(g_nonmet, "vertices") %>%
    filter(category != "Unknown") %>%
    filter(steps %% 1 == 0) %>%
    mutate(core_annotate = paste("Peak", node_id, formula)) %>%
    dplyr::rename(ilp_node_id = name) %>%
    mutate(ilp_node_id = as.integer(ilp_node_id))
  
  ilp_edges_annotate_met = igraph::as_data_frame(g_met)
  ilp_edges_annotate_nonmet = igraph::as_data_frame(g_nonmet)
  
  # Roughly 10s ~ 1000 entries(including skipped ones)
  path_annotation = rep("", nrow(ilp_nodes))
  
  for(i in 1:nrow(ilp_nodes)){
    if(i %% 1000 == 0)print(i)
    # skip ilp_nodes not as 1
    if(ilp_nodes$ilp_result[i] < 0.01){next}
    
    temp_ilp_node_id = ilp_nodes$ilp_node_id[i]
    if(ilp_nodes$class[i] == "Artifact"){
      path_annotation[i] = track_annotation_nonmet(temp_ilp_node_id, ilp_edges_annotate = ilp_edges_annotate_nonmet,
                                                   g_nonmet, graph_path_mode = "out",
                                                   nonmet_dist_mat, core_annotation_nonmet)
    } else if(ilp_nodes$class[i] == "Metabolite"){
      path_annotation[i] = track_annotation_met(temp_ilp_node_id, ilp_edges_annotate = ilp_edges_annotate_met,
                                                g_met, graph_path_mode = "all",
                                                met_dist_mat, core_annotation_unique)
    } else {
      path_annotation[i] = "Unknown"
    }
  }
  
  ilp_nodes_annotation = ilp_nodes %>%
    mutate(path = path_annotation) %>%
    mutate(class = ifelse(class == "Metabolite" & (steps != 0 | transform != ""), "Putative Metabolite", class)) %>%
    # filter(ilp_result != 0) %>%
    filter(T)

  test = ilp_nodes_annotation %>%
    filter(ilp_result == 1) %>%
    filter(class == "Putative Metabolite")
  
  saveRDS(list(ilp_nodes = ilp_nodes_annotation,
               ilp_edges = ilp_edges,
               heterodimer_ilp_edges = heterodimer_ilp_edges,
               FormulaSet_df = FormulaSet_df,
               LibrarySet = LibrarySet),
          paste(timestamp, "output.rds", sep="_"))
}

print("total run time")
save.image()
print(Sys.time()-printtime)


## View Results ####
{
  test = ilp_nodes_annotation %>% 
    filter(ilp_result > 0.01)
  test2 = merge(WL, test, by.x = "Index", by.y = "Input_id", all.x = T, suffixes = c("",".y")) %>%
    dplyr::select(colnames(WL), node_id, class, path, formula, ppm_error, parent_id, parent_formula, transform, category, mass)
    
  test3 = test2 %>%
    # filter(is.na(Background)) %>%
    # filter(is.na(Background)) %>%
    mutate(Formula = ifelse(is.na(Formula), "", Formula)) %>%
    mutate(Formula = check_chemform(isotopes, Formula)$new_formula)
  

  
  test3_filter = test3 %>%
    # filter(grepl("Si", path)) %>%
    # mutate(C_count = sapply(str_extract_all(formula, "(?<=C)([:digit:]|\\.)+"), 
    #                         function(x){sum(as.numeric(x))}),
    #        N_count = sapply(str_extract_all(formula, "(?<=N)([:digit:]|\\.)+"), 
    #                         function(x){sum(as.numeric(x))})) %>%
    # mutate(CN_match = (N_count == N) & (C_count == C)) %>%
    filter(T)
  
  # write_csv(test3_filter, "liver_pos.csv", na="")
  
  
  
  library(janitor)
  crosstable = tabyl(test3, Feature, category)
  crosstable2 = tabyl(test3, category, Feature)
  
  test5 = test3 %>%
    mutate(Feature2 = case_when(
      is.na(Feature) ~ "Unknown",
      Feature == "Metabolite" ~ "Metabolite",
      T ~ "Artifact"
      )) %>%
    # filter(Feature == "Metabolite", category == "Metabolite") %>%
    # filter(Formula == formula) %>%
    filter(T)
  crosstable3 = tabyl(test5, class, Feature2)
  cross_netid = tabyl(test5, class)
  cross_bmw = tabyl(test5, Feature2)
  
  test6 = test3 %>%
    filter(Feature != "Metabolite" | is.na(Feature)) %>%
    filter(category != "Metabolite" ) %>%
    filter(Feature...23 == "Metabolite") %>%
    # mutate(C_count = sapply(formula, elem_num_query, "C"),
    #        N_count = sapply(formula, elem_num_query, "N")) %>%
    # mutate(C_count_BMW = sapply(Formula, elem_num_query, "C"),
    # N_count_BMW = sapply(Formula, elem_num_query, "N")) %>%
    # mutate(CN_match = (N_count == N) & (C_count == C)) %>%
    # mutate(CN_match_BMW = (N_count_BMW == N) & (C_count_BMW == C)) %>%
    # filter(Formula == formula) %>%
    filter(T)
  
  test_7 = tabyl(test3, Feature...23)
  
}

## Query specific nodes or edges ####
{
  Input_id = 603
  Mset$Data$id[Mset$Data$Input_id == Input_id]

  node_id_selected = 7963
  ilp_nodes_selected = ilp_nodes %>%
    filter(node_id %in% node_id_selected)
  ilp_edges_node_related = ilp_edges %>%
    filter(node1 == node_id_selected | node2 == node_id_selected)%>%
    # filter(ilp_nodes1 %in% c(11486, 11488) | ilp_nodes2 %in% c(11486, 11488)) %>%
    arrange(-cplex_score, category) %>%
    filter(ilp_result != 0 | lp_result !=0 ) %>%
    filter(T)
  heterodimer_ilp_edges_node_related = heterodimer_ilp_edges %>%
    filter(node1 == node_id_selected | node2 == node_id_selected | linktype == node_id_selected) %>%
    arrange(-cplex_score, category, edge_id) %>%
    filter(T)
  
  sum(heterodimer_ilp_edges_node_related$lp_result)
  
  edges_selected = EdgeSet_all_df %>%
    filter(node1 == node_id_selected | node2 == node_id_selected)
  
  edge_related_to_node = EdgeSet_all_df %>%
    filter(node1 == node_id_selected | node2 == node_id_selected)
  ilp_nodes_related_to_edge = ilp_nodes_ilp %>%
    filter(node_id %in% c(edge_related_to_node$node1, edge_related_to_node$node2)) 
  
  test_edge = ilp_edges %>%
    filter(category == "Experiment")
  tabyl(ilp_nodes_ilp$category)
  
}


# Evaluate parameter ####
# Quite useful : )
{
  # How RT differ between connected peaks
  node_rt = sapply(NodeSet, "[[", "RT")
  test_rt_correct = ilp_nodes_ilp %>%
    filter(category == "Natural_abundance") %>%
    # filter(Feature == "Isotope") %>%
    mutate(rt_parent = node_rt[parent_id],
           rt_node = node_rt[node_id],
           rt_dif = rt_parent - rt_node) %>%
    filter(!is.na(rt_dif)) %>%
    # filter(quantile(rt_dif, 0.1) < rt_dif,
    #        quantile(rt_dif, 0.9) > rt_dif) %>%
    # filter(abs(rt_dif) > 0.05) %>%
    filter(T)
  fitdist_rt_correct = fitdistrplus::fitdist(test_rt_correct$rt_dif, "norm")
  plot(fitdist_rt_correct)
  print(fitdist_rt_correct)
  shapiro.test(test_rt_correct$rt_dif)
  
  # Test how mz differ 
  node_mz = sapply(NodeSet, "[[", "mz")
  test_mz_correct = ilp_nodes_ilp %>%
    filter(category == "Adduct") %>%
    # filter(Feature != "Adduct") %>%
    mutate(mz_parent = node_mz[parent_id],
           mz_node = node_mz[node_id],
           mz_dif_edge = mz_node - mz_parent - formula_mz(transform),
           mz_node_dif = formula_mz(formula) - mz_node,
           mz_parent_dif = formula_mz(parent_formula) - mz_parent,
           ppm_mz_dif_node = mz_node_dif/mz_node * 1e6,
           ppm_mz_dif_parent = mz_parent_dif/mz_parent * 1e6,
           ppm_mz_dif_edge = mz_dif_edge/mz_node * 1e6) %>%
    filter(!is.na(mz_dif_edge)) %>%
    filter(abs(ppm_mz_dif_edge) > 1) %>%
    filter(T)
  fitdist_mz_correct = fitdistrplus::fitdist(test_mz_correct$ppm_mz_dif_edge, "norm")
  plot(fitdist_mz_correct)
  print(fitdist_mz_correct)
}
# ## Maven output format ####
# {
# 
#   test = ilp_nodes_annotation %>%
#     filter(ilp_result > 0.01) %>%
#     # filter(class == "Metabolite", steps != 0) %>%
#     filter(T)
# 
#   test2 = Mset$Data %>%
#     inner_join(test) %>%
#     dplyr::select(colnames(Mset$Data), node_id, class, formula, path, ppm_error, -id) %>%
#     arrange(Input_id) %>%
#     filter(class == "Unknown")
# 
#   
#   test3 = read_csv("Mcap_NetID.csv") %>%
#     arrange(Input_id) %>%
#     dplyr::select(Input_id, class, formula, path) %>%
#     dplyr::rename(Class = class,
#                   Formula = formula,
#                   Path = path)
#   
#   test4 = test2 %>%
#     inner_join(test3) %>%
#     # filter(log10_inten > 5) %>%
#     # filter(grepl("Br", formula))
#     filter(Path == path)
#   
#   library(janitor)
#   tabyl(test4, Class, class)
#     
#   # write_csv(test2, paste(work_dir, timestamp,"NetID.csv", sep="_"))
# }



## Merge with ANOVA results ####
# {
#   temp = ilp_nodes %>%
#     filter(Input_id == 1)
#   
#   HMDB_result = read.csv("mdata.csv") %>%
#     inner_join(ilp_nodes %>%
#                  filter(ilp_result != 0), 
#                by = c("ID"="Input_id"))
#   
#   write.csv(HMDB_result, "merge.csv")
#   
#   # save.image()
# }

## Get publishable numbers ####
{
  # Seed step
  {
    seed_nodes = ilp_nodes %>%
      filter(steps == 0) %>%
      filter(class != "Unknown")
    table(seed_nodes$class)
    nrow(seed_nodes)
    # number of nodes that are seeds
    length(table(seed_nodes$node_id))
    
  }
  # Propagation step
  {
    ilp_nodes %>% filter(class != "Unknown") %>% nrow()
    length(table(ilp_nodes %>% filter(class != "Unknown") %>% pull(node_id)))
  }
  
  
}
# Deprecated ####
