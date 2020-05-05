# Main Codes ####
## Read files ####

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
printtime = Sys.time()
timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
source("NetID_function.R")

work_dir = "Mcap"
MS2_folder = "MS2"
MS2_library_file = "./dependent/HMDB_pred_MS2_neg3.rds"
MS2_library = readRDS(MS2_library_file)

ion_mode = 1
mass_dist_gamma_rate = 1 # smaller means more penalty on mass error, similar to sd
# LC_method = "Hilic_25min_QE"
LC_method = "Hilic_Rutgers_QEPlus"
Empirical_rules_file = "./dependent/Empirical_rules.csv"
print(work_dir)
print(ion_mode)

# sink(paste(timestamp,"log.txt"))
# sink()

{
  Mset = list()
  # Mset[["Library"]] = read.csv("./dependent/HMDB_CHNOPS_clean.csv", stringsAsFactors = F)
  Mset[["Library_HMDB"]] = read.csv("./dependent/hmdb_formula_order.csv", stringsAsFactors = F)
  Mset[["Library_known"]] = read.csv("./dependent/known_library.csv", stringsAsFactors = F) %>%
    filter(!is.na(.[,eval(LC_method)]))
  
  Mset[["Empirical_rules"]]=Read_rule_table(rule_table_file = Empirical_rules_file)
  
  Mset[["Global_parameter"]]=list(mode = ion_mode,
                                    LC_method = LC_method)
}

## WL data format ####
{
  setwd(work_dir)
  
  filename = "raw_data.csv"
  Mset[["Raw_data"]] <- read_raw_data(filename)
  
  # filename_wl = "Lu-Table-S4-final.xlsx"
  # WL = readxl::read_xlsx(filename_wl) %>%
  #   dplyr::rename(medRt = rt,
  #                 medMz = mz) %>%
  #   dplyr::rename(Formula = Formula...34,
  #                 Feature = Feature...35,
  #                 Background = Background...27) 
  # 
  # raw_data_WL = Mset$Raw_data %>%
  #   merge(WL, all = T) %>%
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
                                mz_tol=1/10^6, 
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
 
  LibrarySet = Initiate_libraryset(Mset)
  
  FormulaSet = Initilize_empty_formulaset(NodeSet)
  
  FormulaSet = Match_library_formulaset(FormulaSet, 
                                        Mset, NodeSet, 
                                        LibrarySet,
                                        expand = T,
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
                                          expand = T,
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
  EdgeSet_oligomer = Oligomer_connection(peak_group, ppm_tol = 10)
  EdgeSet_heterodimer = Heterodimer_connection(peak_group, NodeSet, ppm_tol = 10, inten_threshold = 1e6)
  EdgeSet_experiment_MS2_fragment = Experiment_MS2_fragment_connection(peak_group, NodeSet, ppm_tol = 10, inten_threshold = 1e5)
  # FormulaSet needs to exist before running this line
  EdgeSet_library_MS2_fragment_df = Library_MS2_fragment_connection(peak_group, FormulaSet, MS2_library,
                                                                    inten_threshold = 1e5,
                                                                    ppm_tol = 10, abs_tol = 1e-4)
  EdgeSet_all_df = merge_edgeset(EdgeSet, 
                                 EdgeSet_ring_artifact, 
                                 EdgeSet_oligomer, 
                                 EdgeSet_heterodimer,
                                 EdgeSet_experiment_MS2_fragment,
                                 EdgeSet_library_MS2_fragment_df)
}

## Candidate formula pool
{
  
  FormulaSet = Propagate_formulaset(Mset, 
                                    NodeSet,
                                    FormulaSet,
                                    biotransform_step = 2,
                                    artifact_step = 3,
                                    propagation_ppm_threshold = 5e-6,
                                    propagation_abs_threshold = 2e-4,
                                    record_RT_tol = 0.15,
                                    record_ppm_tol = 5e-6)
  
  all_formulas = bind_rows(FormulaSet) %>%
    distinct(node_id, formula, category, .keep_all = T) 
  
}

print(Sys.time()-printtime)

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
                                            rule_score_biotransform = 0, rule_score_artifact = 0.5, 
                                            rule_score_oligomer = 0.5, rule_score_natural_abundance = 1,
                                            rule_score_fragment = 0.3, rule_score_ring_artifact = 1.5,
                                            rule_score_experiment_MS2_fragment = 1, rule_score_library_MS2_fragment = 0.3,
                                            rt_penalty_artifact_cutoff = 0.05, rt_penalty_artifact_ratio = 10, # This is a proxy for LC correlation
                                            inten_score_isotope = 1, 
                                            MS2_score_similarity = 1, MS2_similarity_cutoff = 0.3,
                                            MS2_score_experiment_fragment = 0.5)
  
  CplexSet[["heterodimer_ilp_edges"]] = score_heterodimer_ilp_edges(CplexSet, rule_score_heterodimer = 1,
                                                                    MS2_score_experiment_fragment = 0.5)
  
  CplexSet[["para"]] = Initiate_cplexset(CplexSet)
}
  
save.image(paste0(timestamp,".RData"))
## Run_cplex ####
{
  CplexSet[["init_solution"]] = list(Run_cplex(CplexSet, obj_cplex = CplexSet$para$obj, 
                                               relative_gap = 1e-4, total_run_time = 1500))
  # Test_para_CPLEX(CplexSet, obj_cplex = CplexSet$para$obj, 
  #                 para = c(0), para2 = NA, 
  #                 relative_gap = 1e-1, total_run_time = 3000)
}
print(Sys.time()-printtime)
## Read out CPLEX results ####
{
  CPLEX_all_x = Read_cplex_result(CplexSet$init_solution)
  # CPLEX_all_x = Read_CPLEX_result(CPLEXset$Pmt_solution)
  CPLEX_x = rowMeans(CPLEX_all_x,na.rm=T)
  
  ilp_nodes = CplexSet$ilp_nodes %>%
    mutate(ilp_result = CPLEX_x[1:nrow(.)]) %>%
    inner_join(Mset$Data %>% dplyr::select(id, Input_id), by = c("node_id"="id")) %>%
    filter(T)
  
  ilp_edges = CplexSet$ilp_edges %>%
    mutate(ilp_result = CPLEX_x[1:nrow(.) + nrow(ilp_nodes)]) %>%
    # filter(ilp_result == 1) %>%
    filter(T)
  
  heterodimer_ilp_edges = CplexSet$heterodimer_ilp_edges %>%
    mutate(ilp_result = CPLEX_x[1:nrow(.) + nrow(ilp_nodes) + nrow(ilp_edges)])
  
  ilp_nodes_ilp = ilp_nodes %>%
    filter(ilp_result > 0.01 | is.na(ilp_result)) %>%
    # filter(grepl("Thiamine", path))
    filter(T)
  
  ilp_edges_ilp = ilp_edges %>%
    filter(ilp_result > 0.01 | is.na(ilp_result)) %>%
    # filter(grepl("Thiamine", path))
    filter(T)
}

## View Results ####
{
  test = ilp_nodes %>% 
    filter(ilp_result > 0.01 | is.na(ilp_result))
  test2 = merge(WL, test, by.x = "Index", by.y = "Input_id", all.x = T) %>%
    dplyr::select(colnames(WL), path, formula, category, ppm_error, parent_formula, transform, node_id, parent_id, steps, temp_cat)
    
  # write_csv(test2, "WL_neg_NetID.csv", na="")
  
  test3 = test2 %>%
    # filter(is.na(Background)) %>%
    filter(is.na(Background)) %>%
    mutate(Formula = ifelse(is.na(Formula), "", Formula)) %>%
    mutate(Formula = check_chemform(isotopes, Formula)$new_formula)
  
  test3_filter = test3 %>%
    # filter(sig > 5) %>%
    # filter(category == "Natural_abundance", is.na(Isotope)) %>%
    filter(category != "Metabolite") %>%
    filter(Feature == "Metabolite") %>%
    # filter(!Feature %in% c("Buffer Sensitive", "Fragment_CID")) %>%
    filter(Feature...23 == "Metabolite") %>%
    filter(category != "Unknown") %>%
    mutate(N_count = sapply(formula, elem_num_query, "N"),
           C_count = sapply(formula, elem_num_query, "C")) %>%
    mutate(CN_match = (N_count == N) & (C_count == C)) %>%
    # filter(!is.na(C)) %>%
    # filter(CN_match) %>%
    filter(T)
  
  library(janitor)
  test4 = tabyl(test3, Feature)
  crosstable = tabyl(test3, Feature, category, Feature...23)
  crosstable = tabyl(test3, Feature...23, category)
  crosstable = tabyl(test3, Feature, temp_cat)
  crosstable = tabyl(test3 
                     %>% mutate(Feature2 = ifelse(is.na(Feature)|Feature == "Metabolite", Feature, "Artifact")), 
                     Feature2, temp_cat)
  # write_csv(crosstable, "crosstable_nonbio2_1e-4.csv")
  
}

# Path annotation
{
  # core_met annotation 
  # Does optimization result impact ranking? 
  # No, it does not. Optimization selects different formulas, the ranking selects structure.
  
  core_annotation = core_annotate(ilp_nodes, FormulaSet_df, LibrarySet)
  core_annotation_unique = core_annotate_unique(core_annotation)
  g_met = initiate_g_met(ilp_nodes, ilp_edges)
  
  met_dist_mat = initiate_met_dist_mat(g_met, ilp_nodes, core_annotation_unique)
  
  g_nonmet = initiate_g_nonmet(ilp_nodes, ilp_edges, heterodimer_ilp_edges)
  
  nonmet_dist_mat = initiate_nonmet_dist_mat(g_nonmet, ilp_nodes, core_annotation_unique)
  
  path_annotation = rep("", nrow(ilp_nodes))
  
  
  
  core_annotation_nonmet = as_data_frame(g_nonmet, "vertices") %>%
    filter(category != "Unknown") %>%
    filter(steps %% 1 == 0) %>%
    mutate(core_annotate = paste("Peak", node_id, formula)) %>%
    dplyr::rename(ilp_node_id = name) %>%
    mutate(ilp_node_id = as.integer(ilp_node_id))
  
  ilp_edges_annotate_met = as_data_frame(g_met)
  ilp_edges_annotate_nonmet = as_data_frame(g_nonmet)
  # options(warn = 2)
  tictoc::tic()
  for(i in 1:nrow(ilp_nodes)){
    print(i)
    # if(i == 1000)break
    if(ilp_nodes$temp_cat[i] == "Artifact"){
      
      
      path_annotation[i] = track_annotation_nonmet(i, ilp_edges_annotate = ilp_edges_annotate_nonmet,
                                                   g_nonmet, graph_path_mode = "out",
                                            nonmet_dist_mat, core_annotation_nonmet)
    } else if(ilp_nodes$temp_cat[i] == "Metabolite"){
      
      path_annotation[i] = track_annotation_met(i, ilp_edges_annotate = ilp_edges_annotate_met,
                                                g_met, graph_path_mode = "all",
                                                met_dist_mat, core_annotation_unique)
    } else {
      path_annotation[i] = "Unknown"
    }
  }
  tictoc::toc()
  
  i=2
  test = ilp_nodes %>%
    mutate(path = path_annotation) %>%
    # filter(path == "No edge connections.") %>%
    filter(ilp_result > 0.5) %>%
    # filter(category == "Metabolite") %>%
    # filter(steps != 0) %>%
    filter(T)
  test2 = FormulaSet_df %>%
    filter(node_id == 135)
  
  test_edge = ilp_edges %>%
    filter(node1 == 14, node2 == 126)
}



print("total run time")
save.image()
print(Sys.time()-printtime)

## Query specific nodes or edges ####
{
  Input_id = 603
  Mset$Data$id[Mset$Data$Input_id == Input_id]
  
  node_id_selected = 996
  ilp_nodes_selected = ilp_nodes %>%
    filter(node_id %in% node_id_selected)
  ilp_edges_node_related = ilp_edges %>%
    filter(node1 == node_id_selected | node2 == node_id_selected)%>%
    # filter(ilp_nodes1 %in% c(11486, 11488) | ilp_nodes2 %in% c(11486, 11488)) %>%
    arrange(-cplex_score, category) %>%
    filter(T)
  heterodimer_ilp_edges_node_related = heterodimer_ilp_edges %>%
    filter(node1 == node_id_selected | node2 == node_id_selected | linktype == node_id_selected)%>%
    arrange(-cplex_score, category) %>%
    filter(T)
  
}


# Evaluate parameter
{
  # How RT differ between connected peaks
  
  node_rt = sapply(NodeSet, "[[", "RT")
  test_rt_correct = test3 %>%
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
  test_mz_correct = test3 %>%
    filter(category == "Adduct") %>%
    filter(Feature != "Adduct") %>%
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
  fitdist_mz_correct = fitdistrplus::fitdist(test_mz_correct$ppm_mz_dif, "norm")
  plot(fitdist_mz_correct)
  print(fitdist_mz_correct)
  
}
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

# Deprecated ####
## Show path ####
# {
#   
#   ilp_nodes_merge = ilp_nodes %>%
#     dplyr::select(node_id, formula, temp_cat, ilp_node_id, ilp_result, score_prior_propagation, cplex_score)
#   
#   result = FormulaSet_df %>%
#     mutate(temp_cat = case_when(
#       category == "Metabolite" ~ "Metabolite", 
#       category == "Unknown" ~ "Unknown",
#       T ~ "Artifact")) %>%
#     merge(ilp_nodes_merge, all.x = T) %>%
#     filter(node_id == 5506)
#   # split(.$node_id)
#   
#   test = result[[5506]]
#   
#   result_ls = vector("list", length(NodeSet))
#   for(i in 1:length(NodeSet)){
#     if(!is.null(result[[as.character(i)]])){
#       result_ls[[i]] = result[[as.character(i)]]
#     } 
#   }
#   
#   node_path = sapply(unique(ilp_nodes$node_id), query_path, result_ls, LibrarySet)
#   path_annotation = data.frame(node_id = unique(ilp_nodes$node_id),
#                                path = node_path)
#   
#   ilp_nodes = ilp_nodes %>%
#     merge(path_annotation) %>%
#     dplyr::select(path, everything()) 
#   
#   # save(ilp_nodes, ilp_edges, heterodimer_ilp_edges, result_ls, file="network.RData")
# }