# Main Codes ####
## Read files ####

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
printtime = Sys.time()
timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
source("NetID_function.R")

work_dir = "WL_liver_neg"
MS2_filepath = "./MS2"
MS2_library_file = "./dependent/HMDB_pred_MS2_neg3.rds"
MS2_library = readRDS(MS2_library_file)

# work_dir = "Lin_Yeast_Neg"
ion_mode = -1
mass_dist_gamma_rate = 1
LC_method = "Hilic_25min_QE"
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
  
  Mset[["Global_parameter"]]=  list(mode = ion_mode,
                                    LC_method = LC_method)
}

## WL data format ####
{
  setwd(work_dir)
  
  filename = "raw_data.csv"
  Mset[["Raw_data"]] <- read_raw_data(filename)
  
  filename_wl = "pks_liver_neg_buff_2020-02-21.xlsx"
  WL = readxl::read_xlsx(filename_wl) %>%
    dplyr::rename(medRt = rt,
                  medMz = mz)
  
  test = Mset$Raw_data %>%
    merge(WL, all = T) %>%
    mutate(id = 1:nrow(.),
           liver = 10^sig,
           liver2 = 10^sig) %>%
    filter(is.na(Background)) %>%
    dplyr::select(colnames(Mset$Raw_data))
  Mset[["Raw_data"]] = test
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
                                high_blank_cutoff = 2,
                                first_sample_col_num = 15)
  
  print(c(nrow(Mset$Raw_data), nrow(Mset$Data)))
  
  manual_library_file = "manual_library.csv"
  Mset[["Manual_library"]] = read_manual_library(manual_library_file)
}



## Initiate nodeset and edgeset ####
{
  NodeSet = Initiate_nodeset(Mset)
  NodeSet = Add_MS2_nodeset(MS2_filepath = MS2_filepath, 
                            NodeSet)
  # node_MS2_logic = sapply(NodeSet, function(x){!is.null(x$MS2)})
  # node_MS2 = (1:length(NodeSet))[node_MS2_logic]
  
  EdgeSet = Initiate_edgeset(Mset, NodeSet, 
                             mz_tol_abs = 0, mz_tol_ppm = 10, 
                             rt_tol_bio = Inf, rt_tol_nonbio = 0.2)
  
  LibrarySet = Initiate_libraryset(Mset)
}


## Extension of EdgeSet ####
{
  peak_group = Peak_grouping(NodeSet, RT_cutoff = 0.2, inten_cutoff = 1e4)
  EdgeSet_ring_artifact = Ring_artifact_connection(peak_group, ppm_range_lb = 50, ppm_range_ub = 1000, ring_fold = 50, inten_threshold = 1e6)
  EdgeSet_oligomer = Oligomer_connection(peak_group, ppm_tol = 10)
  EdgeSet_heterodimer = Heterodimer_connection(peak_group, NodeSet, ppm_tol = 10, inten_threshold = 1e6)
  EdgeSet_experiment_MS2_fragment = Experiment_MS2_fragment_connection(peak_group, NodeSet, ppm_tol = 10, inten_threshold = 1e5)
}

## Candidate formula pool
{
  FormulaSet = Initilize_empty_formulaset(NodeSet)
  
  FormulaSet = Match_library_formulaset(FormulaSet, 
                                        Mset, NodeSet, 
                                        LibrarySet,
                                        expand = T,
                                        ppm_tol = 5e-6)
  
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
  
  FormulaSet = Propagate_formulaset(Mset, 
                                    NodeSet,
                                    FormulaSet,
                                    biotransform_step = 2,
                                    artifact_step = 3,
                                    propagation_ppm_threshold = 2e-6,
                                    record_RT_tol = 0.1,
                                    record_ppm_tol = 5e-6)

  all_bind = bind_rows(FormulaSet) %>%
    # filter(steps != 0) %>%
    distinct(formula, node_id, .keep_all = T)
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
  
  CplexSet[["ilp_nodes"]] = initiate_ilp_nodes(FormulaSet_df) 
  CplexSet[["ilp_nodes"]] = score_ilp_nodes(CplexSet$ilp_nodes, mass_dist_gamma_rate = mass_dist_gamma_rate, 
                                            formula_score = 0, unassigned_penalty = -0.5)
  
  CplexSet[["ilp_edges"]] = initiate_ilp_edges(EdgeSet_all_df, CplexSet)
  
  CplexSet[["heterodimer_ilp_edges"]] = initiate_heterodimer_ilp_edges(EdgeSet_all_df, CplexSet, NodeSet)
  
  CplexSet[["ilp_edges"]] = score_ilp_edges(CplexSet, NodeSet,
                                            rule_score_biotransform = 0, rule_score_artifact = 0.5, 
                                            rule_score_oligomer = 0.5, rule_score_natural_abundance = 1,
                                            rule_score_ring_artifact = 2,
                                            rule_score_experiment_MS2_fragment = 1, rule_score_library_MS2_fragment = 0.3,
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
                                               relative_gap = 1e-4, total_run_time = 2000))
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

## Show path ####
{
  result = FormulaSet_df %>%
    mutate(temp_cat = ifelse(category == "Metabolite", "Metabolite", "Artifact")) %>%
    merge(ilp_nodes %>% 
            mutate(temp_cat = ifelse(category == "Metabolite", "Metabolite", "Artifact")) %>%
            dplyr::select(node_id, formula, temp_cat, ilp_node_id, ilp_result, cplex_score), 
          all.x = T) %>%
    # mutate(node_id = factor(node_id, levels = 1:length(NodeSet))) %>%
    dplyr::select(-temp_cat) %>%
    split(.$node_id)
  
  result_ls = vector("list", length(NodeSet))
  for(i in 1:length(NodeSet)){
    if(!is.null(result[[as.character(i)]])){
      result_ls[[i]] = result[[as.character(i)]]
    } 
  }
  
  node_path = sapply(unique(ilp_nodes$node_id), query_path, result_ls, LibrarySet)
  path_annotation = data.frame(node_id = unique(ilp_nodes$node_id),
                               path = node_path)

  ilp_nodes = ilp_nodes %>%
    merge(path_annotation) %>%
    dplyr::select(path, everything()) 
  
  ilp_nodes2 = ilp_nodes %>%
    filter(ilp_result > 0.01 | is.na(ilp_result)) %>%
    # filter(grepl("Thiamine", path))
    filter(T)
  
  # test = ilp_nodes %>%
  #   filter(node_id == 4292)
  
  # Show all duplicate edges are due to Library_MS2_fragment
  # test = ilp_edges %>%
  #   filter(ilp_result>0.1) %>%
  #   filter(category != "Library_MS2_fragment") %>%
  #   group_by(node1, node2) %>%
  #   filter(n()>1) %>%
  #   ungroup() %>%
  #   group_by(node1, node2, category) %>%
  #   filter(n()==1)
  # save(ilp_nodes, ilp_edges, heterodimer_ilp_edges, result_ls, file="network.RData")
}

## View Results ####
{
  test = ilp_nodes %>% 
    filter(ilp_result > 0.01 | is.na(ilp_result))
  test2 = merge(WL, test, by.x = "Index", by.y = "Input_id", all.x = T) %>%
    dplyr::select(colnames(WL), path, formula, category, ppm_error, parent_formula, node_id, parent_id, steps)
    # dplyr::select(colnames(WL), formula, category, parent_formula, node_id, parent_id, steps)
  
  # write_csv(test2, "WL_neg_NetID.csv", na="")
  
  test3 = test2 %>%
    filter(is.na(Background)) %>%
    # dplyr::select(colnames(WL)[colnames(WL)!="Index"], formula, everything()) %>%
    mutate(Formula = ifelse(is.na(Formula), "", Formula)) %>%
    mutate(Formula = check_chemform(isotopes, Formula)$new_formula)
  
  test3_filter = test3 %>%
    filter(sig > 4) %>%
    filter(category == "Natural_abundance", Feature != "Isotope") %>%
    # filter(steps >= 1) %>%
    # filter(Feature == "BuffSS") %>%
    # filter(Feature == "Metabolite" | is.na(Feature)) %>%
    # filter(is.na(Feature)) %>%
    # filter(Feature == "Metabolite") %>%
    # filter(category == "Library_MS2_fragment") %>%
    # filter(category == "Experirmental_MS2_fragment") %>%
    # filter(category == "Metabolite") %>%
    # filter(!(Feature == "Isotope" & category == "Artifact")) %>%
    # filter(Formula == formula & Formula != "") %>%
    # filter(!is.na(formula)) %>%
    # filter(Formula != formula) %>%
    # filter(parent_formula != formula) %>%
    filter(T)
  test3_filter2 = test3_filter %>%
    dplyr::select(Index, medMz, medRt, sig, Formula, Feature, 
                  path, formula, category, parent_formula, node_id, parent_id, steps) %>%
    filter(node_id %in% c(5727, 4818, 3654, 3304, 7149, 661, 4036, 6302, 8030))
  test3_filter3 = ilp_nodes_ilp %>%
    filter(node_id %in% c(5727, 4818, 3654, 3304, 7149, 661, 4036, 6302, 8030))
  
  test3_filter_fix = test3 %>%
    filter(sig > 5) %>%
    filter(is.na(formula)) %>%
    filter(Feature == "Isotope") %>%
    filter(T)
  
  library(janitor)
  test4 = tabyl(test3, category)
  test4 = tabyl(test3, Feature)
  crosstable = tabyl(test3, Feature, category)
  crosstable = tabyl(test3, category, Feature)
  crosstable_filter = tabyl(test3_filter, Feature, category)
  
  test4 = test3 %>%
    mutate(category2 = case_when(
      is.na(category) ~ "Unannotated",
      category != "Metabolite" ~ "Artifact",
      T ~ "Metabolite"
    )) %>%
    mutate(Feature2 = case_when(
      is.na(Feature) ~ "Unannotated",
      Feature != "Metabolite" ~ "Artifact",
      T ~ "Metabolite"
    ))
  
  tabyl(test4, Feature, Feature2)
  crosstable = tabyl(test4, category2, Feature2)
  
  test5 = test3 %>%
    filter(T)

}

print("total run time")
save.image()
print(Sys.time()-printtime)

# ## PAVE evaluation ####
{
  library(readxl)
  library(janitor)
  groundtruth = read_xlsx("./4-4-Table S10-Annotation of all peaks detected in S. cerevisiae and E. coli-xi_LC.xlsx",
                          sheet = "Yeast-neg-truth") %>%
    arrange(ID...1) %>%
    # dplyr::select(-c(36:67)) %>%
    dplyr::select(ID...1,
                  c(6:7,68:ncol(.)))

  merge_result = ilp_nodes %>%
    merge(groundtruth, by.x="Input_id", by.y="ID...1")
    

  merge_result = merge_result %>%
    mutate(formula_match = formula == Formula...68)

  merge_result_filter = merge_result %>%
    filter(ilp_result > 0.01 | is.na(ilp_result)) %>%
    filter(formula_match) %>%
    # filter(!is.na(Formula_validated) | Formula_validated == "Y") %>%
    filter(!is.na(Formula_validated) | Formula_validated == "Y", Formula_validated != "?") %>%
    # arrange(-ilp_result) %>%
    # distinct(node_id, .keep_all=T) %>% 
    # filter(ilp_result < 0.01) %>% 
    # filter(!formula_match) %>%
    dplyr::select(path, everything()) %>%
    filter(Status_validated == "Adduct", category == "Metabolite")
    filter(T)
  
  # test = tabyl(merge_result_filter, category, Status_validated)
  
  # merge_result_filter2 = merge_result_filter %>%
  #   filter(category == "Metabolite" | Status_validated == "Metabolite", 
  #          category != Status_validated) %>%
  #   filter(T)
}

## Query specific nodes or edges ####
{
  Input_id = 768
  Mset$Data$id[Mset$Data$Input_id == Input_id]
  
  node_id_selected = 688
  ilp_nodes_selected = ilp_nodes %>%
    filter(node_id == node_id_selected)
  ilp_edges_node_related = ilp_edges %>%
    filter(node1 == node_id_selected | node2 == node_id_selected)
  
  test = EdgeSet_all_df %>%
    filter(node2 == 1869)
  
  edge_id_selected = 2191
  EdgeSet_selected = EdgeSet %>%
    filter(edge_id == edge_id_selected)
  ilp_nodes_edge_related = ilp_nodes %>%
    filter(node_id %in% c(EdgeSet_selected$node1, EdgeSet_selected$node2))
}

## deprecated ####
# {
#   ConnectionSet = list()
#   
#   FormulaSet_df = bind_rows(FormulaSet)
#   library_edge = FormulaSet_df %>%
#     filter(RT == -1) %>% # parent RT = -1 means it is library
#     mutate(edge_ID = 1:nrow(.),
#            edge_category = "library_edge") %>%
#     dplyr::rename(node1 = parent_ID, 
#                   node2 = node_id,
#                   formula1 = parent,
#                   formula2 = formula) %>%
#     dplyr::select(edge_ID, node1, node2, formula1, formula2, transform, direction, edge_category) %>%
#     filter(T)
#   
#   data_edge = FormulaSet_df %>%
#     filter(RT != -1) %>% # parent RT != -1 means it is data edge
#     mutate(edge_ID = 1:nrow(.),
#            edge_category = "data_edge") %>%
#     dplyr::rename(node1 = parent_ID, 
#                   node2 = node_id,
#                   formula1 = parent,
#                   formula2 = formula) %>%
#     dplyr::select(edge_ID, node1, node2, formula1, formula2, transform, direction, edge_category) %>%
#     filter(T)
#   
#   ConnectionSet[["library_edge"]] = library_edge
#   ConnectionSet[["data_edge"]] = data_edge
#   
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

