# Main Codes ####
## Read files ####

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("NetID_function.R")

work_dir = "Lin_Yeast_Neg"
ion_mode = -1
Empirical_rules_file = "./dependent/Empirical_rules.csv"

MassDistsigma = 0.5
print(work_dir)
print(ion_mode)

printtime = Sys.time()
timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
# sink(paste(timestamp,"log.txt"))
# sink()


{
  Mset = list()
  Mset[["Library"]] = read.csv("./dependent/HMDB_CHNOPS_clean.csv", stringsAsFactors = F)
  Mset[["Empirical_rules"]]=Read_rule_table(rule_table_file = Empirical_rules_file)
  
  setwd(work_dir)
  filename = "raw_data.csv"
  Mset[["Raw_data"]] <- read_raw_data(filename)
}

## Initialise ####
{
  Mset[["Global_parameter"]]=  list(mode = ion_mode,
                                    normalized_to_col_median = F)
  Mset[["Cohort"]]=Cohort_Info(Mset, first_sample_col_num = 15)
  print(Mset$Cohort)
  
  #Clean-up duplicate peaks 
  Mset[["Data"]] = Peak_cleanup(Mset,
                                ms_dif_ppm=0/10^6, 
                                rt_dif_min=0.1,
                                detection_limit=0,
                                remove_high_blank_ratio = 0,
                                first_sample_col_num = 15)
  print(c(nrow(Mset$Raw_data), nrow(Mset$Data)))
  
  Mset[["ID"]] = Mset$Data$ID
}

## Initiate nodeset and edgeset
{
  NodeSet = Initiate_nodeset(Mset)
  
  EdgeSet = Initiate_edgeset(Mset, NodeSet, 
                             mass_abs = 0, mass_ppm = 10, 
                             nonbio_RT_tol = 0.1)
  
  LibrarySet = Initiate_libraryset(Mset)
  
  
}

## Candidate edge pool ####
{
  
  peak_group = Peak_grouping(NodeSet, RT_cutoff = 0.2, inten_cutoff = 1e4)
  Ring_artifact = Ring_artifact_connection(peak_group, 
                                           ppm_range_lb = 50, ppm_range_ub = 1000, ring_fold = 50, inten_threshold = 1e6)
  Oligomer = Oligomer_connection(peak_group, ppm_tol = 10)
  Heterodimer = Heterodimer_connection(peak_group, 
                                       ppm_tol = 10, inten_threshold = 1e5)
  
}

## Candidate formula pool
{
  FormulaSet = Initilize_formulaset(Mset, NodeSet, 
                                    LibrarySet,
                                    ppm_tol = 5e-6)
  
  FormulaSet = Propagate_formulaset(Mset, 
                                    NodeSet,
                                    FormulaSet,
                                    biotransform_step = 5,
                                    artifact_step = 5,
                                    propagation_ppm_threshold = 2,
                                    record_RT_tol = 0.1,
                                    record_ppm_tol = 5e-6,
                                    propagation_intensity_threshold = 2e4,
                                    max_formula_num = 1e6,
                                    top_n = 50)

  all_bind = bind_rows(FormulaSet) %>%
    # mutate(origin = 1) %>%
    # distinct(node_id,formula, .keep_all=T) %>%
    filter(steps != 0) %>%
    filter(T)
  
  # all_bind2 = bind_rows(sf) %>% 
  #   mutate(origin = 2) %>%
  #   bind_rows(all_bind) %>%
  #   distinct(node_id,formula, .keep_all=T) %>%
  #   filter(!node_id %in% all_bind$node_id)
}

Sys.time()-printtime



## CplexSet & Scoring ####
{
  profvis::profvis({
    
  
  CplexSet = list()
  CplexSet[["ilp_nodes"]] = initiate_ilp_nodes(FormulaSet) 
  CplexSet[["ilp_nodes"]] = score_ilp_nodes(CplexSet$ilp_nodes, MassDistsigma = MassDistsigma, 
                                            rdbe_score=F, step_score=F)
  
  CplexSet[["ilp_edges"]] = initiate_ilp_edges(EdgeSet, CplexSet$ilp_nodes)
  CplexSet[["ilp_edges"]] = score_ilp_edges(CplexSet$ilp_edges, NodeSet, MassDistsigma = MassDistsigma, 
                                            rule_score_biotransform = 0.5, rule_score_artifact = 1.5, inten_score_isotope = 1.5)
  CplexSet[["para"]] = Initiate_cplexset(CplexSet)
  })
  CplexSet[["init_solution"]] = list(Run_CPLEX(CplexSet, obj_cplex = CplexSet$para$obj))
}

## Results ####
{
  CPLEX_all_x = Read_cplex_result(CplexSet$init_solution)
  # CPLEX_all_x = Read_CPLEX_result(CPLEXset$Pmt_solution)
  CPLEX_x = rowMeans(CPLEX_all_x,na.rm=T)
  
  ilp_nodes = CplexSet$ilp_nodes %>%
    mutate(ilp_result = CPLEX_x[1:nrow(.)])
  
  ilp_edges = CplexSet$ilp_edges %>%
    mutate(ilp_result = CPLEX_x[(nrow(ilp_nodes)+1):length(CPLEX_x)])
}

## PAVE evaluation ####
{
  library(readxl)
  library(janitor)
  groundtruth = read_xlsx("./4-4-Table S10-Annotation of all peaks detected in S. cerevisiae and E. coli-xi_LC.xlsx",
                          sheet = "Yeast-neg-truth") %>%
    arrange(ID...1) %>%
    # dplyr::select(-c(36:67)) %>%
    dplyr::select(ID...1,
                  c(6:7,68:ncol(.)))
  
  ID_input_id_map = Mset$Data %>%
    dplyr::select(ID, Input_id)
  
  merge_result = ilp_nodes %>%
    merge(ID_input_id_map, by.x="node_id", by.y="ID") %>%
    merge(groundtruth, by.x="Input_id", by.y="ID...1")
  
  merge_result = merge_result %>%
    mutate(formula_match = formula == Formula...68)
  
  merge_result_filter = merge_result %>%
    filter(ilp_result != 0 | is.na(ilp_result)) %>%
    filter(!formula_match) %>%
    filter(!is.na(Formula_validated)) %>%
    filter(T)
  
  
  tabyl(merge_result, ilp_result, formula_match)
  
  test = ilp_edges %>%
    group_by(ilp_nodes1, ilp_nodes2) %>%
    filter(n()>1)
  
}


## Query specific nodes or edges
{
  node_id_selected = 951
  ilp_nodes_selected = ilp_nodes %>%
    filter(node_id == node_id_selected)
  ilp_edges_node_related = ilp_edges %>%
    filter(node1 == 871 | node2 == 871)
  
  edge_id_selected = 2191
  EdgeSet_selected = EdgeSet %>%
    filter(edge_id == edge_id_selected)
  ilp_nodes_edge_related = ilp_nodes %>%
    filter(node_id %in% c(EdgeSet_selected$node1, EdgeSet_selected$node2))
}

save.image(paste0(timestamp,".RData"))

print("total run time")
print(Sys.time()-printtime)



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