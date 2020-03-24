# Main Codes ####
## Read files ####

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("NetID_function.R")

work_dir = "Lin_Yeast_Neg"
ion_mode = -1
Empirical_rules_file = "./dependent/Empirical_rules.csv"

MassDistsigma = 0.5
LC_method = "Hilic_25min_QE"
print(work_dir)
print(ion_mode)

printtime = Sys.time()
timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
# sink(paste(timestamp,"log.txt"))
# sink()

{
  Mset = list()
  # Mset[["Library"]] = read.csv("./dependent/HMDB_CHNOPS_clean.csv", stringsAsFactors = F)
  Mset[["Library_HMDB"]] = read.csv("./dependent/hmdb_database.csv", stringsAsFactors = F)
  Mset[["Library_known"]] = read.csv("./dependent/known_library.csv", stringsAsFactors = F) %>%
    filter(!is.na(.[,eval(LC_method)]))
  
  Mset[["Empirical_rules"]]=Read_rule_table(rule_table_file = Empirical_rules_file)
  
  setwd(work_dir)
  filename = "raw_data.csv"
  Mset[["Raw_data"]] <- read_raw_data(filename)
  manual_library_file = "manual_library.csv"
  
  Mset[["Manual_library"]] = read_manual_library(manual_library_file)
}

## Initialise ####
{
  Mset[["Global_parameter"]]=  list(mode = ion_mode,
                                    LC_method = LC_method)
  Mset[["Cohort"]]=Cohort_Info(Mset, first_sample_col_num = 15)
  print(Mset$Cohort)
  
  #Clean-up duplicate peaks 
  Mset[["Data"]] = Peak_cleanup(Mset,
                                mz_tol=1/10^6, 
                                rt_tol=0.1,
                                inten_cutoff=0,
                                high_blank_cutoff = 0,
                                first_sample_col_num = 15)
  print(c(nrow(Mset$Raw_data), nrow(Mset$Data)))
}

## Initiate nodeset and edgeset ####
{
  NodeSet = Initiate_nodeset(Mset)
  
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
  EdgeSet_all_df = merge_edgeset(EdgeSet, EdgeSet_ring_artifact, EdgeSet_oligomer, EdgeSet_heterodimer)
}

## Candidate formula pool
{
  FormulaSet = Initilize_empty_formulaset(NodeSet)
  
  FormulaSet = Match_library_formulaset(FormulaSet, 
                                        Mset, NodeSet, 
                                        LibrarySet,
                                        expand = F,
                                        ppm_tol = 5e-6)
  
  FormulaSet = Propagate_formulaset(Mset, 
                                    NodeSet,
                                    FormulaSet,
                                    biotransform_step = 3,
                                    artifact_step = 4,
                                    propagation_ppm_threshold = 1e-6,
                                    record_RT_tol = 0.1,
                                    record_ppm_tol = 5e-6)

  all_bind = bind_rows(FormulaSet) %>%
    # filter(steps != 0) %>%
    distinct(formula, node_id, .keep_all = T) 
  
}
print(Sys.time()-printtime)
save.image()

## CplexSet & Scoring ####
{
  CplexSet = list()
  FormulaSet_df = Score_formulaset(FormulaSet,
                                   database_match = 0.6, 
                                   manual_match = 0.5,
                                   rt_match = 1, 
                                   known_rt_tol = 0.5,
                                   bio_decay = -1,
                                   artifact_decay = -0.5)
  
  CplexSet[["ilp_nodes"]] = initiate_ilp_nodes(FormulaSet_df) 
  CplexSet[["ilp_nodes"]] = score_ilp_nodes(CplexSet$ilp_nodes, MassDistsigma = MassDistsigma, 
                                            formula_score = 0)
  
  CplexSet[["ilp_edges"]] = initiate_ilp_edges(EdgeSet_all_df, CplexSet)
  
  CplexSet[["heterodimer_ilp_edges"]] = initiate_heterodimer_ilp_edges(EdgeSet_all_df, CplexSet, NodeSet)
  
  CplexSet[["ilp_edges"]] = score_ilp_edges(CplexSet, NodeSet,
                                            rule_score_biotransform = 0.05, rule_score_artifact = 0.5, 
                                            rule_score_oligomer = 0.5, rule_score_ring_artifact = 2,
                                            inten_score_isotope = 1)
  
  CplexSet[["heterodimer_ilp_edges"]] = score_heterodimer_ilp_edges(CplexSet, rule_score_heterodimer = 2)
  
  CplexSet[["para"]] = Initiate_cplexset(CplexSet)
}
  
save.image(paste0(timestamp,".RData"))

print("total run time")
print(Sys.time()-printtime)

## Run_cplex ####
{
  CplexSet[["init_solution"]] = list(Run_cplex(CplexSet, obj_cplex = CplexSet$para$obj))
}

## Results ####
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
  
  test = ilp_nodes %>%
    filter(ilp_result > 0.01 | is.na(ilp_result))
}

## View Results ####
{
  result = FormulaSet_df %>%
    merge(ilp_nodes %>% 
            dplyr::select(node_id, formula, category, ilp_node_id, ilp_result, cplex_score), 
          all.x = T) %>%
    # mutate(node_id = factor(node_id, levels = 1:length(NodeSet))) %>%
    split(.$node_id)
  
  result_ls = vector("list", length(NodeSet))
  for(i in 1:length(NodeSet)){
    if(!is.null(result[[as.character(i)]])){
      result_ls[[i]] = result[[as.character(i)]]
    } 
  }
  

  tictoc::tic()
  for(i in unique(ilp_nodes$node_id)){
    # test = sapply(, query_path, result_ls, LibrarySet)
    test = query_path(i, result_ls, LibrarySet)
    print(i)
  }

  tictoc::toc()
  
  
}


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

  id_Input_id_map = Mset$Data %>%
    dplyr::select(id, Input_id)

  merge_result = ilp_nodes %>%
    merge(groundtruth, by.x="Input_id", by.y="ID...1")

  merge_result = merge_result %>%
    mutate(formula_match = formula == Formula...68)

  merge_result_filter = merge_result %>%
    filter(ilp_result > 0.01 | is.na(ilp_result)) %>%
    filter(formula_match) %>%
    filter(!is.na(Formula_validated) | Formula_validated == "Y", Formula_validated != "?") %>%
    # arrange(-ilp_result) %>%
    # distinct(node_id, .keep_all=T) %>% 
    # filter(ilp_result < 0.01) %>% 
    # filter(!formula_match) %>%
    filter(T)
  
  test = tabyl(merge_result_filter, category, Status_validated)
  
  merge_result_filter2 = merge_result_filter %>%
    filter(category == "Metabolite" | Status_validated == "Metabolite", 
           category != Status_validated) %>%
    filter(T)
}

## Query specific nodes or edges
{
  node_id_selected = 6
  ilp_nodes_selected = ilp_nodes %>%
    filter(node_id == node_id_selected)
  ilp_edges_node_related = ilp_edges %>%
    filter(node1 == node_id_selected | node2 == node_id_selected)
  
  test = ilp_edges %>%
    group_by(node1, node2, formula1, formula2) %>%
    filter(n()>1)

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
