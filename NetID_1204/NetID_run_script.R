# Main Codes ####
## Read files ####

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("NetID_function.R")

work_dir = "Lin_Yeast_Neg"
ion_mode = -1
Empirical_rules_file = "./dependent/Empirical_rules.csv"

sigma = 0.5
print(work_dir)
print(ion_mode)

# printtime = Sys.time()
# timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
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
  
  test = EdgeSet %>% filter(linktype == "H-2K2")
  test = sapply(NodeSet, "[[", "mz")
  
}


## Candidate formula pool
{
  FormulaSet = Initilize_formulaset(Mset, NodeSet, 
                                    LibrarySet,
                                    ppm_tol = 5e-6)
  
  FormulaSet = Propagate_formulaset(Mset, 
                                    FormulaSet,
                                    ppm_tol = 5e-6)
  
  all_bind = bind_rows(FormulaSet)
  
}



# Network ####
{
  EdgeSet = list()
  
  Mset[["NodeSet"]]=Form_node_list(Mset)
  
  EdgeSet[["Biotransform"]] = Edge_biotransform(Mset, 
                                                mass_abs = 0, 
                                                mass_ppm = 5)
  
  adjust = Check_sys_measure_error(EdgeSet$Biotransform, inten_threshold=1e5)
  if(abs(adjust[1])>0.5 | abs(adjust[2])> 0.0001){
    Mset$Data$medMz = Mset$Data$medMz*(1+adjust[1]/10^6)+adjust[2]
    Mset[["NodeSet"]]=Form_node_list(Mset)
    EdgeSet[["Biotransform"]] = Edge_biotransform(Mset,
                                                  mass_abs = 0,
                                                  mass_ppm = 5)
  }
  
  mass_dist_sigma = sigma
  EdgeSet[["Biotransform"]] = Edge_score(EdgeSet$Biotransform, mass_dist_sigma = mass_dist_sigma)
  
  EdgeSet[["Peak_inten_correlation"]] = Peak_variance(Mset,
                                                      time_cutoff=0.1,
                                                      TIC_cutoff = 10000,
                                                      correlation_cutoff = -1)
  
  EdgeSet[["Artifacts"]] = Artifact_prediction(Mset, 
                                               EdgeSet$Peak_inten_correlation, 
                                               search_ms_cutoff=0,
                                               search_ppm_cutoff=10)
  EdgeSet[["Artifacts"]] = Edge_score(EdgeSet$Artifacts, mass_dist_sigma = mass_dist_sigma)
  
  #heterodimer  
  EdgeSet[["Heterodimer"]] = Hetero_dimer(EdgeSet$Peak_inten_correlation, ppm_tolerance = 5, inten_threshold = 1e6)
  EdgeSet[["Heterodimer"]] = Edge_score(EdgeSet$Heterodimer, mass_dist_sigma = mass_dist_sigma)
  
  # Mass_ring_artifact
  EdgeSet[["Ring_artifact"]] = Ring_artifact(Peak_inten_correlation = EdgeSet$Peak_inten_correlation, 
                                             ppm_range_lb = 50, 
                                             ppm_range_ub = 1000, 
                                             ring_fold = 50, 
                                             inten_threshold = 1e6)
  
  EdgeSet[["Merge"]] = Merge_edgeset(EdgeSet) 
  
  Mset[["NodeSet_network"]] = Formula_propagate(Mset, 
                                                 EdgeSet,
                                                 biotransform_step = 6,
                                                 artifact_step = 6,
                                                 propagation_score_threshold = 0.2,
                                                 propagation_intensity_threshold = 2e4,
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




# ## Feature generation ####
# {
#   #Identify peaks with high blanks
#   Mset[["High_blanks"]]=High_blank(Mset, fold_cutoff = 2)
#   
#   #library_match
#   Mset[["library_match"]] = library_match(Mset, ppm=5/10^6)
#   
#   #Metaboanalyst_Statistic
#   if(length(unique(Mset$Cohort$sample_cohort))>1 & 
#      min(table(Mset$Cohort$sample_cohort))>2){
#     Mset[["Metaboanalyst_Statistic"]]=Metaboanalyst_Statistic(Mset)
#   }
#   
#   # output assigned formula
#   Mset[["Summary"]] = Summary_Mset(Mset)
# }


save.image(paste0(timestamp,".RData"))

print("total run time")
print(Sys.time()-printtime)



