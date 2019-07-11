# Main Codes ####
## Read files ####
{
  Mset = list()
  Mset[["Library"]] = read_library("hmdb_CHNOPS.csv")
  Mset[["Biotransform"]]=Read_rule_table(rule_table_file = "biotransform.csv")
  Mset[["Artifacts"]]=Read_rule_table(rule_table_file = "artifacts.csv")
  
  setwd("./10tissues_pos_mode")
  filename = c("raw_data.csv")
  Mset[["Raw_data"]] <- read_csv(filename)
}

## Initialise ####
{
  Mset[["Global_parameter"]]=  list(mode = 1,
                                    normalized_to_col_median = F)
  Mset[["Cohort"]]=Cohort_Info(Mset)
  print(Mset$Cohort)
  
  #Clean-up duplicate peaks 
  Mset[["Data"]] = Peak_cleanup(Mset,
                                ms_dif_ppm=1/10^6, 
                                rt_dif_min=0.01,
                                detection_limit=500)
  
  Mset[["ID"]] = Mset$Data$ID
}

# Network ####
{
  read_from_csv = F
  EdgeSet = list()
  
  Mset[["NodeSet"]]=Form_node_list(Mset)
  
  
  EdgeSet[["Biotransform"]] = Edge_biotransform(Mset, 
                                                mass_abs = 0.001, 
                                                mass_ppm = 5/10^6, 
                                                read_from_csv = read_from_csv)
  
  adjust = Check_sys_measure_error(EdgeSet$Biotransform, inten_threshold=1e5)
  if(abs(adjust[1])>0.5 | abs(adjust[2])> 0.0001){
    Mset$Data$medMz = Mset$Data$medMz*(1+adjust[1]/10^6)+adjust[2]
    Mset[["NodeSet"]]=Form_node_list(Mset)
    EdgeSet[["Biotransform"]] = Edge_biotransform(Mset,
                                                  mass_abs = 0.001,
                                                  mass_ppm = 5/10^6,
                                                  read_from_csv = read_from_csv)
  }
  EdgeSet[["Biotransform"]] = Edge_score(EdgeSet$Biotransform)
  
  EdgeSet[["Peak_inten_correlation"]] = Peak_variance(Mset,
                                                      time_cutoff=0.1,
                                                      TIC_cutoff = 10000,
                                                      correlation_cutoff = -1)
  
  EdgeSet[["Artifacts"]] = Artifact_prediction(Mset, 
                                               EdgeSet$Peak_inten_correlation, 
                                               search_ms_cutoff=0.002,
                                               search_ppm_cutoff=10,
                                               read_from_csv = read_from_csv)
  EdgeSet[["Artifacts"]] = Edge_score(EdgeSet$Artifacts)
  
  #heterodimer  
  EdgeSet[["Heterodimer"]] = Hetero_dimer(EdgeSet$Peak_inten_correlation)
  
  
  EdgeSet[["Merge"]] = Merge_edgeset(EdgeSet)
  
  
  
  Mset[["NodeSet_network"]] = Network_prediction(Mset, 
                                                 edge_biotransform = EdgeSet$Biotransform, 
                                                 edge_artifact = EdgeSet$Artifacts,
                                                 biotransform_step = 4,
                                                 artifact_step = 5,
                                                 propagation_score_threshold = 0.5,
                                                 top_n = 50,
                                                 read_from_csv = read_from_csv)
  
  CPLEXset = Prepare_CPLEX(Mset, EdgeSet, read_from_csv = read_from_csv)
}

# Run CPLEX ####
{
  save.image("temp.RData")
  CPLEXset$data$unknown_formula = Score_formula(CPLEXset,
                                                rdbe=T, step_score=T, iso_penalty_score=F)
  edge_info_sum = Score_edge_cplex(CPLEXset, edge_bonus = 0.1)
  obj_cplex = c(CPLEXset$data$unknown_formula$cplex_score, edge_info_sum$edge_score)
}


## Feature generation ####
{
  #Identify peaks with high blanks
  Mset[["High_blanks"]]=High_blank(Mset, fold_cutoff = 2)
  
  #library_match
  Mset[["library_match"]] = library_match(Mset, ppm=5/10^6)
  
  #Metaboanalyst_Statistic
  if(length(unique(Mset$Cohort$sample_cohort))>1 & 
     min(table(Mset$Cohort$sample_cohort))>2){
    Mset[["Metaboanalyst_Statistic"]]=Metaboanalyst_Statistic(Mset)
  }
  
  # output assigned formula
  Mset[["Summary"]] = Summary_Mset(Mset)
}
save.image("final.RData")