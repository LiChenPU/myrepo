


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("mz_calculator_function.R")

# pseudo codes

# build edges
# generate artificial formula (no isotopes), score based on mass accuracy
# propagate formula 
# optimize via CPLEX
# return list of mz + formula


ion_mode = 1
{
  Mset = list()
  Mset[["Library"]] = read_library("./dependent/library.csv")
  Mset[["Connect_rules"]]=Read_rule_table(rule_table_file = "./dependent/connect_rules.csv")
  
  Mset[["Raw_data"]] <- read.csv("raw_data.csv")
  
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
                                                 biotransform_step = 10,
                                                 propagation_score_threshold = 0.25,
                                                 top_n = 50
                                                 )
  
  CPLEXset = Prepare_CPLEX(Mset, EdgeSet)
}

# Run CPLEX ####
{
  CPLEXset$data$unknown_formula = Score_formula(CPLEXset,
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



test = Mset$NodeSet_network
edgetest = EdgeSet$Connect_rules




