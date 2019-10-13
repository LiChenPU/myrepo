source("../NetID_function.R")
library(cplexAPI)

{
  CPLEXset[["Init_solution"]] = list(Run_CPLEX(CPLEXset, obj_cplex = CPLEXset$para$obj))
  # CPLEXset[["Screen_solution"]] = CPLEX_screen_edge(CPLEXset, edge_bonus_range = seq(-.6, -0.9, by=-0.1))
  # CPLEXset[["Pmt_solution"]] = CPLEX_permutation(CPLEXset, n_pmt = 20, sd_rel_max = 0.2)
  
  # Test_para_CPLEX(CPLEXset, obj_cplex = CPLEXset$para$obj, test_para =c(1,3))
}

# Read CPLEX result ####
{
  CPLEX_all_x = Read_CPLEX_result(CPLEXset$Init_solution)
  # CPLEX_all_x = Read_CPLEX_result(CPLEXset$Pmt_solution)
  
  CPLEX_x = rowMeans(CPLEX_all_x,na.rm=T)
  #CPLEX_x = result_solution$x
}

{
  unknown_nodes = CPLEXset$formula$unknown_nodes[,1:3]
  unknown_formula = CPLEXset$formula$unknown_formula %>% mutate(ILP_result = CPLEX_x[1:nrow(.)])
  unknown_formula_CPLEX = unknown_formula %>% filter(ILP_result!=0 )
  edge_info_sum = CPLEXset$edge
  edge_info_sum["ILP_result"] = 0
  edge_info_sum = edge_info_sum %>%
    filter(edge_score >=0) %>% 
    mutate(ILP_result = CPLEX_x[(nrow(unknown_formula)+1):length(CPLEX_x)]) %>%
    rbind(edge_info_sum) %>%
    distinct(edge_ilp_id, .keep_all = T) %>%
    arrange(edge_ilp_id)
  print(paste("pred formula num =", nrow(unknown_formula_CPLEX)))
  
  # unknown_node_CPLEX = merge(unknown_nodes,unknown_formula_CPLEX,by.x = "ID", by.y = "id",all=T)
}

{
  select_node_id = 909
  select_edge_info_sum = edge_info_sum %>%
    filter(node1 == select_node_id | node2 == select_node_id)
  
}

{
  select_ilp_id = 3555
  select_edge_info_sum = edge_info_sum %>%
    filter(ILP_id1 == select_ilp_id | ILP_id2 == select_ilp_id)
  
}

{
  pred_formula_ls[[1797]]
  pred_formula_ls[[5109]]
}


g_vertex_edge = determine_is_metabolite()
formula_list = g_vertex_edge$formula_list
relation_list = g_vertex_edge$relation_list

test = formula_list %>% 
  filter(ILP_result !=0) %>%
  filter(ILP_id <= nrow(unknown_formula))  %>%
  filter(is_metabolite == "NA")
table()

edge_info_sum["ILP_result"] = 0
edge_info_sum = edge_info_sum %>%
  filter(edge_score >=0) %>% 
  # filter(!category == "Heterodimer") %>%
  # filter(!category == "Ring_artifact") %>%
  mutate(ILP_result = CPLEX_x[(nrow(unknown_formula)+1):length(CPLEX_x)]) %>%
  rbind(edge_info_sum) %>%
  distinct(edge_ilp_id, .keep_all = T) %>%
  arrange(edge_ilp_id)
write.csv(edge_info_sum, paste(timestamp, "edge.csv"), row.names = F)
write.csv(formula_list2, paste0(timestamp, ".csv"), row.names = F)
save.image(paste0(timestamp,".RData"))






# # Graphic analysis ####
# {
#   g_vertex = formula_list2[!is.na(formula_list2$ILP_id),]
#   g_vertex = g_vertex[with(g_vertex, order(ID, -ILP_result)),]
#   g_edge = relation_list2
#   g <- graph_from_data_frame(d = g_edge, vertices = g_vertex, directed = T)
#   
#   write.csv(g_vertex, "g_vertex.txt", row.names = F)
#   write.csv(g_edge, "g_edge.txt", row.names = F)
# }
# 
# Mset[["Summary"]] = Summary_Mset(Mset)
# 
# mdata = Mset$Summary
# mdata = data.frame(ID=mdata[,1],formula=NA, ILP_result=NA,mdata[2:ncol(mdata)])
# if(!is.null(g_vertex)){
#   ilp_formula = g_vertex[g_vertex$ILP_result !=0,c("ID","formula","ILP_result", "is_metabolite")]
#   mdata = merge(ilp_formula, Mset$Summary, all.y = T)
# } 
# 
# write_csv(mdata, "mdata.csv")
# save.image("optimized.RData")
