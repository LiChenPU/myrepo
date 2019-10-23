
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))


NetID_files = list.files(pattern = ".RData")
for(repeating_n in 1:length(NetID_files)){
  NetID_files = list.files(pattern = ".RData")
  if(file.exists(sub(".RData", ".csv", NetID_files[repeating_n]))){
    next
  } else {
    print(NetID_files[repeating_n])
    load(NetID_files[repeating_n])
    # source("../NetID_CPLEX.R")
  }
}



for(repeating_n2 in 1:15){
  NetID_files = list.files(pattern = ".RData")
  if(file.exists(sub(".RData", " edge.csv", NetID_files[repeating_n2]))){
    next
  } else {
    print(NetID_files[repeating_n2])
    load(NetID_files[repeating_n2])
    edge_info_sum["ILP_result"] = CPLEX_x[(nrow(unknown_formula)+1):length(CPLEX_x)]
    write.csv(edge_info_sum, paste(timestamp, "edge.csv"), row.names = F)
  }
}

# 
# 
# library(dplyr)
# NetID_files = list.files(pattern = "[0-9].csv")
# formula_list_ls = list()
# for(i in 1:length(NetID_files)){
#   formula_list_ls[[length(formula_list_ls)+1]] = read.csv( NetID_files[i], stringsAsFactors = F) 
# }
# 
# NetID_files = list.files(pattern = "edge.csv")
# edge_info_sum_ls = list()
# for(i in 1:length(NetID_files)){
#   edge_info_sum_ls[[length(edge_info_sum_ls)+1]] = read.csv( NetID_files[i], stringsAsFactors = F)
# }

# 
# CPLEX_summary = formula_list_ls[[1]] %>%
#   filter(ILP_result!=0) %>%
#   dplyr::select(ID, formula)
# for(i in 2:length(formula_list_ls)){
#   CPLEX_summary = merge(CPLEX_summary,formula_list_ls[[i]]%>%
#                           filter(ILP_result!=0) %>%
#                           dplyr::select(ID, formula), all = T, by.x = "ID", by.y = "ID")
# }
# colnames(CPLEX_summary) = c("ID", paste("run",1:6,sep=""))
# 
# CPLEX_summary = CPLEX_summary %>%
#   mutate(all_same = run1==run2 & run1==run3 & run1==run4 & run1==run5 & run1==run6)
# 
# CPLEX_summary_inconsistent = CPLEX_summary %>%
#   filter(!all_same)
#   
#   
edge_info_sum = edge_info_sum_ls[[2]]
formula_list = formula_list_ls[[2]]

# edge_info_sum = g_edge
# formula_list = g_vertex

{
  selected_ID = 1743

  edge_selected = edge_info_sum %>% filter(node1==selected_ID|node2==selected_ID) #%>% filter(category != "biotransform")
  
  formula_selected = formula_list %>% filter(ID == selected_ID)
  formula_selected_in_edges = formula_list %>% filter(ID %in% c(edge_selected$node1, edge_selected$node2))
}


{
  test = edge_info_sum %>%
     filter(linktype == "H2O1") %>%
    distinct(edge_id, .keep_all = T) %>%
    filter(ILP_result ==1)
}

















