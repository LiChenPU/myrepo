
source("../NetID_function.R")

{
  CPLEXset[["Init_solution"]] = list(Run_CPLEX(CPLEXset, obj_cplex))
  # CPLEXset[["Screen_solution"]] = CPLEX_screen_edge(CPLEXset, edge_bonus_range = seq(-.6, -0.9, by=-0.1))
  # CPLEXset[["Pmt_solution"]] = CPLEX_permutation(CPLEXset, n_pmt = 20, sd_rel_max = 0.2)
}

# Read CPLEX result ####
{
  CPLEX_all_x = Read_CPLEX_result(CPLEXset$Init_solution)
  # CPLEX_all_x = Read_CPLEX_result(CPLEXset$Pmt_solution)
  
  CPLEX_x = rowMeans(CPLEX_all_x,na.rm=T)
  #CPLEX_x = result_solution$x
}

{
  unknown_nodes = CPLEXset$data$unknown_nodes[,1:3]
  unknown_formula = CPLEXset$data$unknown_formula %>% mutate(ILP_result = CPLEX_x[1:nrow(.)])
  unknown_formula_CPLEX = unknown_formula %>% filter(ILP_result!=0 )

  print(paste("pred formula num =", nrow(unknown_formula_CPLEX)))
  
  # unknown_node_CPLEX = merge(unknown_nodes,unknown_formula_CPLEX,by.x = "ID", by.y = "id",all=T)
}

## determine_is_metabolite - A messy function so far, probably need to clean up
determine_is_metabolite = function(){
  data_peak_num = nrow(Mset$Data)
  
  {
    formula_list = merge(Mset$NodeSet, unknown_formula,by.x = "ID", by.y = "id",all=T)
    formula_list$formula[formula_list$ID>nrow(Mset$Data)] = formula_list$MF[formula_list$ID>nrow(Mset$Data)]
    formula_list$rdbe.y[formula_list$ID>nrow(Mset$Data)] = formula_list$rdbe.x[formula_list$ID>nrow(Mset$Data)]
    formula_list=formula_list[,!(colnames(formula_list) %in% c("MF", "rdbe.x", "is_metabolite"))]
    
    formula_list["ILP_id"]=NA
    formula_list$ILP_id[!is.na(formula_list$formula)] = 1: sum(!is.na(formula_list$formula))
    colors <- c("grey", "white", "red", "yellow", "green")
    formula_list["color"] = colors[formula_list$category+2]
    
    edge_info_sum["ILP_result"] = CPLEX_x[(nrow(unknown_formula)+1):length(CPLEX_x)]
    
    ilp_id_non0 = formula_list %>%
      filter(ILP_result!=0) %>%
      pull(ILP_id)
    
    
    edge_ilp_id_non0 = edge_info_sum %>% 
      filter(ILP_result==0) %>%
      filter((node1 > nrow(Mset$Data) & node2 %in% ilp_id_non0) | 
               (node2 > nrow(Mset$Data) & node1 %in% ilp_id_non0)) %>%
      pull(edge_ilp_id)
    
    relation_list = edge_info_sum %>%
      filter(ILP_result!=0 | edge_ilp_id %in% edge_ilp_id_non0)
  }

  ## Define if formula comes from artifact or biotransform
  {
    formula_list["is_artifact"]=FALSE
    artifact_edgeset = relation_list[relation_list$category!=1,]
    artifact_formula = unique(c(artifact_edgeset$ILP_id1[artifact_edgeset$direction==-1], 
                                artifact_edgeset$ILP_id2[artifact_edgeset$direction!=-1]))
    formula_list[formula_list$ILP_id %in% artifact_formula,"is_artifact"]=TRUE
    
    
    
    formula_list2 = cbind(formula_list[,c("ILP_id")], formula_list[,-which(colnames(formula_list) =="ILP_id")])
    colnames(formula_list2)[1] = "ILP_id"
    relation_list2 = cbind(relation_list[,c("ILP_id1","ILP_id2")], relation_list[,-which(colnames(relation_list) %in% c("ILP_id1","ILP_id2"))] )
    
    formula_list2["is_biotransform"]=FALSE
    biotranform_edgeset = relation_list2[relation_list2$category==1,]
    
    g_bio = graph_from_data_frame(d = biotranform_edgeset, 
                                  vertices =  formula_list2[formula_list2$ILP_id %in% 
                                                              c(biotranform_edgeset$ILP_id1, 
                                                                biotranform_edgeset$ILP_id2),],
                                  directed = F)
    clu=components(g_bio)
    #subnetwork criteria 
    # g_bio_subnetwork = igraph::groups(clu)[table(clu$membership)<10000]
    g_bio_subnetwork = igraph::groups(clu)
    test2 = as.data.frame(clu$membership)
    test3 = test2[as.numeric(row.names(test2))>nrow(unknown_formula),]
    biotranform_nodes_id = c()
    for(i in unique(test3)){
      biotranform_nodes_id = c(biotranform_nodes_id, as.numeric(g_bio_subnetwork[[i]]))
    }
    
    #biotranform_nodes = unique(c(biotranform_edgeset$node1, biotranform_edgeset$node2))
    formula_list2[!is.na(formula_list2$formula)& formula_list2$ILP_id %in% biotranform_nodes_id, "is_biotransform"]=TRUE
    
    formula_list2["is_metabolite"]=NA
    formula_list2$is_metabolite[formula_list2$is_artifact & formula_list2$is_biotransform] = "Maybe"
    formula_list2$is_metabolite[formula_list2$is_artifact & !formula_list2$is_biotransform] = "No"
    formula_list2$is_metabolite[!formula_list2$is_artifact & formula_list2$is_biotransform] = "Yes"
    formula_list2$is_metabolite[!formula_list2$is_artifact & !formula_list2$is_biotransform] = "NA"
  }
  return(list(formula_list2 = formula_list2, relation_list2 = relation_list2))
}
g_vertex_edge = determine_is_metabolite()
formula_list2 = g_vertex_edge$formula_list2
relation_list2 = g_vertex_edge$relation_list2

printtime = Sys.time()
timestamp = paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
write.csv(formula_list2, paste(timestamp, ".csv"), row.names = F)
save.image(paste(timestamp,".RData"))



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
