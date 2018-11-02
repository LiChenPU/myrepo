#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
#import library
{
  library(readr)
  library(rcdk)
  library(Matrix)
  library(dplyr)
  #install.packages("slam")
  library(slam)
  library(Rglpk)
  library(cplexAPI)
}

#Read data
{
  #setwd("C:/Users/lc8/Dropbox/HRMass_ID/")
  setwd("C:/Users/Li Chen/Dropbox/HRMass_ID/180629 Triwizard Tournament/integer linear programing/ilp_known_metabolite")
  raw_node_list=read_csv("merge_node_list.csv")
  raw_edge_list=read_csv("formula_network_edge_list.csv")
  raw_pred_formula=read_csv("All_formula_predict.csv")
  
}

#Clean up
{
  pred_formula = raw_pred_formula[!grepl("Illegal",raw_pred_formula$formula),]
  pred_formula = pred_formula[!duplicated(pred_formula[,c("id","formula")]),]

  #Select node that has only 1 formula as library nodes, 
  #For future applicatoin, should select nodes where degree==0
  lib_nodes = raw_node_list[raw_node_list$category==0, ]
  unknown_nodes = raw_node_list[raw_node_list$category!=0,]
  unknown_nodes = unknown_nodes[!is.na(unknown_nodes$Predict_formula),]
  
  lib_formula = pred_formula[pred_formula$id %in% lib_nodes$ID,]
  lib_formula = lib_formula[lib_formula$degree==0,]
  unknown_formula = pred_formula[pred_formula$id %in% unknown_nodes$ID,]
  
  merge_formula = rbind(lib_formula,unknown_formula)
  
  #Select edge list where it relates only to unknown nodes that lie in the main network
  edge_list = raw_edge_list[raw_edge_list$node1 %in% unknown_nodes$ID |
                              raw_edge_list$node2 %in% unknown_nodes$ID,]
  
  pred_formula_ls = list()
  for(n in 1: max(merge_formula$id)){
    pred_formula_ls[[n]]=merge_formula[merge_formula$id==n,]
  }
}



  

##Core codes

#Construct constraint matrix 
read_from_csv = 0
if(!read_from_csv)
{
  #Library nodes 
  
  # num_lib_nodes = nrow(lib_nodes)
  # triplet_library_node = data.frame(i=1:num_lib_nodes, j=1:num_lib_nodes, v=1)
  
  #Unknown nodes
  num_unknown_nodes = nrow(unknown_nodes)
  triplet_unknown_nodes_ls = list()
  temp_j=1
  
  for(n in 1:nrow(unknown_nodes)){
    temp_i = n
    temp_unknown_formula = unknown_formula[unknown_formula$id==unknown_nodes$ID[n],]
    if(nrow(temp_unknown_formula)==0){next()}
    triplet_unknown_nodes_ls[[n]]=data.frame(i=rep(temp_i, nrow(temp_unknown_formula)), j=temp_j:(temp_j+nrow(temp_unknown_formula)-1), v=1)
    temp_j = temp_j+nrow(temp_unknown_formula)
  }
  triplet_unknown_node = bind_rows(triplet_unknown_nodes_ls)
  
  
  #Edge list

  #function group information
  {
    fun_group_1 = data.frame(c("H2",
                               "C2H2O", #ACETYL
                               "C10H12N5O6P", #AMP
                               "C16H30O", #Palmityl
                               "C6H10O5", #Glycosyl
                               "C11H17NO9", #Sialic acid
                               "HPO3",
                               # "CO", #FORMYL
                               # "NH3",
                               # "NH3-O" ,#transaminase
                               # "OH-NH2", #amino group
                               # "C",
                               # "CH2O",
                               "C9H13N3O10P2", #CDP
                               "S",
                               "SO3",
                               "C2H4",
                               "O",
                               "NH",
                               "CH2",
                               "CO2",
                               "H2O"
    ),
    stringsAsFactors=F)
    for (i in 1:nrow(fun_group_1)){
      fun_group_1[i,2]=get.formula(fun_group_1[i,1])@mass
    }
    colnames(fun_group_1) = c("fun_group", "mass")
    fun_group_1[fun_group_1$fun_group=="NH3-O",2]=get.formula("NH3")@mass-get.formula("O")@mass
    fun_group_1[fun_group_1$fun_group=="OH-NH2",2]=get.formula("OH")@mass-get.formula("NH2")@mass
    fun_group_1=rbind(c("Same", 0),fun_group_1)
    fun_group_1$mass=as.numeric(fun_group_1$mass)
  }
  
  
  time = Sys.time()
  temp_i=temp_j=1
  triplet_edge_ls_edge=triplet_edge_ls_node=list()
  edge_info = list()
  n=1
  for(n in 1:nrow(edge_list)){
  #for(n in 1:2000){
    temp_edge = edge_list[n,]
    node_1 = temp_edge$node1
    node_2 = temp_edge$node2
    formula_1 = pred_formula_ls[[node_1]]
    formula_2 = pred_formula_ls[[node_2]]
    temp_fg = fun_group_1$fun_group[temp_edge$linktype]
    temp_formula = unique(formula_1$formula[1])
    for(temp_formula in unique(formula_1$formula)){
      #Assuming formula in node_1 is always smaller than node_2
      if(temp_fg=="Same"){
        temp_formula_2=temp_formula
      }
      else{
        #temp_formula_2 = formula_manipulate(temp_formula, temp_fg, +1)
        temp_formula_2 = get.formula(paste(temp_formula,"+", temp_fg))@string
      }
      
      #Write triplet for edge and corresponding 2 nodes
      if(temp_formula_2 %in% formula_2$formula){
        temp_j1 = which(merge_formula$formula==temp_formula&merge_formula$id==node_1)
        temp_j2 = which(merge_formula$formula==temp_formula_2&merge_formula$id==node_2)

        #if one node is library node,
        if(temp_j1<=nrow(lib_nodes)){
          triplet_edge_ls_edge[[temp_i]] = data.frame(i=temp_i,
                                                      j=temp_j,
                                                      v=1)
          triplet_edge_ls_node[[temp_i]] = data.frame(i=temp_i,
                                                      j=temp_j2-nrow(lib_nodes),
                                                      v=-1)
          edge_info[[temp_i]] = data.frame(edge_id=n,
                                           edge_score=temp_edge$edge_massdif_score,
                                           temp_formula,
                                           temp_formula_2,
                                           temp_j1,
                                           temp_j2,
                                           stringsAsFactors=F
          )
        }
        if(temp_j2<=nrow(lib_nodes)){
          triplet_edge_ls_edge[[temp_i]] = data.frame(i=temp_i,
                                                      j=temp_j,
                                                      v=1)
          triplet_edge_ls_node[[temp_i]] = data.frame(i=temp_i,
                                                      j=temp_j1-nrow(lib_nodes),
                                                      v=-1)
          edge_info[[temp_i]] = data.frame(edge_id=n,
                                           edge_score=temp_edge$edge_massdif_score,
                                           temp_formula,
                                           temp_formula_2,
                                           temp_j1,
                                           temp_j2,
                                           stringsAsFactors=F
          )
        }
        
        #if both nodes are unknown nodes
        if(temp_j2>nrow(lib_nodes)&temp_j1>nrow(lib_nodes)){
          triplet_edge_ls_edge[[temp_i]] = data.frame(i=temp_i,
                                            j=temp_j,
                                            v=2)
          triplet_edge_ls_node[[temp_i]] = data.frame(i=rep(temp_i,2),
                                            j=c(temp_j1-nrow(lib_nodes),temp_j2-nrow(lib_nodes)),
                                            v=-1)
          edge_info[[temp_i]] = data.frame(edge_id=n,
                                           edge_score=temp_edge$edge_massdif_score,
                                           temp_formula,
                                           temp_formula_2,
                                           temp_j1,
                                           temp_j2,
                                           stringsAsFactors=F
                                           )
        }
        temp_i = temp_i+1
        temp_j = temp_j+1
      }
    } 
    #which(merge_formula$formula==temp_formula_2&merge_formula$id==node_2)
  }
  triplet_edge_ls_edge_sum = bind_rows(triplet_edge_ls_edge)
  triplet_edge_ls_edge_sum$i=triplet_edge_ls_edge_sum$i+num_unknown_nodes
  triplet_edge_ls_edge_sum$j=triplet_edge_ls_edge_sum$j+nrow(unknown_formula)
  triplet_edge_ls_node_sum = bind_rows(triplet_edge_ls_node)
  triplet_edge_ls_node_sum$i=triplet_edge_ls_node_sum$i+num_unknown_nodes
  edge_info_sum = bind_rows(edge_info)
  
  test=edge_list[-edge_info_sum$edge_id,]
  
  
  time = Sys.time()-time
  print(time)
  
  
  

  
  #Generate sparse matrix on left hand side
  triplet_df = rbind(
                     triplet_unknown_node,
                     triplet_edge_ls_edge_sum,
                     triplet_edge_ls_node_sum
                     )
  #write.csv(triplet_df,"triplet_df.csv",row.names = F)
  mat = simple_triplet_matrix(i=triplet_df$i,
                              j=triplet_df$j,
                              v=triplet_df$v)
}else
{
  triplet_df = read_csv("triplet_df.csv")
  mat = simple_triplet_matrix(i=triplet_df$i,
                              j=triplet_df$j,
                              v=triplet_df$v)
  
}

#Other parameters
{
  #should have number of columns, i.e. j
  obj <- c(rep(0, nrow(unknown_formula)), edge_info_sum$edge_score)
  types <- c(rep("B",nrow(unknown_formula)+nrow(edge_info_sum)))
  
  #should have number of rows, i.e. i
  dir <- c(rep("==",num_unknown_nodes), rep("<=", nrow(edge_info_sum)))
  rhs = c(rep(1,num_unknown_nodes),rep(0,nrow(edge_info_sum)))
  max <- TRUE
}  


# #R glpk Solver
# ILP_result = Rglpk_solve_LP(obj, mat, dir, rhs, types = types, max = max, 
#                             control = list(tm_limit=60*1000,canonicalize_status=F,verbose=T))
# para = data.frame(ILP_result$solution)
# unknown_formula["ILP_result"] = para$ILP_result.solution[1:nrow(unknown_formula)]
# edge_info_sum["ILP_result"] = para$ILP_result.solution[(nrow(unknown_formula)+1):nrow(para)]
# 
# # test=data.frame(ILP_result$solution)
# # test=cbind(test,ILP_result$solution)
# # test2=test[test$ILP_result.solution!=test$`ILP_result$solution`,]
# 
# time = Sys.time()
# print(time)



#CPLEX solver
{
  env <- openEnvCPLEX()
  prob <- initProbCPLEX(env)
  chgProbNameCPLEX(env, prob, "sample")
  nc <- max(mat$j)
  nr <- max(mat$i)
  obj <- c(rep(0, nrow(unknown_formula)), edge_info_sum$edge_score)
  rhs = c(rep(1,nrow(unknown_nodes)),rep(0,nrow(edge_info_sum)))
  sense <- c(rep("E",nrow(unknown_nodes)), rep("L", nrow(edge_info_sum)))
  lb <- rep(0, nc)
  ub <- rep(1, nc)
  # cn <- c("x1", "x2", "x3")
  # rn <- c("q1", "q2", "q3")
  
  
  triplet_df=triplet_df[with(triplet_df,order(j)),]
  cnt=as.vector(table(triplet_df$j))
  beg=vector()
  beg[1]=0
  for(i in 2:length(cnt)){beg[i]=beg[i-1]+cnt[i-1]}
  ind=triplet_df$i-1
  val = triplet_df$v
  
  
  #xctype=rep(CPX_BINARY,nc)
  #  checkCopyColTypeCPLEX(env, prob, xctype)
  
  
  copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj, rhs, sense,
                    beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
  writeProbCPLEX(env, prob, "prob.lp")
  #copyColTypeCPLEX(env, prob, xctype)  
  
}
{
  #primoptCPLEX(env, prob)
  lpoptCPLEX(env, prob)
  #dualoptCPLEX(env, prob)
  #mipoptCPLEX(env, prob)
  
  test=solutionCPLEX(env, prob)
  writeProbCPLEX(env, prob, "prob.lp")
  
  lp <- initProbCPLEX(env)
  readCopyProbCPLEX(env, lp, "prob.lp")
  
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
}



