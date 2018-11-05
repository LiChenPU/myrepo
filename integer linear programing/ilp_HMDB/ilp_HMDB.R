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
}

#Read data
{
  #setwd("C:/Users/lc8/Dropbox/HRMass_ID/")
  setwd("C:/Users/Li Chen/Dropbox/HRMass_ID/180629 Triwizard Tournament/integer linear programing")
  raw_node_list=read_csv("HMDB_node_list2.csv")
  raw_edge_list=read_csv("HMDB_edge_list.csv")
  raw_pred_formula=read_csv("HMDB_predicted_formula.csv")
  
}

#Clean up
{
  
  pred_formula = raw_pred_formula[!grepl("Illegal",raw_pred_formula$formula),]
  
  #Select the main network where formula is matched
  temp_pred_formula = pred_formula[!duplicated(pred_formula$id),c("id","formula")]
  node_list=merge(raw_node_list,temp_pred_formula,by.x="ID",by.y = "id",all=F)
  node_list_match = node_list[node_list$MF==node_list$formula,]
  #pred_dif = node_list[node_list$MF!=node_list$formula,]
  
  #Select formula that belong to main network nodes
  pred_formula = pred_formula[pred_formula$id %in% node_list_match$ID,]
  pred_formula = pred_formula[!duplicated(pred_formula[,c("id","formula")]),]
  
  
  #Select node that has only 1 formula as library nodes, 
  #For future applicatoin, should select nodes where degree==0
  lib_nodes = node_list_match[!(node_list_match$ID %in% pred_formula$id[duplicated(pred_formula$id)]), ]
  unknown_nodes = node_list_match[!(node_list_match$ID %in% lib_nodes$ID),]
  lib_formula = pred_formula[pred_formula$id %in% lib_nodes$ID,]
  unknown_formula = pred_formula[pred_formula$id %in% unknown_nodes$ID,]
  merge_nodes = rbind(lib_nodes,unknown_nodes)
  merge_formula = rbind(lib_formula,unknown_formula)
  
  #Select edge list where it relates only to unknown nodes that lie in the main network
  edge_list = raw_edge_list[raw_edge_list$node1 %in% unknown_nodes$ID |
                              raw_edge_list$node2 %in% unknown_nodes$ID,]
  
  pred_formula_ls = list()
  for(n in 1: max(pred_formula$id)){
    pred_formula_ls[[n]]=pred_formula[pred_formula$id==n,]
  }
}



  

##Core codes

#Construct constraint matrix 
read_from_csv = 0
if(!read_from_csv)
{
  #Library nodes 
  
  num_lib_nodes = nrow(lib_nodes)
  triplet_library_node = data.frame(i=1:num_lib_nodes, j=1:num_lib_nodes, v=1)
  
  #Unknown nodes
  num_unknown_nodes = nrow(unknown_nodes)
  triplet_unknown_nodes_ls = list()
  temp_j=num_lib_nodes+1
  for(n in 1:nrow(unknown_nodes)){
    temp_i = num_lib_nodes+n
    temp_unknown_formula = unknown_formula[unknown_formula$id==unknown_nodes$ID[n],]
    triplet_unknown_nodes_ls[[n]]=data.frame(i=rep(temp_i, nrow(temp_unknown_formula)), j=temp_j:(temp_j+nrow(temp_unknown_formula)-1), v=1)
    temp_j = temp_j+nrow(temp_unknown_formula)
  }
  triplet_unknown_node = bind_rows(triplet_unknown_nodes_ls)
  
  
  #Edge list
  
  #formula_manipulate, functions for process edge list
  {
    formula_table = function(X){
      Xformula = get.formula(X)
      temp=Xformula@isotopes
      df = data.frame(as.numeric(temp[,2]))
      rownames(df)=temp[,1]
      colnames(df)=Xformula@string
      return(df)
    }
    table_to_formula = function(result_formula){
      r = as.character()
      for(i in 1:nrow(result_formula)){
        if(result_formula$result[i]!=0){
          if(result_formula$result[i]==1){r=paste(r, result_formula$Row.names[i], sep="")}
          else{r=paste(r, result_formula$Row.names[i],  result_formula$result[i], sep="")}
        }
      }
      if(sum(result_formula$result<0)!=0){
        r=paste(r, "Illegal_formula")
      }
      return (r)
    }
    formula_manipulate = function(F1, F2, plusminus){
      temp_1 = formula_table(F1)
      temp_2 = formula_table(F2)
      result_formula = merge(temp_1,temp_2, all=T,  by="row.names")
      result_formula[is.na(result_formula)]=0
      result_formula["result"]=result_formula[,2]+result_formula[,3]*plusminus
      return (table_to_formula(result_formula))
    }
  }
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
  for(n in 1:nrow(edge_list)){
  #for(n in 1:2000){
    temp_edge = edge_list[n,]
    node_1 = temp_edge$node1
    node_2 = temp_edge$node2
    formula_1 = pred_formula_ls[[node_1]]
    formula_2 = pred_formula_ls[[node_2]]
    temp_fg = fun_group_1$fun_group[temp_edge$linktype]
    
    for(temp_formula in formula_1$formula){
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
        triplet_edge_ls_edge[[temp_i]] = data.frame(i=temp_i,
                                          j=temp_j,
                                          v=2)
        triplet_edge_ls_node[[temp_i]] = data.frame(i=rep(temp_i,2),
                                          j=c(temp_j1,temp_j2),
                                          v=-1)
        edge_info[[temp_i]] = data.frame(edge_id=n,
                                         edge_score=temp_edge$edge_massdif_score
                                         )
        temp_i = temp_i+1
        temp_j = temp_j+1
      }
    } 
    #which(merge_formula$formula==temp_formula_2&merge_formula$id==node_2)
  }
  triplet_edge_ls_edge_sum = bind_rows(triplet_edge_ls_edge)
  triplet_edge_ls_edge_sum$i=triplet_edge_ls_edge_sum$i+num_lib_nodes+num_unknown_nodes
  triplet_edge_ls_edge_sum$j=triplet_edge_ls_edge_sum$j+nrow(pred_formula)
  triplet_edge_ls_node_sum = bind_rows(triplet_edge_ls_node)
  triplet_edge_ls_node_sum$i=triplet_edge_ls_node_sum$i+num_lib_nodes+num_unknown_nodes
  edge_info_sum = bind_rows(edge_info)
  
  
  time = Sys.time()-time
  print(time)
  
  
  

  
  #Generate sparse matrix on left hand side
  triplet_df = rbind(triplet_library_node,
                     triplet_unknown_node,
                     triplet_edge_ls_edge_sum,
                     triplet_edge_ls_node_sum
                     )
  write.csv(triplet_df,"triplet_df.csv",row.names = F)
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
  obj <- c(rep(0,nrow(lib_formula)), rep(-0.05, nrow(unknown_formula)), edge_info_sum$edge_score)
  dir <- c(rep("==",num_lib_nodes+num_unknown_nodes), rep("<=", nrow(edge_info_sum)))
  rhs = c(rep(1,num_lib_nodes+num_unknown_nodes),rep(0,nrow(edge_info_sum)))
  types <- c(rep("B",nrow(lib_formula)+nrow(unknown_formula)+nrow(edge_info_sum)))
  max <- TRUE
}  
  
#Solver
ILP_result = Rglpk_solve_LP(obj, mat, dir, rhs, types = types, max = max)


time = Sys.time()
print(time)
