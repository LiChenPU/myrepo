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
  library(enviPat)
}



time = Sys.time()

#Read data
{
  setwd("C:/Users/Li Chen/Desktop/Github local/myrepo/Network_ID/yeast neg test set")
  getwd()
  raw_node_list=read_csv("merge_node_list.csv")
  raw_edge_list=read_csv("merge_edge_list.csv")
  #raw_pred_formula=read_csv("pred_formula_prune.csv")
  raw_pred_formula=read_csv("All_formula_predict.csv")
  lin_result = read.csv("yeast neg.csv", stringsAsFactors = F)
  
  
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}


#Function handle chemical formula
{
  My_mergefrom=function (formula1, formula2)
  {
    formula2 <- gsub("D", "[2]H", formula2)
    ende2 <- nchar(formula2)
    element2 <- c()
    number2 <- c()
    j <- c(1)
    while (j <= ende2) {
      #browser()
      if (substr(formula2, j, j) == c("[")) {
        b <- j
        while (any(substr(formula2, j, j) == c("]")) != 
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula2, j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element2 <- c(element2, substr(formula2, b, m))
      }
      if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                          "4", "5", "6", "7", "8", "9")) != TRUE) {
        k <- j
        while (any(substr(formula2, j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element2 <- c(element2, substr(formula2, k, m))
      }
      if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                          "4", "5", "6", "7", "8", "9")) == TRUE) {
        k <- j
        while (any(substr(formula2, j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) == 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number2 <- c(number2, as.numeric(substr(formula2, 
                                                k, m)))
      }
      j <- j + 1
    }
    formulas <- c()
    for (i in 1:length(formula1)) {
      #i=1
      formula1[i] <- gsub("D", "[2]H", formula1[i])
      ende1 <- nchar(formula1[i])
      element1 <- c()
      number1 <- c()
      j <- c(1)
      while (j <= ende1) {
        if (substr(formula1[i], j, j) == c("[")) {
          b <- j
          while (any(substr(formula1[i], j, j) == c("]")) != 
                 TRUE) {
            j <- c(j + 1)
          }
          k <- j
          while (any(substr(formula1[i], j, j) == c("0", 
                                                    "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          element1 <- c(element1, substr(formula1[i], 
                                         b, m))
        }
        if (any(substr(formula1[i], j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) != 
            TRUE) {
          k <- j
          while (any(substr(formula1[i], j, j) == c("0", 
                                                    "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          element1 <- c(element1, substr(formula1[i], 
                                         k, m))
        }
        if (any(substr(formula1[i], j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) == 
            TRUE) {
          k <- j
          while (any(substr(formula1[i], j, j) == c("0", 
                                                    "1", "2", "3", "4", "5", "6", "7", "8", "9")) == 
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          number1 <- c(number1, as.numeric(substr(formula1[i], 
                                                  k, m)))
        }
        j <- j + 1
      }
      both <- unique(c(element1, element2))
      both = both[order(both)]
      counts <- c()
      for (i in 1:length(both)) {
        if (any(element1 == both[i])) {
          it1 <- c(number1[element1 == both[i]])
        }
        else {
          it1 <- c(0)
        }
        if (any(element2 == both[i])) {
          it2 <- c(number2[element2 == both[i]])
        }
        else {
          it2 <- c(0)
        }
        counts <- c(counts, it1 + it2)
      }
      formula_all <- ""
      for (i in 1:length(both)) {
        if(counts[i]==0){next}
        formula_all <- paste(formula_all, both[i], counts[i], 
                             sep = "")
      }
      formulas <- c(formulas, formula_all)
    }
    return(formulas)
  }
}

#Clean up
{
  pred_formula = raw_pred_formula[!grepl("-",raw_pred_formula$formula),]
  pred_formula = pred_formula[!duplicated(pred_formula[,c("id","formula")]),]
  
  #Select node that has only 1 formula as library nodes, 
  #For future applicatoin, should select nodes where degree==0
  lib_nodes = raw_node_list[raw_node_list$category==0,]
  lib_nodes_cutoff = nrow(raw_node_list)-nrow(lib_nodes)
  unknown_nodes = raw_node_list[raw_node_list$category!=0,]
  unknown_nodes = unknown_nodes[!is.na(unknown_nodes$Predict_formula),]
  num_unknown_nodes = nrow(unknown_nodes)
  
  lib_formula = pred_formula[pred_formula$id %in% lib_nodes$ID,]
  lib_formula = lib_formula[lib_formula$steps==0,]
  unknown_formula = pred_formula[pred_formula$id %in% unknown_nodes$ID,]
  
  merge_formula = rbind(unknown_formula,lib_formula)
  merge_formula["ilp_index"]=1:nrow(merge_formula)
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
read_from_csv = F
if(!read_from_csv)
{
  #Unknown nodes
  triplet_unknown_nodes_ls = list()
  temp_j=1
  
  for(n in 1:nrow(unknown_nodes)){
    temp_i = n
    temp_unknown_formula = unknown_formula[unknown_formula$id==unknown_nodes$ID[n],]
    if(nrow(temp_unknown_formula)==0){next()}
    triplet_unknown_nodes_ls[[n]]=list(i=rep(temp_i, nrow(temp_unknown_formula)), 
                                       j=temp_j:(temp_j+nrow(temp_unknown_formula)-1), 
                                       v=rep(1,nrow(temp_unknown_formula)))
    temp_j = temp_j+nrow(temp_unknown_formula)
  }
  triplet_unknown_node = bind_rows(triplet_unknown_nodes_ls)
  rm(triplet_unknown_nodes_ls)
  
  
  #Edge list
  
  temp_i=temp_j=1
  triplet_edge_ls_edge=triplet_edge_ls_node=list()
  edge_info = list()
  timer=Sys.time()
  n=1
  for(n in 1:nrow(edge_list)){
    #for(n in 1:10000){
    if(n%%1000==0){
      print(paste("n=",n,"elapsed="))
      print(Sys.time()-timer)
    }
    
    temp_edge = edge_list[n,]
    node_1 = temp_edge$node1
    node_2 = temp_edge$node2
    formula_1 = pred_formula_ls[[node_1]]
    formula_2 = pred_formula_ls[[node_2]]
    temp_fg = temp_edge$linktype
    temp_formula = unique(formula_1$formula[2])
    for(temp_formula in unique(formula_1$formula)){
      #Assuming formula in node_1 is always smaller than node_2
      if(temp_fg=="Same"){
        temp_formula_2=temp_formula
      }
      else{
        #temp_formula_2 = formula_manipulate(temp_formula, temp_fg, +1)
        
        temp_formula_2 = My_mergefrom(temp_formula, temp_fg)
      }
      #Write triplet for edge and corresponding 2 nodes
      if(temp_formula_2 %in% formula_2$formula){
        temp_j1 = formula_1$ilp_index[which(formula_1$formula==temp_formula )]
        temp_j2 = formula_2$ilp_index[which(formula_2$formula==temp_formula_2 )]
        #temp_j1 = which(merge_formula$formula==temp_formula&merge_formula$id==node_1)
        #temp_j2 = which(merge_formula$formula==temp_formula_2&merge_formula$id==node_2)
        
        #if one node is library node,
        if(node_1>lib_nodes_cutoff){
          triplet_edge_ls_edge[[temp_i]] = list(i=temp_i,
                                                j=temp_j,
                                                v=1)
          triplet_edge_ls_node[[temp_i]] = list(i=temp_i,
                                                j=temp_j2,
                                                v=-1)
          edge_info[[temp_i]] = list(edge_id=n,
                                     edge_score=temp_edge$edge_massdif_score,
                                     formula1 = temp_formula,
                                     formula2 = temp_formula_2,
                                     ilp_index1 = temp_j1,
                                     ilp_index2 = temp_j2
          )
        }
        if(node_2>lib_nodes_cutoff){
          triplet_edge_ls_edge[[temp_i]] = list(i=temp_i,
                                                j=temp_j,
                                                v=1)
          triplet_edge_ls_node[[temp_i]] = list(i=temp_i,
                                                j=temp_j1,
                                                v=-1)
          edge_info[[temp_i]] = list(edge_id=n,
                                     edge_score=temp_edge$edge_massdif_score,
                                     formula1 = temp_formula,
                                     formula2 = temp_formula_2,
                                     ilp_index1 = temp_j1,
                                     ilp_index2 = temp_j2
          )
        }
        
        #if both nodes are unknown nodes
        if(node_1<=lib_nodes_cutoff&node_2<=lib_nodes_cutoff){
          triplet_edge_ls_edge[[temp_i]] = list(i=temp_i,
                                                j=temp_j,
                                                v=2)
          triplet_edge_ls_node[[temp_i]] = list(i=c(temp_i,temp_i),
                                                j=c(temp_j1,temp_j2),
                                                v=c(-1,-1))
          edge_info[[temp_i]] = list(edge_id=n,
                                     edge_score=temp_edge$edge_massdif_score,
                                     formula1 = temp_formula,
                                     formula2 = temp_formula_2,
                                     ilp_index1 = temp_j1,
                                     ilp_index2 = temp_j2
          )
        }
        temp_i = temp_i+1
        temp_j = temp_j+1
      }
    } 
  }
  
  
  
  triplet_edge_ls_edge_sum = bind_rows(triplet_edge_ls_edge)
  triplet_edge_ls_edge_sum$i=triplet_edge_ls_edge_sum$i+num_unknown_nodes
  triplet_edge_ls_edge_sum$j=triplet_edge_ls_edge_sum$j+nrow(unknown_formula)
  triplet_edge_ls_node_sum = bind_rows(triplet_edge_ls_node)
  triplet_edge_ls_node_sum$i=triplet_edge_ls_node_sum$i+num_unknown_nodes
  
  #Generate sparse matrix on left hand side
  triplet_df = rbind(
    triplet_unknown_node, 
    triplet_edge_ls_edge_sum,
    triplet_edge_ls_node_sum
  )
  
  
  ##Objective parameter 
  {
    
    edge_info_sum = bind_rows(edge_info)
    edge_info_sum["edge_ilp_id"]=1:nrow(edge_info_sum)
    test = edge_info_sum
    
    
    
    test1 = test[test$ilp_index2>nrow(unknown_formula)&
                   test$ilp_index1<=nrow(unknown_formula),]
    test2 = test[test$ilp_index1>nrow(unknown_formula)&
                   test$ilp_index2<=nrow(unknown_formula),]
    colnames(test2)=sub(1,3,colnames(test2))
    colnames(test2)=sub(2,1,colnames(test2))
    colnames(test2)=sub(3,2,colnames(test2))
    
    test1 = merge(test1,test2,all=T)
    
    test1 = test1[duplicated(test1[,c("formula1","ilp_index1")]) | 
                    duplicated(test1[,c("formula1","ilp_index1")], fromLast=TRUE),]
    test1 = test1[order(test1$formula1,test1$edge_score,decreasing = T),]
    test1$edge_ilp_id[duplicated(test1[,c("ilp_index1","formula1")])]
    edge_info_sum$edge_score[test1$edge_ilp_id[duplicated(test1[,c("ilp_index1","formula1")])]]=0
    
    # test2 = test[test$ilp_index1>nrow(unknown_formula)&
    #                test$ilp_index2<=nrow(unknown_formula),]
    # test2 = test2[duplicated(test2[,c("formula2","ilp_index2")]) | 
    #                 duplicated(test2[,c("formula2","ilp_index2")], fromLast=TRUE),]
    # test2 = test2[order(test2$formula1,test2$edge_score,decreasing = T),]
    # edge_info_sum$edge_score[test2$edge_ilp_id[duplicated(test2[,c("ilp_index2","formula2")])]]=0
    
    
    
    # 
    # 
    # test6 = test1[duplicated(test1[,c("edge_id")]) | 
    #                 duplicated(test1[,c("edge_id")], fromLast=TRUE),]
  }
  
  write_csv(triplet_df,"triplet_df.csv")
  write_csv(edge_info_sum,"edge_info_sum.csv")
  mat = simple_triplet_matrix(i=triplet_df$i,
                              j=triplet_df$j,
                              v=triplet_df$v)
}else
{
  #triplet_df = read_csv("triplet_df.csv")
  triplet_df = read.csv("triplet_df.csv")
  edge_info_sum = read_csv("edge_info_sum.csv")
  mat = simple_triplet_matrix(i=triplet_df$i,
                              j=triplet_df$j,
                              v=triplet_df$v)
  
}



print(Sys.time()-time)



#CPLEX solver parameter
{

  
  nc <- max(mat$j)
  obj <- c(rep(0, nrow(unknown_formula)), edge_info_sum$edge_score)
  lb <- rep(0, nc)
  ub <- rep(1, nc)
  ctype <- rep("B",nc)
  
  nr <- max(mat$i)
  rhs = c(rep(1,nrow(unknown_nodes)),rep(0,nrow(edge_info_sum)))
  sense <- c(rep("L",nrow(unknown_nodes)), rep("L", nrow(edge_info_sum)))
  
  # cn <- c("x1", "x2", "x3")
  # rn <- c("q1", "q2", "q3")
  
  
  triplet_df=triplet_df[with(triplet_df,order(j)),]
  cnt=as.vector(table(triplet_df$j))
  beg=vector()
  beg[1]=0
  for(i in 2:length(cnt)){beg[i]=beg[i-1]+cnt[i-1]}
  ind=triplet_df$i-1
  val = triplet_df$v
}


##Permutation to edge score 

solution_ls = list()
for(i in 1:1){
  
{
  edge_score_permutated = edge_info_sum$edge_score * (1+ rnorm(length(edge_info_sum$edge_score),
                                                              mean = 0,
                                                              sd = 0.1))
  #obj <- c(rep(0, nrow(unknown_formula)), edge_score_permutated)
}
  
{
  env <- openEnvCPLEX()
  prob <- initProbCPLEX(env)
  copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj, rhs, sense,
                    beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
  copyColTypeCPLEX(env, prob, ctype)
  
  #primoptCPLEX(env, prob)
  #lpoptCPLEX(env, prob)
  #dualoptCPLEX(env, prob)
  mipoptCPLEX(env, prob)
  
  result_solution=solutionCPLEX(env, prob)
  #writeProbCPLEX(env, prob, "prob.lp")
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
}
  
  print(i)
  solution_ls[[i]]=result_solution
}


test_x = as.data.frame(lapply(solution_ls, function(x){return (x[[3]])}))
test_x_rowmean = rowMeans(test_x)


unknown_formula["ILP_result"] = test_x_rowmean[1:nrow(unknown_formula)]
test46=unknown_formula[unknown_formula$ILP_result!=0&unknown_formula$ILP_result!=1,]

##Evaluation
{
  unknown_formula["ILP_result"] = result_solution$x[1:nrow(unknown_formula)]
  edge_info_sum["ILP_result"] = result_solution$x[(nrow(unknown_formula)+1):length(result_solution$x)]
  table(result_solution$x[(nrow(unknown_formula)+1):length(result_solution$x)])
  
  table(result_solution$x)
  #merge_formula$id[merge_formula$ilp_index==584]
  
  unknown_formula_CPLEX = unknown_formula[unknown_formula$ILP_result==1,]
  unknown_node_CPLEX = merge(unknown_nodes,unknown_formula_CPLEX,by.x = "ID", by.y = "id",all=T)
  
  edge_info_CPLEX = edge_info_sum[edge_info_sum$ILP_result==1,]
  
  
  
  lin_result = lin_result[lin_result$feature=="Metabolite"|lin_result$feature=="Adduct"|lin_result$feature=="Fragment"|
                            lin_result$feature=="[]",]
  lin_result = lin_result[order(lin_result$mz),]
  
  data("isotopes")
  lin_result["formula_check"]=check_chemform(isotopes,lin_result$formula)$new_formula
  lin_result["ilp_id"]=1:nrow(lin_result)
  
  
  #test_formula = unknown_node_CPLEX[unknown_node_CPLEX$formula%in%lin_result$formula_check,]
  #merge_Lin_ILP = merge(lin_result, test_formula, by.x="formula_check",by.y="formula",all=T)
  #merge_Lin_ILP["mz_diff"]=merge_Lin_ILP$mz.x-merge_Lin_ILP$mz.y
  #merge_Lin_ILP["rt_diff"]=merge_Lin_ILP$RT.x-merge_Lin_ILP$RT.y
  
    
    
  merge_Lin_ILP = merge(lin_result, unknown_node_CPLEX, by.x="ilp_id",by.y="ID",all=T)
  merge_Lin_ILP_metabolite = merge_Lin_ILP[merge_Lin_ILP$feature=="Metabolite",]
  merge_Lin_ILP_metabolite = merge_Lin_ILP_metabolite[!is.na(merge_Lin_ILP_metabolite$formula.y),]
  merge_Lin_ILP_metabolite_diff = merge_Lin_ILP_metabolite[merge_Lin_ILP_metabolite$formula_check!=
                                                             merge_Lin_ILP_metabolite$formula.y,]
  withoutILP=withILP=0
  for(i in 1:nrow(merge_Lin_ILP_metabolite)){
    if(merge_Lin_ILP_metabolite$formula_check[i]==merge_Lin_ILP_metabolite$Predict_formula[i])
    {withoutILP=withoutILP+1}
    if(merge_Lin_ILP_metabolite$formula_check[i]==merge_Lin_ILP_metabolite$formula.y[i])
    {withILP=withILP+1}
  }
  withoutILP = sum(merge_Lin_ILP_metabolite$formula_check==merge_Lin_ILP_metabolite$Predict_formula)
  
  unknown_node_CPLEX_diff = unknown_node_CPLEX[unknown_node_CPLEX$ID %in% merge_Lin_ILP_metabolite_diff$ilp_id,]
  
  edge_info_sum_debug = edge_info_sum[unique(c(which(edge_info_sum$formula1=="C16H12O8"),
                                               which(edge_info_sum$formula2=="C16H12O8"))),]
  unknown_formula[unknown_formula$id==merge_formula$id[which(merge_formula$ilp_index==3584)],]
  edge_info_sum_debug2 = edge_info_sum[unique(c(which(edge_info_sum$formula1=="C20H25N1S1"),
                                               which(edge_info_sum$formula2=="C20H25N1S1"))),]
  edge_info_sum_debug3 = edge_info_sum[unique(c(which(edge_info_sum$formula1=="C11H14O6"),
                                                which(edge_info_sum$formula2=="C11H14O6"))),]
  
  ## C5H7NOS 10 ppm off
  
  
  
  edge_info_sum_debug = edge_info_sum[unique(c(which(edge_info_sum$formula1=="C19H18O13"),
                                               which(edge_info_sum$formula2=="C19H18O13"))),]
  # raw_merge = merge(raw, unknown_formula_CPLEX, by.x = "ID", by.y = "id",all=T)
  # raw_merge_diff = raw_merge[raw_merge$formula.x!=raw_merge$formula.y,]
}


{
  param <- getChgParmCPLEX(env)
  solnInfoCPLEX(env, prob)
  
  getStatCPLEX(env, prob)
  
}


# 
# #Other parameters
# {
#   #should have number of columns, i.e. j
#   obj <- c(rep(0, nrow(unknown_formula)), edge_info_sum$edge_score)
#   types <- c(rep("B",nrow(unknown_formula)+nrow(edge_info_sum)))
#   
#   #should have number of rows, i.e. i
#   dir <- c(rep("==",num_unknown_nodes), rep("<=", nrow(edge_info_sum)))
#   rhs = c(rep(1,num_unknown_nodes),rep(0,nrow(edge_info_sum)))
#   max <- TRUE
# }  
# 
# 
# #R glpk Solver
# ILP_result = Rglpk_solve_LP(obj, mat, dir, rhs, types = types, max = max,
#                             control = list(tm_limit=10*60*1000,canonicalize_status=F,verbose=T))
# para = data.frame(ILP_result$solution)
# unknown_formula["ILP_result"] = para$ILP_result.solution[1:nrow(unknown_formula)]
# edge_info_sum["ILP_result"] = para$ILP_result.solution[(nrow(unknown_formula)+1):nrow(para)]

# test=data.frame(ILP_result$solution)
# test=cbind(test,ILP_result$solution)
# test2=test[test$ILP_result.solution!=test$`ILP_result$solution`,]
