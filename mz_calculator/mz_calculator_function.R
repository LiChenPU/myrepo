

library(lc8)
library(enviPat)
library(dplyr)
library(fitdistrplus)
library(slam)
library(cplexAPI)

## read_library ####
read_library = function(library_file){
  data(isotopes)
  hmdb_lib = read.csv(library_file, stringsAsFactors = F)
  hmdb_lib$MF = check_chemform(isotopes, hmdb_lib$MF)$new_formula
  hmdb_lib$MF = sapply(hmdb_lib$MF, my_calculate_formula,"C1")
  hmdb_lib$MF = sapply(hmdb_lib$MF, my_calculate_formula,"C1",-1)
  hmdb_lib$Exact_Mass = formula_mz(hmdb_lib$MF)
  hmdb_lib["rdbe"]=formula_rdbe(hmdb_lib$MF)
  return(hmdb_lib)
}

## Read_rule_table - for Connect_rules ####
Read_rule_table = function(rule_table_file, extend_rule = F){
  data("isotopes")
  Connect_rules = read.csv(rule_table_file,stringsAsFactors = F)
  for(i in 1: nrow(Connect_rules)){
    if(Connect_rules$Formula[i]==""){next}
    Connect_rules$Formula[i] = check_chemform(isotopes,Connect_rules$Formula[i])$new_formula
    Connect_rules$Formula[i] = my_calculate_formula(Connect_rules$Formula[i], "C1",Is_valid = F)
    Connect_rules$Formula[i] = my_calculate_formula(Connect_rules$Formula[i], "C1", -1 ,Is_valid = F)
    Connect_rules$mass[i] = formula_mz(Connect_rules$Formula[i])
  }
  
  if(extend_rule){
    extend_rules = list()
    i=3
    for(i in 1: nrow(Connect_rules)){
      if(Connect_rules$allow_rep[i] <= 1){next}
      for(j in 2:Connect_rules$allow_rep[i]){
        temp_rule = Connect_rules[i,]
        temp_rule$Symbol = paste(temp_rule$Symbol, "x", j, sep="")
        temp_rule$Formula = my_calculate_formula(temp_rule$Formula, temp_rule$Formula, 
                                                 sign = j-1,Is_valid = F)
        temp_rule$mass = temp_rule$mass * j
        temp_rule$rdbe = temp_rule$rdbe * j
        
        extend_rules[[length(extend_rules)+1]] = temp_rule
      }
      
    }
    
    Connect_rules = rbind(Connect_rules, bind_rows(extend_rules))
    
    
    Connect_rules = Connect_rules %>%
      arrange(mass) %>%
      dplyr::select(-allow_rep)
  }

   
  return(Connect_rules)
}

## Peak_cleanup - Clean up duplicate peaks from peak picking ####
Peak_cleanup = function(Mset)
{
  raw = Mset$Raw_data
  
  # colnames(raw)[colnames(raw)=="ID"] = "oldID"
  raw["ID"] = 1:nrow(raw)
  
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ion_mode = Mset$Global_parameter$mode
  raw$mz = raw$mz - (H_mass-e_mass)*ion_mode
  
 
  return(raw)
}

## Form_node_list ####
Form_node_list = function(Mset)
{
  NodeSet = list()
  NodeSet[["Expe"]] = Mset$Data[,c("ID","mz")]
  NodeSet$Expe["formula"]=NA
  NodeSet$Expe["category"]=1
  NodeSet$Expe["compound_name"]=NA
  NodeSet$Expe["rdbe"]=NA
  colnames(NodeSet$Expe) = c("ID","mz","MF", "category","compound_name","rdbe")
  
  NodeSet[["Library"]] = Mset$Library
  colnames(NodeSet$Library) = c("ID","compound_name","MF","mz","category","rdbe")
  NodeSet$Library$ID = 1:nrow(NodeSet$Library)+nrow(NodeSet$Expe)
  NodeSet$Library$category=0
  
  merge_node_list = rbind(NodeSet$Expe,NodeSet$Library)
  
  
  return(merge_node_list)
}

## Edge_Connect_rules ####
Edge_Connect_rules = function(Mset, mass_abs = 0.001, mass_ppm = 5/10^6)
{
  
  
  merge_node_list = Mset$NodeSet[Mset$NodeSet$category!=-1,]
  merge_node_list = merge_node_list[with(merge_node_list, order(mz)),]
  merge_nrow = nrow(merge_node_list)
  
  timer=Sys.time()
  {
    edge_ls = list()
    temp_mz_list=merge_node_list$mz
    
    for (k in 1:nrow(Mset$Connect_rules)){
      # print(k)
      temp_fg=Mset$Connect_rules$mass[k]
      temp_direction = Mset$Connect_rules$direction[k]
      temp_rdbe = Mset$Connect_rules$rdbe[k]
      i=j=1
      temp_edge_list = data.frame(node1=as.numeric(), 
                                  node2=as.numeric(), 
                                  linktype=as.numeric(), 
                                  mass_dif=as.numeric(), 
                                  direction = as.numeric(),
                                  rdbe = as.numeric())
      while(i<=merge_nrow){
        temp_ms=0
        mass_tol = max(temp_mz_list[i]*mass_ppm,mass_abs)
        while(1){
          j=j+1
          if(j>merge_nrow){break}
          temp_ms = temp_mz_list[j]-temp_mz_list[i]
          
          if(abs(temp_ms-temp_fg)<mass_tol){
            temp_edge_list[nrow(temp_edge_list)+1,]=list(merge_node_list$ID[i], 
                                                         merge_node_list$ID[j], 
                                                         k, 
                                                         (temp_ms-temp_fg)/temp_mz_list[j]*1E6, 
                                                         temp_direction,
                                                         temp_rdbe
            )
          }
          if((temp_ms-temp_fg)>mass_tol){break}
        }
        i=i+1
        
        while(j>i){
          j=j-1
          temp_ms = temp_mz_list[j]-temp_mz_list[i]
          if((temp_ms-temp_fg)<(-mass_tol)){break}
        }
      }
      edge_ls[[k]]=temp_edge_list
      
    }
  }
  print("Finish Connect_rules.")
  # print(Sys.time()-timer)
  edge_list = bind_rows(edge_ls)
  
  edge_list$linktype=Mset$Connect_rules$Formula[edge_list$linktype]
  
  edge_list_sub = subset(edge_list, 
                         (edge_list$node1<=nrow(Mset$Data)|
                            edge_list$node2<=nrow(Mset$Data))&
                           edge_list$node1!=edge_list$node2
  )
  
  edge_list_sub["category"]=1
  # write_csv(edge_list_sub,"edge_list_sub.txt")
    
  
  return(edge_list_sub)
}

## Merge_edgeset ####
Merge_edgeset = function(EdgeSet){
  edge_merge = EdgeSet$Connect_rules
  edge_merge["edge_id"]=1:nrow(edge_merge)
  
  
  return(edge_merge)
}


## Edge_score
Edge_score = function(Biotransform, fix_distribution_sigma = F ,plot_graph = F){
  #Biotransform = EdgeSet$Biotransform
  if(fix_distribution_sigma){
    temp_sigma = fix_distribution_sigma
  } else {
    if(nrow(Biotransform)>10000){
      edge_mzdif_FIT <- fitdist(as.numeric(Biotransform$mass_dif[base::sample(nrow(Biotransform),10000)]), "norm")    
    } else {
      edge_mzdif_FIT <- fitdist(as.numeric(Biotransform$mass_dif), "norm")    
    }
    
    if(plot_graph){
      plot(edge_mzdif_FIT)
      print(summary(edge_mzdif_FIT))
    }
    temp_sigma = edge_mzdif_FIT$estimate[2]
  }
 
  
  Biotransform["edge_massdif_score"]=dnorm(Biotransform$mass_dif, 0, temp_sigma)
  Biotransform["edge_massdif_score"]=Biotransform["edge_massdif_score"]/max(Biotransform["edge_massdif_score"])
  return(Biotransform)
}

## Network_prediction ####
Network_prediction = function(Mset, 
                              edge_biotransform, 
                              biotransform_step = 10,
                              propagation_score_threshold = 0.25,
                              top_n = 50,
                              print_step = F
)
{
  
  mnl=Mset$NodeSet
  # edge_biotransform = EdgeSet$Connect_rules
  # top_formula_n=1
  # edge_list_sub = EdgeSet$Merge
  #Initialize predict_formula from HMDB known formula
  {
    data_str = data.frame(id=as.numeric(),
                          formula=as.character(), 
                          steps=as.numeric(), 
                          parent=as.numeric(), 
                          is_metabolite = as.logical(),
                          score=as.numeric(),
                          rdbe=as.numeric(),
                          stringsAsFactors = F)
    #sf stands for summary_formula
    
    sf = lapply(1:nrow(mnl),function(i)(data_str))
    
    for(i in mnl$ID[mnl$category==0]){
      Initial_formula =sf[[i]]
      Initial_formula[1,]= list(i,mnl$MF[i],0,0,T,1,mnl$rdbe[i])
      sf[[i]]=Initial_formula
    }
  }
  
  #while loop to Predict formula based on known formula and edgelist 
  nrow_experiment = nrow(Mset$Data)
  step=0
  timer = Sys.time()
  timer=Sys.time()
  while(step <= biotransform_step){
    
    all_nodes_df = bind_rows(sf)
    
    
    # Handle biotransform
    {
      all_nodes_df = bind_rows(sf)
      new_nodes_df = all_nodes_df[all_nodes_df$steps==step
                                  &all_nodes_df$score>propagation_score_threshold,]
      if(print_step){
        print(paste("nrow",nrow(all_nodes_df),"in step",step,"elapsed="))
        print((Sys.time()-timer))
      }
      step = step + 1
      if(nrow(new_nodes_df)==0 | step > biotransform_step){break}
      edge_biotransform_sub = edge_biotransform[edge_biotransform$node1 %in% new_nodes_df$id | 
                                                  edge_biotransform$node2 %in% new_nodes_df$id,]
      
      for(n in 1:nrow(new_nodes_df)){
        temp_new_node = new_nodes_df[n,]
        flag_id = temp_new_node$id
        flag_formula = temp_new_node$formula
        flag_is_metabolite = temp_new_node$is_metabolite
        flag_score = temp_new_node$score
        flag_rdbe = temp_new_node$rdbe
        
        # temp_edge_list=subset(edge_biotransform_sub, edge_biotransform_sub$node1==flag_id |
        #                         edge_biotransform_sub$node2==flag_id)
        flag_id_in_node1_or_node2 = edge_biotransform_sub$node1==flag_id | edge_biotransform_sub$node2==flag_id
        temp_edge_list = edge_biotransform_sub[flag_id_in_node1_or_node2,]
        
        # #If head signal is < defined cutoff, then prevent it from propagating out, but it can still get formula from others.
        # if(flag_id <= nrow_experiment){
        #   if(Mset$Data$mean_inten[flag_id]< 2e4){next}
        # }
        # 
        # #If flag is an isotopic peak, then only look for isotopic peaks
        # if(grepl("\\[",flag_formula)){
        #   temp_edge_list = temp_edge_list[grepl("\\[",temp_edge_list$category),]
        #   if(nrow(temp_edge_list)==0){next}
        # }
        # 
        # #Filter edge direction
        # temp_edge_list = temp_edge_list[(temp_edge_list$node1 == flag_id & temp_edge_list$direction != -1) 
        #                                 |(temp_edge_list$node2 == flag_id & temp_edge_list$direction != 1),]
        if(nrow(temp_edge_list)==0){next}
        
        i=1
        for(i in 1:nrow(temp_edge_list)){
          temp_rdbe = temp_edge_list$rdbe[i]
          #If flag is head
          if(temp_edge_list$node1[i] == flag_id){
            partner_id = temp_edge_list$node2[i]
            if(partner_id>nrow_experiment){next}
            if(grepl("x", temp_rdbe)){
              fold = as.numeric(gsub("x","",temp_rdbe))
              partner_rdbe = flag_rdbe * fold
            } else{
              partner_rdbe = flag_rdbe + as.numeric(temp_rdbe)
            }
            
            temp_fg = temp_edge_list$linktype[i]
            if(temp_fg==""){
              partner_formula = flag_formula
            }else if (grepl("x",temp_fg)){
              fold = as.numeric(gsub("x","",temp_fg))
              partner_formula=my_calculate_formula(flag_formula,flag_formula,fold-1,Is_valid = T)
            }else {
              partner_formula=my_calculate_formula(flag_formula,temp_fg,1,Is_valid = T)
            }
          }
          
          #If flag is tail
          if(temp_edge_list$node2[i] == flag_id){
            partner_id = temp_edge_list$node1[i]
            if(partner_id>nrow_experiment){next}
            
            if(grepl("x", temp_rdbe)){
              fold = as.numeric(gsub("x","",temp_rdbe))
              partner_rdbe = flag_rdbe / fold
            } else{
              partner_rdbe = flag_rdbe - as.numeric(temp_rdbe)
            }
            temp_fg = temp_edge_list$linktype[i]
            if(temp_fg==""){
              partner_formula = flag_formula
            }else if (grepl("x",temp_fg)){
              fold = as.numeric(gsub("x","",temp_fg))
              partner_formula=my_calculate_formula(flag_formula,flag_formula,-(fold-1)/fold,Is_valid = T)
              if(grepl(".", partner_formula)){partner_formula=F}
            }else {
              partner_formula=my_calculate_formula(flag_formula,temp_fg,-1,Is_valid = T)
            }
          }
          
          #function return false if not valid formula
          if(is.logical(partner_formula)){next}
          
          #Criteria to enter new entry into formula list
          #1. If it is a new formula or a new metabolite status, then record
          
          partner_score = flag_score*temp_edge_list$edge_massdif_score[i]
          partner_parent = flag_id
          partner_steps = step
          partner_is_metabolite = T
          
          temp = sf[[partner_id]]
          
          formula_metabolite_status_matched = temp$formula==partner_formula & temp$is_metabolite==partner_is_metabolite
          temp_subset = temp[formula_metabolite_status_matched,]
          # temp_subset=subset(temp, temp$formula==partner_formula & temp$is_metabolite==partner_is_metabolite)
          
          if(nrow(temp_subset)!=0){
            #2. if not much higher scores, then next
            if(partner_score<=(1.2*max(temp_subset$score))){
              next
            }
          }
          
          #New entry
          temp[nrow(temp)+1, ] = list(partner_id, 
                                      partner_formula, 
                                      partner_steps, 
                                      partner_parent,
                                      partner_is_metabolite, 
                                      partner_score,
                                      partner_rdbe)
          temp = temp[with(temp, order(-score)),]
          sf[[partner_id]]=temp
        }
      }
    }
  }
  
  
  # Pruning formula  #
  {
    pred_formula = bind_rows(lapply(sf, head, top_n))
    #pred_formula = pred_formula[!grepl("-",pred_formula$formula),]
    pred_formula = pred_formula[!duplicated(pred_formula[,c("id","formula")]),]
    
    merge_formula = pred_formula
    sf = list()
    for(n in 1: max(merge_formula$id)){
      sf[[n]]=merge_formula[merge_formula$id==n,]
    }
    
    # write_csv(merge_formula,"All_formula_predict.txt")
  }
  
  print("Finish Network_prediction.")
  return(sf)
}



## Prepare_CPLEX parameter ####
Prepare_CPLEX = function(Mset, EdgeSet){
  raw_pred_formula=bind_rows(Mset[["NodeSet_network"]])
  raw_node_list = Mset$NodeSet
  raw_edge_list = EdgeSet$Merge
  
  #Clean up
  {
    pred_formula = raw_pred_formula
    
    #1 is measured peaks; 0 is library; -1 is system adduct
    lib_nodes = raw_node_list[raw_node_list$category!=1,]
    lib_nodes_cutoff = nrow(Mset$Data)
    unknown_nodes = raw_node_list[raw_node_list$category==1,]
    unknown_nodes = unknown_nodes[unknown_nodes$ID %in% unique(pred_formula$id),]
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
    merge_formula_id=merge_formula$id
    for(n in 1: max(merge_formula$id)){
      pred_formula_ls[[n]]=merge_formula[merge_formula_id==n,]
    }
  }
  
  
  
  ##Core codes
  
  #Construct constraint matrix 

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
    gc()
    temp_i=temp_j=1
    triplet_edge_ls_edge=triplet_edge_ls_node=list()
    edge_info = list()
    timer=Sys.time()
    n=1
    for(n in 1:nrow(edge_list)){
      # for(n in 1:20000){
      if(n%%10000==0){
        # print(paste("n=",n,"elapsed="))
        # print(Sys.time()-timer)
        
      }
      
      
      temp_edge = edge_list[n,]
      node_1 = temp_edge$node1
      node_2 = temp_edge$node2
      formula_1 = pred_formula_ls[[node_1]]
      formula_2 = pred_formula_ls[[node_2]]
      temp_fg = temp_edge$linktype
      
      
      #temp_formula = formula_1$formula[2]
      for(temp_formula in unique(formula_1$formula)){
        temp_score = temp_edge$edge_massdif_score
        
        
        #Assuming formula in node_1 is always smaller than node_2
        if(temp_fg==""){
          temp_formula_2 = temp_formula
        }else if (grepl("x",temp_fg)){
          fold = as.numeric(gsub("x","",temp_fg))
          temp_formula_2=my_calculate_formula(temp_formula,temp_formula,fold-1,Is_valid = T)
        }else {
          temp_formula_2=my_calculate_formula(temp_formula,temp_fg,1,Is_valid = T)
        }
        
        #Write triplet for edge and corresponding 2 nodes
        if(temp_formula_2 %in% formula_2$formula){
          temp_j1 = formula_1$ilp_index[which(formula_1$formula==temp_formula )]
          temp_j2 = formula_2$ilp_index[which(formula_2$formula==temp_formula_2 )]
          
          
          #if one node is library node,
          if(node_1>lib_nodes_cutoff){
            triplet_edge_ls_edge[[temp_i]] = list(i=temp_i,
                                                  j=temp_j,
                                                  v=1)
            triplet_edge_ls_node[[temp_i]] = list(i=temp_i,
                                                  j=temp_j2,
                                                  v=-1)
            edge_info[[temp_i]] = list(edge_id=temp_edge$edge_id,
                                       edge_score= temp_score,
                                       formula1 = temp_formula,
                                       formula2 = temp_formula_2,
                                       ILP_id1 = temp_j1,
                                       ILP_id2 = temp_j2
            )
          }
          if(node_2>lib_nodes_cutoff){
            triplet_edge_ls_edge[[temp_i]] = list(i=temp_i,
                                                  j=temp_j,
                                                  v=1)
            triplet_edge_ls_node[[temp_i]] = list(i=temp_i,
                                                  j=temp_j1,
                                                  v=-1)
            edge_info[[temp_i]] = list(edge_id=temp_edge$edge_id,
                                       edge_score=temp_score,
                                       formula1 = temp_formula,
                                       formula2 = temp_formula_2,
                                       ILP_id1 = temp_j1,
                                       ILP_id2 = temp_j2
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
            edge_info[[temp_i]] = list(edge_id=temp_edge$edge_id,
                                       edge_score=temp_score,
                                       formula1 = temp_formula,
                                       formula2 = temp_formula_2,
                                       ILP_id1 = temp_j1,
                                       ILP_id2 = temp_j2
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
      #edge_info_sum = CPLEXset$data$edge_info_sum
      gc()
      edge_info_sum = bind_rows(edge_info)
      
      
      edge_info_sum["edge_ilp_id"]=1:nrow(edge_info_sum)
      test = edge_info_sum
      test1 = test[test$ILP_id2>nrow(unknown_formula)&
                     test$ILP_id1<=nrow(unknown_formula),]
      test2 = test[test$ILP_id1>nrow(unknown_formula)&
                     test$ILP_id2<=nrow(unknown_formula),]
      colnames(test2)=sub(1,3,colnames(test2))
      colnames(test2)=sub(2,1,colnames(test2))
      colnames(test2)=sub(3,2,colnames(test2))
      
      test1 = merge(test1,test2,all=T)
      
      test1 = test1[duplicated(test1[,c("formula1","ILP_id1")]) | 
                      duplicated(test1[,c("formula1","ILP_id1")], fromLast=TRUE),]
      test1 = test1[order(test1$formula1,test1$edge_score,decreasing = T),]
      #test1$edge_ilp_id[duplicated(test1[,c("ILP_id1","formula1")])
      edge_info_sum$edge_score[test1$edge_ilp_id[duplicated(test1[,c("ILP_id1","formula1")])]]=1e-10
      
    }
    
    
    
    # write_csv(triplet_df,"triplet_df.txt")
    # write_csv(edge_info_sum,"edge_info_sum.txt")
    mat = simple_triplet_matrix(i=triplet_df$i,
                                j=triplet_df$j,
                                v=triplet_df$v)

  
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
    
    triplet_df=triplet_df[with(triplet_df,order(j)),]
    cnt=as.vector(table(triplet_df$j))
    beg=vector()
    beg[1]=0
    for(i in 2:length(cnt)){beg[i]=beg[i-1]+cnt[i-1]}
    ind=triplet_df$i-1
    val = triplet_df$v
  }
  CPX_MAX = -1
  CPLEX_para = list(nc = nc,
                    nr = nr,
                    CPX_MAX = CPX_MAX,
                    obj = obj,
                    rhs = rhs,
                    sense = sense,
                    beg = beg,
                    cnt = cnt,
                    ind = ind, 
                    val = val,
                    lb = lb,
                    ub = ub,
                    ctype = ctype
  )
  CPLEX_data = list(edge_info_sum = edge_info_sum,
                    pred_formula_ls = pred_formula_ls,
                    unknown_nodes = unknown_nodes,
                    unknown_formula = unknown_formula
  )
  print("Finish CPLEXset.")
  return(CPLEX = list(data = CPLEX_data,
                      para = CPLEX_para)
  )
}
### Score_formula ####
Score_formula = function(Mset, CPLEXset, rdbe=T, step_score=T)
{
  unknown_formula = CPLEXset$data$unknown_formula
  
  #when measured and calculated mass differ, score based on normal distirbution with mean=0 and sd=1e-3
  unknown_formula["msr_mass"] = Mset$Data$mz[unknown_formula$id]
  unknown_formula["cal_mass"] = formula_mz(unknown_formula$formula)
  
  unknown_formula["msr_cal_mass_dif"] = unknown_formula["msr_mass"]-unknown_formula["cal_mass"]
  # edge_mzdif_FIT <- fitdist(unknown_formula$msr_cal_mass_dif*1000, "norm")    
  # summary(edge_mzdif_FIT)
  unknown_formula["Mass_score"] = dnorm(unknown_formula$msr_cal_mass_dif, 0, 1e-3)/dnorm(0, 0, 1e-3)
  
  #when rdbe < -1, penalty to each rdbe less
  #unknown_formula["rdbe"] = formula_rdbe(unknown_formula$formula)
  if(rdbe==T){
    unknown_formula["rdbe_score"] = sapply((0.2*(unknown_formula$rdbe)), min, 0)
  } else{
    unknown_formula["rdbe_score"]=0
  }
  
  #when step is large, the likelihood of the formula is true decrease from its network score
  #Penalty happens when step > 5 on the existing score
  if(step_score==T){
    unknown_formula["step_score"] = sapply(-0.1*(unknown_formula$steps-3), min, 0) 
  } else{
    unknown_formula["step_score"]=0
  }
  
  
  # the mass score x step score evaluate from mass perspective how likely the formula fits the peak
  # the rdbe score penalizes unsaturation below -1
  #Each node should be non-positive, to avoid node formula without edge connection
  unknown_formula["cplex_score"] = log10(unknown_formula["Mass_score"]) +
    log10(unknown_formula["score"]) +
    unknown_formula["step_score"] + 
    unknown_formula["rdbe_score"] 
  
  
  # hist(unknown_formula$cplex_score)
  # length(unknown_formula$cplex_score[unknown_formula$cplex_score<1])
  # print("Finish scoring formula.")
  return(unknown_formula)
}
### Score_edge_cplex ####
Score_edge_cplex = function(CPLEXset, edge_bonus = -log10(0.5))
{
  edge_info_sum = CPLEXset$data$edge_info_sum %>%
    arrange(edge_ilp_id)
  unknown_formula = CPLEXset$data$unknown_formula
  
  edge_info_sum$edge_score = log10(edge_info_sum$edge_score) + edge_bonus
  
  test3 = edge_info_sum[edge_info_sum$ILP_id2<=nrow(unknown_formula)&
                          edge_info_sum$ILP_id1<=nrow(unknown_formula),]
  
  test3_same12 = test3[test3$formula1==test3$formula2,]
  df_same12 = table(test3_same12$formula1)
  
  test3_dif12 = test3[test3$formula1!=test3$formula2,]
  test3_dif12 = test3_dif12[duplicated(test3_dif12[,c("formula1","formula2")]) | 
                              duplicated(test3_dif12[,c("formula1","formula2")], fromLast=TRUE),]
  temp_merge = with(test3_dif12, paste0(formula1, formula2))
  df_dif12 = table(temp_merge)
  
  # Scenerio: Each duplicated mz will generate duplicated edge connectoin, exaggerating the score
  # Solution: when n duplicated mz exist, the intra mz edge number are n*(n+1)/2, but effectively they should only get n/2
  
  sol_mat = data.frame(n=1:max(df_same12, df_dif12), div = NA)
  sol_mat$div = (-1+sqrt(1+8*sol_mat$n))/2
  
  test3_same12$edge_score = test3_same12$edge_score / (sol_mat$div[df_same12[test3_same12$formula1]])
  test3_dif12$edge_score = test3_dif12$edge_score / (sol_mat$div[df_dif12[temp_merge]])
  
  temp_edge_info_sum = rbind(test3_same12,test3_dif12,edge_info_sum)
  temp_edge_info_sum = temp_edge_info_sum[!duplicated(temp_edge_info_sum$edge_ilp_id),]
  
  temp_edge_info_sum = temp_edge_info_sum[with(temp_edge_info_sum, order(edge_ilp_id)),]
  
  # print("Finish scoring edges.")
  return(temp_edge_info_sum)
  
}

## Run_CPLEX ####
Run_CPLEX = function(CPLEXset, obj){
  # obj = obj_cplex
  env <- openEnvCPLEX()
  prob <- initProbCPLEX(env)
  
  nc = CPLEXset$para$nc
  nr = CPLEXset$para$nr
  CPX_MAX = CPLEXset$para$CPX_MAX
  rhs = CPLEXset$para$rhs
  sense = CPLEXset$para$sense
  beg = CPLEXset$para$beg
  cnt = CPLEXset$para$cnt
  ind = CPLEXset$para$ind
  val = CPLEXset$para$val
  lb = CPLEXset$para$lb
  ub = CPLEXset$para$ub
  ctype = CPLEXset$para$ctype
  
  
  
  copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj = obj, rhs, sense,
                    beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
  
  
  copyColTypeCPLEX(env, prob, ctype)
  
  # Conserve memory true
  setIntParmCPLEX(env, CPX_PARAM_MEMORYEMPHASIS, CPX_ON)
  # setDefaultParmCPLEX(env)
  # getChgParmCPLEX(env)
  
  # tictoc::tic()
  
  return_code = mipoptCPLEX(env, prob)
  result_solution=solutionCPLEX(env, prob)
  
  # print(paste(return_codeCPLEX(return_code),"-",
  #             status_codeCPLEX(env, getStatCPLEX(env, prob)),
  #             " - OBJ_value =", result_solution$objval))
  
  # tictoc::toc()
  
  # writeProbCPLEX(env, prob, "prob.lp")
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
  
  return(list(obj = obj, result_solution = result_solution))
}

## Read CPLEX result ####
Read_CPLEX_result = function(solution){
  CPLEX_all_x=list()
  for(i in 1:length(solution)){
    CPLEX_all_x[[i]] =  solution[[i]]$result_solution$x
  }
  CPLEX_all_x = bind_cols(CPLEX_all_x)
  return(CPLEX_all_x)
}

## CPLEX_permutation ####
CPLEX_permutation = function(CPLEXset, n_pmt = 5, sd_rel_max = 0.5){
  unknown_formula = CPLEXset$data$unknown_formula
  obj = CPLEXset$Init_solution[[1]]$obj
  obj_node = obj[1:nrow(unknown_formula)]
  obj_edge = obj[(nrow(unknown_formula)+1):length(obj)]
  solution_ls = list()
  for(i in 1:n_pmt){
    
    temp_obj_edge = obj_edge + rnorm(length(obj_edge), mean = 0, sd = max(obj_edge) * sd_rel_max)
    temp_obj <- c(obj_node, temp_obj_edge)
    
    result_solution = Run_CPLEX(CPLEXset, temp_obj)
    solution_ls[[length(solution_ls)+1]] = result_solution
    
  }
  return(solution_ls)
}

## expand_formula_to_library ####
expand_formula_to_library = function(formula = "C2H4O2"){
  data(isotopes)
  formula = check_chemform(isotopes,formula)$new_formula
  data.frame(ID = 1:length(formula),
             Name = formula,
             MF = formula,
             Exact_Mass = formula_mz(formula),
             category = 0, 
             rdbe = formula_rdbe(formula),
             stringsAsFactors = F
  )
}


# !!Run_function!! ####
mz_calculator = function(raw_data, 
                         library_data,
                         connect_depth = 5,
                         ion_mode = 1,
                         rule_table_file = "~/myrepo/mz_calculator/dependent/connect_rules.csv"
                         
)
{
  # raw_data = test_rawdata
  # library_data = expand_formula_to_library("C5H7N3O")
  {
    Mset = list()
    Mset[["Raw_data"]] = raw_data
    Mset[["Library"]] = library_data
    Mset[["Connect_rules"]]=Read_rule_table(rule_table_file = rule_table_file)
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
                                                   biotransform_step = connect_depth,
                                                   propagation_score_threshold = 0.25,
                                                   top_n = 50)
    
    CPLEXset = Prepare_CPLEX(Mset, EdgeSet)
  }
  
  # Run CPLEX ####
  {
    CPLEXset$data$unknown_formula = Score_formula(Mset, CPLEXset,
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
  
  {
    unknown_nodes = CPLEXset$data$unknown_nodes[,1:3]
    unknown_formula = CPLEXset$data$unknown_formula
    
    unknown_formula["ILP_result"] = CPLEX_x[1:nrow(unknown_formula)]
    unknown_formula_CPLEX = unknown_formula[unknown_formula$ILP_result !=0,]
    print(paste("pred formula num =", nrow(unknown_formula_CPLEX)))
    
    unknown_node_CPLEX = merge(unknown_nodes,unknown_formula_CPLEX,by.x = "ID", by.y = "id",all=T)
  }
  
  return(unknown_node_CPLEX)
  
} 

