#显示中文 
#Sys.setlocale(category = "LC_ALL", locale = "Chinese")


#import library
{
  library(readr)
  #install.packages("rcdk")
  #library(rcdk)
  library(igraph)
  #install.packages("fitdistrplus")
  library("fitdistrplus")
  library(tidyr)
  library(dplyr)
  library(enviPat)
  # setwd("C:/Users/Li Chen/Desktop/Github local")
  # devtools::install("lc8")
  library(lc8)
  
  #install.packages("stringi")
  library(stringi)
  #install.packages("matrixStats")
  library(matrixStats)
}

# Fucntions ####
# Function for parsing#### 
## Data name and cohorts####
Cohort_Info = function(mset)
{
  raw = mset$Raw_data
  all_names=colnames(raw)[15:ncol(raw)]
  
  if(length(grep("blank|blk", all_names, ignore.case = T))!=0){
    sample_names=all_names[-grep("blank|blk", all_names, ignore.case = T)]
  } else {
    sample_names=all_names
  }
  blank_names=all_names[grep("blank|blk", all_names, ignore.case = T)]
  sample_cohort=stri_replace_last_regex(sample_names,'-1|-2|-3|-4|-a|-b|-c|_mean', '',stri_opts_regex(case_insensitive=T))
  
  return(list("sample_names"=sample_names,"blank_names"=blank_names, "sample_cohort"=sample_cohort))
}

## Clean up duplicate peaks from peak picking####
Peak_cleanup = function(mset,
                        ms_dif_ppm=5/10^6, 
                        rt_dif_min=0.2,
                        detection_limit=2500
                        )
{
  raw = mset$Raw_data
  raw = raw[complete.cases(raw[, 14:ncol(raw)]),]
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  raw$medMz = raw$medMz*abs(mset$Global_parameter$mode) - (H_mass-e_mass)*mset$Global_parameter$mode
  
  #Group MS groups
  {
    s = raw[with(raw, order(medMz, medRt)),]
    s["merge_group"]=NA
    mgMS_count = 1 
    i_max=i_min=1
    
    #
    while (i_min <= nrow(raw)){
      while(s$medMz[i_max]-s$medMz[i_min]< (s$medMz[i_min]*ms_dif_ppm)){
        i_max = i_max+1
        if(i_max>nrow(raw)){
          break
        }
      }
      if(i_max-i_min == 1){
        s$merge_group[i_min] = 0
        i_min = i_min + 1
        next
      }
      while(i_min < i_max){
        s$merge_group[i_min] = mgMS_count
        i_min=i_min+1
      }  
      mgMS_count=mgMS_count+1
    }
  }
  
  #Group RT similar groups based on MS groups
  {
    s2 = s[with(s, order(merge_group, medRt)),]
    s2["RTmerge_group"]=NA
    
    mgRT_count=1
    
    j_min=length(which(s2$merge_group == 0))
    
    for (n in 1:j_min){
      s2$RTmerge_group[n] = 0
    }
    
    j_min=j_min+1
    if(j_min <=nrow(raw)){s2$RTmerge_group[j_min] = 0}
    j_max=j_min+1
    
    while (j_max <= nrow(raw)){
      while (s2$merge_group[j_min] != s2$merge_group[j_max] |
             s2$medRt[j_max] - s2$medRt[j_min] >= rt_dif_min){
        s2$RTmerge_group[j_max] = 0
        j_max = j_max+1
        j_min = j_min+1
        if(j_max > nrow(raw)){break}
      }
      if(j_max > nrow(raw)){break}
      
      while (s2$medRt[j_max] - s2$medRt[j_max-1] < rt_dif_min){
        j_max = j_max+1
        if(j_max>nrow(raw) | s2$merge_group[j_min] != s2$merge_group[j_max]) {break}
      }
      j_max = j_max-1
      
      for(j in j_min:j_max){
        s2$RTmerge_group[j] = mgRT_count
      }
      
      mgRT_count = mgRT_count+1
      j_min=j_max
      j_max=j_max+1
    }
  }
  
  #Flag groups for deletion & combine signal
  {
    s3 = s2[with(s2, order(RTmerge_group)),]
    s3[["flag"]]=NA
    k_max=k_min=length(which(s3$RTmerge_group==0))
    for(k in 1:k_min){
      s3$flag[k] = T
    }
    
    while (k_max <= nrow(raw)){
      k_min = k_max
      s3$flag[k_min] = T
      while (s3$RTmerge_group[k_min] == s3$RTmerge_group[k_max]){
        k_max = k_max+1
        if(k_max > nrow(raw)){break}
        s3$flag[k_max] = F
      }
      
      if(k_min == k_max-1){
        next
      }
      
      s3$medMz[k_min]=mean(s3$medMz[k_min:(k_max-1)])
      s3$medRt[k_min]=mean(s3$medRt[k_min:(k_max-1)])
      s3$goodPeakCount[k_min]=max(s3$goodPeakCount[k_min:(k_max-1)])
      s3$groupId[k_min]=min(s3$groupId[k_min:(k_max-1)])
      for (n in 14:ncol(raw)){
        s3[k_min,n]=max(s3[k_min:(k_max-1),n])
      }
    }
  }
  
  #intermediate files, replace below detection number to random small number
  {
    
    s4=s3
    # s4[,4:ncol(s4)][s4[,4:ncol(s4)]<detection_limit]=sample(1:detection_limit, 
    #                                                         size=sum(s4[,4:ncol(s4)]<detection_limit), 
    #                                                         replace=T)
    
    s4["mean_inten"]=NA
    s4["mean_inten"]=rowMeans(s4[,mset$Cohort$sample_names])
    s4$flag[s4$mean_inten<detection_limit]=F
  }
  
  
  s5 = s4[s4$flag, 1:ncol(raw)]
  s5 = s5[with(s5, order(groupId)),]
  s5$groupId = 1:nrow(s5)
  return(s5)
}

## Identify peaks with high blanks####
High_blank = function(mset, fold_cutoff = 2)
{
  s7 = mset$Data
  
  s7["high_blank"]=NA
  if(length(mset$Cohort$blank_names)>0){
    s7["high_blank"]= rowMeans(s7[,mset$Cohort$sample_names]) < fold_cutoff*rowMeans(s7[,mset$Cohort$blank_names])
  }
  
  result = s7[,c("groupId","high_blank")]
  return(result)  
}

#Annontate base on library mz
library_match = function(mset, ppm=5/10^6, library_file = "hmdb_unique.csv")
{
  s5 = mset$Data
  s5 = s5[with(s5,order(medMz)),]
  
  
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  mode = mset$Global_parameter$mode
  
  hmdb_df = read_csv("hmdb_unique.csv")
  hmdb_df["adjustmz"] = hmdb_df$Exact_Mass + (H_mass-e_mass)*mode
  hmdb_df=hmdb_df[with(hmdb_df, order(adjustmz)),]
  
  i_min=1
  j=1
  
  library_match = list()
  
  while(j<=nrow(s5)){
    while(hmdb_df$adjustmz[i_min+1] < (s5$medMz[j]*(1-ppm))){
      i_min=i_min+1
    }
    temp_metabolite=as.character()
    temp_formula=as.character()
    k=1
    while(hmdb_df$adjustmz[i_min+k]<(s5$medMz[j]*(1+ppm))){
      k=k+1
    }
    if(k>1){
      library_match[[s5$groupId[j]]]=hmdb_df[(i_min+1):(i_min+k-1),]
    }
    else{library_match[[s5$groupId[j]]]=hmdb_df[0,]}
    j=j+1
  }
  
  
  
  return (library_match)
}

## Statistical analysis with MetaboAnalyst####
Metaboanalyst_Statistic = function(mset)
{
  library(MetaboAnalystR)
  
  MA_output = mset$Data[,c("groupId",  mset$Cohort$sample_names)]
  MA_output = rbind(c("cohort",mset$Cohort$sample_cohort),MA_output)

  write.csv(MA_output, file="MetaboAnalyst_file.csv", row.names=F)
  
  mSet <- InitDataObjects("pktable", "stat", FALSE)
  mSet<-Read.TextData(mSet, "MetaboAnalyst_file.csv", "colu", "disc");
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet);
  mSet<-FilterVariable(mSet, "iqr", "F", 25)
  mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", "a11", ratio=FALSE, ratioNum=20)
  mSet<-ANOVA.Anal(mSet, F, 0.05, "fisher")
  
  ANOVA_file = "anova_posthoc.csv"
  ANOVA_raw <- read_csv(ANOVA_file)
  
  ANOVA_FDR = ANOVA_raw[,c("X1","FDR")]
  ANOVA_FDR$FDR=-log10(ANOVA_FDR$FDR)
  colnames(ANOVA_FDR)=c("ID","_log10_FDR")
  
  return(ANOVA_FDR)
  
}



# Function for network ######

## Merge experiment and library nodes####
Form_node_list = function(mset, library_file = "hmdb_unique.csv")
{
  NodeSet = list()
  NodeSet[["Expe"]] = mset$Data[,c("groupId","medMz","medRt","formula")]
  NodeSet$Expe["origin"]="Experiment"
  NodeSet$Expe["compound_name"]=NA
  colnames(NodeSet$Expe) = c("ID","mz","RT","MF", "origin","compound_name")
  
  NodeSet[["Library"]] = read.csv(library_file,stringsAsFactors = F)
  NodeSet$Library = cbind(NodeSet$Library, RT=NA)
  colnames(NodeSet$Library) = c("ID","compound_name","MF","mz","origin", "RT")
  NodeSet$Library$ID = 1:nrow(NodeSet$Library)+nrow(NodeSet$Expe)
  data("isotopes")
  NodeSet$Library$MF = check_chemform(isotopes, NodeSet$Library$MF)$new_formula
  NodeSet$Library$origin="Library"
  
  merge_node_list = rbind(NodeSet$Expe,NodeSet$Library )

  return(merge_node_list)
}

## Define transformation and artifact rules####
Read_rule_table = function(rule_table_file = "biotransform.csv"){
  library(enviPat)
  data("isotopes")
  biotransform = read.csv(rule_table_file,stringsAsFactors = F)
  biotransform$Formula = check_chemform(isotopes,biotransform$Formula)$new_formula
  biotransform = biotransform[with(biotransform, order(mass)),]
  return(biotransform)
}
## Edge_list for biotransformation ####
Edge_biotransform = function(mset, mass_abs = 0.001, mass_ppm = 5/10^6)
{
  
  
  
  merge_node_list = mset$NodeSet[with(mset$NodeSet, order(mz)),]
  merge_nrow = nrow(merge_node_list)
  
  timer=Sys.time()
  {
    edge_ls = list()
    temp_mz_list=merge_node_list$mz
    
    for (k in 1:nrow(mset$Biotransform)){
      temp_fg=mset$Biotransform$mass[k]
      i=j=1
      temp_edge_list = data.frame(node1=as.numeric(), node2=as.numeric(), linktype=as.numeric(), mass_dif=as.numeric())
      while(i<=merge_nrow){
        temp_ms=0
        mass_tol = max(temp_mz_list[i]*mass_ppm,0.001)
        while(1){
          j=j+1
          if(j>merge_nrow){break}
          temp_ms = temp_mz_list[j]-temp_mz_list[i]
          
          if(abs(temp_ms-temp_fg)<mass_tol){
            temp_edge_list[nrow(temp_edge_list)+1,]=list(merge_node_list$ID[i], merge_node_list$ID[j], k, (temp_ms-temp_fg)/temp_mz_list[j]*1E6)
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
      print(paste("Biotransform", mset$Biotransform$category[k], nrow(temp_edge_list),"found."))
    }
  }
  
  print(Sys.time()-timer)
  edge_list = bind_rows(edge_ls)
  
  edge_list$linktype=mset$Biotransform$Formula[edge_list$linktype]
  
  edge_list_sub = subset(edge_list, 
                         (edge_list$node1<=nrow(mset$Data)|
                            edge_list$node2<=nrow(mset$Data))&
                           edge_list$node1!=edge_list$node2
  )
  
  
  edge_list_sub["category"]=1
  return(edge_list_sub)
}
### Scoring edge based on mass accuracy####
Edge_score = function(Biotransform)
{
  edge_mzdif_FIT <- fitdist(as.numeric(Biotransform$mass_dif), "norm")    
  hist(Biotransform$mass_dif)
  
  plot(edge_mzdif_FIT)  
  Biotransform["edge_massdif_score"]=dnorm(Biotransform$mass_dif, edge_mzdif_FIT$estimate[1], edge_mzdif_FIT$estimate[2])
  Biotransform["edge_massdif_score"]=Biotransform["edge_massdif_score"]/max(Biotransform["edge_massdif_score"])
  return(Biotransform)
}
## Network_prediction used to connect nodes to library and predict formula####
Network_prediction = function(mset, edge_list_sub, 
                              top_formula_n = 1
                              )
{
  
  mnl=mset$NodeSet
  edge_list_sub=test3
  
  #Initialize predict_formula from HMDB known formula
  {
    data_str = data.frame(id=as.numeric(),
                          formula=as.character(), 
                          steps=as.numeric(), 
                          parent=as.numeric(), 
                          score=as.numeric(), stringsAsFactors = F)
    #sf stands for summary_formula
    
    
    
    sf = lapply(1:nrow(mnl),function(i)(data_str))
    
    for(i in (1+nrow(mset$Data)):nrow(mnl)){
      Initial_formula =sf[[i]]
      Initial_formula[1,]= list(i,mnl$MF[i],0,1, 1)
      sf[[i]]=Initial_formula
      
    }
  }
  
  #while loop to Predict formula based on known formula and edgelist 
  
  
  step=0
  timer=Sys.time()
  New_nodes_in_network = 1
  while(New_nodes_in_network==1){
    
    all_nodes_df = bind_rows(lapply(sf, head,top_formula_n))
    new_nodes_df = all_nodes_df[all_nodes_df$steps==step 
                                &all_nodes_df$score>0,]
    print(paste("nrow",nrow(all_nodes_df),"in step",step,"elapsed="))
    print((Sys.time()-timer))
    
    if(nrow(new_nodes_df)==0){break}
    
    
    
    #Handling head
    {
      edge_list_node1 = edge_list_sub[edge_list_sub$node1 %in% new_nodes_df$id,]
      head_list = edge_list_node1$node1
      n=1
      for (n in 1: nrow(new_nodes_df)){
        
        head = new_nodes_df$id[n]
        head_formula = new_nodes_df$formula[n]
        
        temp_edge_list=subset(edge_list_node1, edge_list_node1$node1==head)
        if(nrow(temp_edge_list)==0){next}
        
        head_score = sf[[head]]$score[match(head_formula,sf[[head]]$formula)]
        
        if(head_score==0){next}
        i=1
        for(i in 1:nrow(temp_edge_list)){
          tail=temp_edge_list$node2[i]
          temp_fg = temp_edge_list$linktype[i]
          if(temp_fg==""){
            temp_formula = head_formula
          }else{
            temp_formula=my_calculate_formula(head_formula,temp_fg,1,Is_valid = T)
          }
          
          if(is.logical(temp_formula))#function return false if not valid formula
          { temp_score=0
          temp_formula=paste(head_formula,temp_fg,sep="+")
          }else{
            temp_score = head_score*temp_edge_list$edge_massdif_score[i]
          }
          
          temp_parent = head
          temp_steps = step+1
          
          #Criteria to enter new entry into formula list
          #1. new formula
          temp = sf[[tail]]
          temp_subset=subset(temp, temp$formula==temp_formula)
          if(nrow(temp_subset)!=0){
            #2. much higher scores
            if(temp_score<=(1.2*max(temp_subset$score))){
              next
            }
          }
          #Enter new entry
          temp[nrow(temp)+1, ] = list(tail, temp_formula, temp_steps, temp_parent, temp_score)
          temp = temp[with(temp, order(-score)),]
          sf[[tail]]=temp
        }
      }
    }
    
    #Handling tail
    {
      edge_list_node2 = edge_list_sub[edge_list_sub$node2 %in% new_nodes_df$id,]
      tail_list = edge_list_node2$node2
      
      for (n in 1: nrow(new_nodes_df)){
        
        tail = new_nodes_df$id[n]
        tail_formula = new_nodes_df$formula[n]
        
        temp_edge_list=subset(edge_list_node2, edge_list_node2$node2==tail)
        if(nrow(temp_edge_list)==0){next}
        
        tail_score = sf[[tail]]$score[match(tail_formula,sf[[tail]]$formula)]
        
        if(tail_score==0){next}
        i=1
        for(i in 1:nrow(temp_edge_list)){
          head=temp_edge_list$node1[i]
          temp_fg = temp_edge_list$linktype[i]
          if(temp_fg==""){
            temp_formula = tail_formula
          }else{
            temp_formula=my_calculate_formula(tail_formula,temp_fg,-1,Is_valid = T)
          }
          
          if(is.logical(temp_formula))#function return false if not valid formula
          { temp_score=0
          temp_formula=paste(tail_formula,temp_fg,sep="-")
          }else{
            temp_score = tail_score*temp_edge_list$edge_massdif_score[i]
          }
          
          temp_parent = tail
          temp_steps = step+1
          
          #Criteria to enter new entry into formula list
          #1. new formula
          temp = sf[[head]]
          temp_subset=subset(temp, temp$formula==temp_formula)
          if(nrow(temp_subset)!=0){
            #2. much higher scores
            if(temp_score<=(1.2*max(temp_subset$score))){
              next
            }
          }
          #Enter new entry
          temp[nrow(temp)+1, ] = list(head, temp_formula, temp_steps, temp_parent, temp_score)
          temp = temp[with(temp, order(-score)),]
          sf[[head]]=temp
        }
      }
    }
    
    step=step+1
    
    
  }
  return(sf)
}


## Variance between peaks####
Peak_variance = function(mset, 
                         time_cutoff=0.1,
                         mass_cutoff=0,
                         correlation_cutoff = 0.7)
{
  df_raw = mset$Data[,c("groupId","medMz","medRt",mset$Cohort$sample_names)]
  df_raw["mean_inten"]=rowMeans(df_raw[,mset$Cohort$sample_names])
  df_raw["log10_inten"]=log10(df_raw$mean_inten)
  
  {
    df_raw = df_raw[with(df_raw, order(medRt)),]
    i_min = i_max =1
    
    edge_list = list()
    timer = Sys.time()
    while(i_min!= nrow(df_raw)){
      if(i_min %% 1000 ==0 ){
        print(paste("nrow",i_min,"elapsed="))
        print((Sys.time()-timer))
      }
      
      if(df_raw$mean_inten[i_min]<mass_cutoff){
        i_min = i_min+1
        next
      }

      
      
      
      temp_t = df_raw$medRt[i_min]
      temp_t_end = temp_t+time_cutoff
      
      while(df_raw$medRt[i_max]<temp_t_end){
        if(i_max == nrow(df_raw)){break}
        i_max = i_max+1
      }
      i_max = i_max-1
      
      temp_df_raw = df_raw[i_min:i_max,]
      
      temp_df_raw$time_dif=temp_df_raw$medRt-temp_df_raw$medRt[1]
      temp_df_raw$mz_dif = round(temp_df_raw$medMz-temp_df_raw$medMz[1], digits=5)
      
      temp_x = t(temp_df_raw[1,mset$Cohort$sample_names])
      temp_y = t(temp_df_raw[,mset$Cohort$sample_names])
      temp_df_raw$correlation = as.numeric(t(cor((temp_x), (temp_y))))

      temp_df_raw["node1"]=temp_df_raw$groupId[1]
      temp_df_raw["node2"]=temp_df_raw$groupId
      temp_df_raw["mz_node1"] = temp_df_raw$medMz[1]
      temp_df_raw["mz_node2"] = temp_df_raw$medMz
      temp_df_raw["log10_inten_node1"] = temp_df_raw$log10_inten[1]
      temp_df_raw["log10_inten_node2"] = temp_df_raw$log10_inten
      temp_df_raw["log10_inten_ratio"] = temp_df_raw["log10_inten_node2"]-temp_df_raw["log10_inten_node1"]
      
      temp_edge_ls=temp_df_raw[,c("node1","node2","time_dif", "mz_dif", "correlation","mz_node1",
                                  "mz_node2","log10_inten_node1","log10_inten_node2",
                                  "log10_inten_ratio")]
      
      # cor(temp_x, temp_y, method = "kendall")
      # cor(temp_x, temp_y, method = "spearman")
      # cosine(tem#p_x, temp_y)
      temp_edge_ls=temp_edge_ls[temp_edge_ls$correlation>0.8,]
      edge_list[[i_min]]=temp_edge_ls
      i_min = i_min+1
    }
  }
  
  edge_ls = bind_rows(edge_list)
  edge_ls=edge_ls[edge_ls$node1!=edge_ls$node2,]
  rm(edge_list)
  
  #Switch directionality
  {
    edge_ls = edge_ls[with(edge_ls, order(mz_dif)),]
    temp_data = edge_ls[edge_ls$mz_dif<0,]
    
    temp_data$mz_node2=temp_data$mz_node1
    temp_data$mz_node1=temp_data$mz_node1+temp_data$mz_dif
    
    temp_node=temp_data$node1
    temp_data$node1=temp_data$node2
    temp_data$node2=temp_node
    temp_data$time_dif=-temp_data$time_dif
    temp_data$mz_dif=-temp_data$mz_dif
    temp_data$log10_inten_ratio=-temp_data$log10_inten_ratio
    
    edge_ls[edge_ls$mz_dif<0,] = temp_data 
    rm(temp_data)
  }
  
  edge_ls = edge_ls[with(edge_ls, order(mz_dif)),]
  return(edge_ls)
}








### Edge_list for artifacts####
Artifact_prediction = function(mset, Peak_inten_correlation, search_ms_cutoff=0.001){
  edge_ls_highcor=Peak_inten_correlation
  edge_ls_highcor = edge_ls_highcor[with(edge_ls_highcor, order(mz_dif)),]
  junk_df = mset$Artifacts
  
  
  i=j=1
  temp_ls = list()
  temp_df = edge_ls_highcor[1,]
  temp_df["category"]=NA
  temp_df["linktype"]=NA
  temp_df["mass_dif"]=NA
  
  temp_df = temp_df[0,]
  
  while (i <= nrow(edge_ls_highcor)){
    i=i+1
    if(edge_ls_highcor$mz_dif[i]<(junk_df$mass[j]-search_ms_cutoff)){
      next
    }
    if(edge_ls_highcor$mz_dif[i]>(junk_df$mass[j]+search_ms_cutoff)){
      temp_ls[[j]]=temp_df
      temp_df = temp_df[0,]
      j=j+1
      if(j>nrow(junk_df)){break}
      #search_ms_cutoff = 10 * initial_FIT$estimate[2]/1000 
      while(edge_ls_highcor$mz_dif[i]>(junk_df$mass[j]-search_ms_cutoff)){i=i-1}
      next
    }
    temp_df[(nrow(temp_df)+1),]=c(edge_ls_highcor[i,],
                                  as.character(junk_df$Symbol[j]), 
                                  junk_df$Formula[j], 
                                  (edge_ls_highcor$mz_dif[i]-junk_df$mass[j])/edge_ls_highcor$mz_node1[i]*10^6
    )
  }
  
  
  temp_df_isotope = bind_rows(temp_ls)
  
  junk_summary=table(temp_df_isotope$category)
  print(junk_summary)
  
  
  
  
  
  #Oligomers. Note that it is indistinguishable between 2-charge-parent pair and parent-dimer pair
  
  test_time = Sys.time()
  {
    H_mass = 1.00782503224
    e_mass = 0.00054857990943
    mode=-1
    ppm=5/10^6
    temp_df_oligo = temp_df[0,]
    
    df_raw = mset$Data
    
    temp_edge_ls = data.frame(ratio=edge_ls_highcor$mz_dif/edge_ls_highcor$mz_node1)
    temp_edge_ls["rounding"]=round(temp_edge_ls[,1],digit=0)
    temp_edge_ls["dif"]=temp_edge_ls["rounding"]-temp_edge_ls[,1]
    temp_edge_ls = temp_edge_ls[(temp_edge_ls$rounding!=0)
                                &(abs(temp_edge_ls$dif)<0.05),]
    
    # for(i in 1:nrow(temp_edge_ls)){
    #   temp_data = edge_ls_highcor[rownames(temp_edge_ls)[i],]
    #   temp_mz1 = df_raw$medMz[which(df_raw$groupId==temp_data$node1)]
    #   temp_mz2 = temp_mz1 + temp_data$mz_dif
    #   for(j in 2:10){
    #     #if charge data
    #     #if(abs(temp_mz1*j-(H_mass-e_mass)*mode*(j-1)-temp_mz2)<temp_mz2*ppm){
    #     #if neutral data
    #     if(abs(temp_mz1*j-temp_mz2)<search_ms_cutoff){
    #       temp_df_oligo[(nrow(temp_df_oligo)+1),]=c(temp_data,"oligomer",paste("x",j,sep=""), (temp_mz1*j-temp_mz2)/temp_mz1*10^6)
    #     }
    #   }
    # }
  }
  rm(temp_edge_ls)
  temp_df_oligo
  nrow(temp_df_oligo)
  
  oligo_summary=table(temp_df_oligo$linktype)
  print(oligo_summary)
  
  #data$parent[data$groupId==6]
  test_time = Sys.time()-test_time
  
  
  edge_ls_annotate=rbind(temp_df_isotope,temp_df_oligo)
  
  
  edge_ls_annotate_network = edge_ls_annotate[,c("node1","node2","linktype","mass_dif","category")]
  
  
  
  return(edge_ls_annotate_network)
}
# Main Codes ####
## Read files ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  filename = c("Yeast-Ecoli-neg-peakpicking_blank.csv")
  mset = list()
  mset[["Raw_data"]] <- read_csv(filename)
  mset[["Raw_data"]] = mset$Raw_data[base::sample(nrow(mset$Raw_data),10000),]
}

## Initialise ####
{
  mset[["Global_parameter"]]=  list(mode = -1,
                                    normalized_to_col_median = F)
  mset[["Cohort"]]=Cohort_Info(mset)
  
  #Clean-up duplicate peaks 
  mset[["Data"]] = Peak_cleanup(mset,
                                ms_dif_ppm=5/10^6, 
                                rt_dif_min=0.2,
                                detection_limit=2500)
  #View(mset$Data)
  mset[["ID"]]=mset$Data$groupId
}

## Feature generation ####
{
  #Identify peaks with high blanks
  mset[["High_blanks"]]=High_blank(mset, fold_cutoff = 2)
  #View(mset$High_blanks)
  
  #library_match
  mset[["library_match"]] = library_match(mset, ppm=5/10^6, library_file = "hmdb_unique.csv")
  counts=unlist(lapply(mset$library_match,nrow))
  num_of_library_match = cbind(ID=mset$ID,counts)
  mset[["library_match"]][["num_of_library_match"]]=num_of_library_match
  
  #Metaboanalyst_Statistic
  mset[["Metaboanalyst_Statistic"]]=Metaboanalyst_Statistic(mset)
} 





# Network ####

{
  
  EdgeSet = list()
  
  mset[["NodeSet"]]=Form_node_list(mset, library_file = "hmdb_unique.csv")
  mset[["Biotransform"]]=Read_rule_table(rule_table_file = "biotransform.csv")
  
  EdgeSet[["Biotransform"]] = Edge_biotransform(mset)
  EdgeSet[["Biotransform"]] = Edge_score(EdgeSet$Biotransform)
  
  
  
  mset[["Artifacts"]]=Read_rule_table(rule_table_file = "artifacts.csv")
  EdgeSet[["Peak_inten_correlation"]] = Peak_variance(mset,
                                                      time_cutoff=0.1,
                                                      mass_cutoff = 2e4,
                                                      correlation_cutoff = 0.8)
  EdgeSet[["Artifacts"]] = Artifact_prediction(mset, 
                                               EdgeSet$Peak_inten_correlation, 
                                               search_ms_cutoff=0.001)
  EdgeSet[["Artifacts"]] = Edge_score(EdgeSet$Artifacts)
  
  

  mset[["NodeSet_network"]] = Network_prediction(mset, 
                                                 rbind(EdgeSet$Artifacts,EdgeSet$Biotransform), 
                                                 top_formula_n = 2)
  
  
  sf=mset[["NodeSet_network"]]
  All_formula_predict=bind_rows(mset[["NodeSet_network"]])
  test = test = All_formula_predict[grepl("\\[", All_formula_predict$formula) & All_formula_predict$score>0.4,]
  #write.csv(All_formula_predict, "All_formula_predict.csv",row.names = F)
  
}






