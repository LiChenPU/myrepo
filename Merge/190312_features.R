# 显示中文
# Sys.setlocale(category = "LC_ALL", locale = "Chinese")
# !diagnostics off

# Import library ####
{
  library(readr)
  library(igraph)
  library(fitdistrplus)
  library(tidyr)
  library(dplyr)
  library(enviPat)
  library(stringi)
  library(matrixStats)
  library(Matrix)
  library(slam)
  library(cplexAPI)
  
  #devtools::install_github("LiChenPU/Formula_manipulation")
  library(lc8)
  
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

# Fucntions ####
# Function for parsing#### 
## read_library ####
read_library = function(library_file = "hmdb_unique.csv"){
  data(isotopes)
  hmdb_lib = read_csv(library_file)
  hmdb_lib$MF = check_chemform(isotopes, hmdb_lib$MF)$new_formula
  hmdb_lib$Exact_mass = formula_mz(hmdb_lib$MF)

  return(hmdb_lib)
}
  
## Cohort_Info - Data name and cohorts ####
Cohort_Info = function(Mset)
{
  raw = Mset$Raw_data
  all_names=colnames(raw)[15:ncol(raw)]
  
  if(length(grep("blank|blk", all_names, ignore.case = T))!=0){
    sample_names=all_names[-grep("blank|blk", all_names, ignore.case = T)]
  } else {
    sample_names=all_names
  }
  blank_names=all_names[grep("blank|blk", all_names, ignore.case = T)]
  sample_cohort=stri_replace_last_regex(sample_names,'_\\d+|-\\d+|-a|-b|-c|_mean', '',stri_opts_regex(case_insensitive=T))
  
  return(list("sample_names"=sample_names,"blank_names"=blank_names, "sample_cohort"=sample_cohort))
}

## Peak_cleanup - Clean up duplicate peaks from peak picking ####
Peak_cleanup = function(Mset,
                        ms_dif_ppm=1/10^6, 
                        rt_dif_min=0.02,
                        detection_limit=2500
                        )
{
  
  raw = Mset$Raw_data
  raw = raw[complete.cases(raw[, (1+which(colnames(raw)=="parent")):ncol(raw)]),]
  colnames(raw)[colnames(raw)=="groupId"] = "ID"
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  raw$medMz = raw$medMz*abs(Mset$Global_parameter$mode) - (H_mass-e_mass)*Mset$Global_parameter$mode
  
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
      s3$ID[k_min]=min(s3$ID[k_min:(k_max-1)])
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
    s4["mean_inten"]=rowMeans(s4[,Mset$Cohort$sample_names])
    s4$flag[s4$mean_inten<detection_limit]=F
  }
  
  
  s5 = s4[s4$flag, 1:ncol(raw)]
  s5 = s5[with(s5, order(ID)),]
  s5$ID = 1:nrow(s5)
  
  s5["mean_inten"]=rowMeans(s5[,Mset$Cohort$sample_names])
  s5["log10_inten"]=log10(s5$mean_inten)
  
  return(s5)
}

## High_blank - Identify peaks with high blanks ####
High_blank = function(Mset, fold_cutoff = 2)
{
  s7 = Mset$Data
  
  s7["high_blank"]=NA
  if(length(Mset$Cohort$blank_names)>0){
    s7["high_blank"]= rowMeans(s7[,Mset$Cohort$sample_names]) < fold_cutoff*rowMeans(s7[,Mset$Cohort$blank_names])
  }
  
  result = s7[,c("ID","high_blank")]
  return(result)  
}

## library_match - Annontate base on library mz ####
library_match = function(Mset, ppm=5/10^6, library_file = "hmdb_unique.csv")
{
  s5 = Mset$Data
  s5 = s5[with(s5,order(medMz)),]
  
  
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  mode = Mset$Global_parameter$mode
  
  hmdb_df = Mset$Library
  data(isotopes)
  hmdb_df$MF = check_chemform(isotopes, hmdb_df$MF)$new_formula
  hmdb_df=hmdb_df[with(hmdb_df, order(Exact_Mass)),]
  
  i_min=1
  j=1
  
  library_match = list()
  
  while(j<=nrow(s5)){
    while(hmdb_df$Exact_Mass[i_min+1] < (s5$medMz[j]*(1-ppm))){
      i_min=i_min+1
    }
    temp_metabolite=as.character()
    temp_formula=as.character()
    k=1
    while(hmdb_df$Exact_Mass[i_min+k]<(s5$medMz[j]*(1+ppm))){
      k=k+1
    }
    if(k>1){
      library_match[[s5$ID[j]]]=hmdb_df[(i_min+1):(i_min+k-1),]
    }
    else{library_match[[s5$ID[j]]]=hmdb_df[0,]}
    j=j+1
  }
  
  
  counts=unlist(lapply(library_match,nrow))
  num_of_library_match = cbind(ID=Mset$ID,counts)
  
  
  
  
  df = as.data.frame(num_of_library_match)
  
  df["library_match_formula"]=NA
  df["library_match_name"]=NA
  for(i in 1:length(library_match)){
    
    if(nrow(library_match[[i]]) == 0){next
    }else{
      df$library_match_formula[i] = paste(library_match[[i]]$MF, collapse = " | ")
      df$library_match_name[i] = paste(library_match[[i]]$Name, collapse = " | ")
    }
    
  }
  
  
  return (list(library_match_list = library_match,
               library_match_formula = df
               ))
}
## 
## Metaboanalyst_Statistic - Statistical analysis with MetaboAnalyst ####
Metaboanalyst_Statistic = function(Mset){
  library(MetaboAnalystR)
  
  MA_output = Mset$Data[,c("ID",  Mset$Cohort$sample_names)]
  MA_output = rbind(c("cohort",Mset$Cohort$sample_cohort),MA_output)

  write.csv(MA_output, file="MetaboAnalyst_file.csv", row.names=F)
  
  mSet <- InitDataObjects("pktable", "stat", FALSE)
  mSet<-Read.TextData(mSet, "MetaboAnalyst_file.csv", "colu", "disc")
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet);
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
  # mSet<-FilterVariable(mSet, "iqr", "F", 25)
  # mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", "a11", ratio=FALSE, ratioNum=20)
  mSet<-ANOVA.Anal(mSet, F, 0.05, "fisher")
  
  mSet<-PlotHeatMap(mSet, "full_", "png", 600, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview", T, T, NA, T, F)
  mSet<-my_PlotSubHeatMap(mSet, "top50_", "png", 600, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 50, "overview", T, T, T, F)
  
  ANOVA_file = "anova_posthoc.csv"
  ANOVA_raw <- read_csv(ANOVA_file)
  
  ANOVA_FDR = ANOVA_raw[,c("X1","FDR")]
  ANOVA_FDR$FDR=-log10(ANOVA_FDR$FDR)
  colnames(ANOVA_FDR)=c("ID","_log10_FDR")
  
  fn <- "t_test.csv"
  if (file.exists(fn)) file.remove(fn)
  fn <-"MetaboAnalyst_file.csv"
  if (file.exists(fn)) file.remove(fn)
  fn <-ANOVA_file
  if (file.exists(fn)) file.remove(fn)
  
  return(ANOVA_FDR)
}


### get PlotSubHeatMap function work ####
my_PlotSubHeatMap = function (mSetObj = NA, imgName, format = "png", dpi = 72, width = NA, 
                              dataOpt, scaleOpt, smplDist, clstDist, palette, method.nm, 
                              top.num, viewOpt, rowV = T, colV = T, border = T, grp.ave = F) 
{
  var.nms = colnames(mSetObj$dataSet$norm)
  if (top.num < length(var.nms)) {
    if (method.nm == "tanova") {
      if (mSetObj$dataSet$cls.num == 2) {
        if (is.null(mSetObj$analSet$tt)) {
          Ttests.Anal(mSetObj)
          #mSetObj <- .get.mSet(mSetObj)
        }
        #var.nms <- names(sort(mSetObj$analSet$tt$p.value))[1:top.num]
        var.nms <- names(sort(mSetObj$analSet$aov$p.value))[1:top.num]
      }
      else {
        if (is.null(mSetObj$analSet$aov)) {
          ANOVA.Anal(mSetObj)
          #mSetObj <- .get.mSet(mSetObj)
        }
        var.nms <- names(sort(mSetObj$analSet$aov$p.value))[1:top.num]
      }
    }
    else if (method.nm == "cor") {
      if (is.null(mSetObj$analSet$cor.res)) {
        Match.Pattern(mSetObj)
        #mSetObj <- .get.mSet(mSetObj)
      }
      cor.res <- mSetObj$analSet$cor.res
      ord.inx <- order(cor.res[, 3])
      cor.res <- cor.res[ord.inx, ]
      ord.inx <- order(cor.res[, 1])
      cor.res <- cor.res[ord.inx, ]
      var.nms <- rownames(cor.res)[1:top.num]
    }
    else if (method.nm == "vip") {
      if (is.null(mSetObj$analSet$plsda)) {
        PLSR.Anal(mSetObj)
        PLSDA.CV(mSetObj)
        #mSetObj <- .get.mSet(mSetObj)
      }
      vip.vars <- mSetObj$analSet$plsda$vip.mat[, 1]
      var.nms <- names(rev(sort(vip.vars)))[1:top.num]
    }
    else if (method.nm == "rf") {
      if (is.null(analSet$rf)) {
        RF.Anal(mSetObj)
        #mSetObj <- .get.mSet(mSetObj)
      }
      var.nms <- GetRFSigRowNames()[1:top.num]
    }
  }
  var.inx <- match(var.nms, colnames(mSetObj$dataSet$norm))
  PlotHeatMap(mSetObj, imgName, format, dpi, width, dataOpt, 
              scaleOpt, smplDist, clstDist, palette, viewOpt, rowV, 
              colV, var.inx, border, grp.ave)
}





## Summary_Mset ####
Summary_Mset = function(Mset){
  Mdata = Mset$Data[,c(3:7,14:ncol(Mset$Data))]
  High_blanks = Mset$High_blanks
  Mdata = merge(High_blanks,Mdata, all=T)
  HMDB = Mset$library_match$library_match_formula[,c("ID","library_match_formula","library_match_name")]
  Mdata = merge(HMDB,Mdata, all=T)
  ANOVA_FDR = Mset$Metaboanalyst_Statistic
  Mdata = merge(ANOVA_FDR,Mdata, all=T)
  
  return(Mdata)
}
# Function for network ######

## Form_node_list - Merge experiment and library nodes ####
Form_node_list = function(Mset)
{
  NodeSet = list()
  NodeSet[["Expe"]] = Mset$Data[,c("ID","medMz","medRt","formula")]
  NodeSet$Expe["category"]=1
  NodeSet$Expe["compound_name"]=NA
  colnames(NodeSet$Expe) = c("ID","mz","RT","MF", "category","compound_name")
  
  NodeSet[["Library"]] = Mset$Library
  NodeSet$Library = cbind(NodeSet$Library, RT=NA)
  colnames(NodeSet$Library) = c("ID","compound_name","MF","mz","category", "RT")
  NodeSet$Library$ID = 1:nrow(NodeSet$Library)+nrow(NodeSet$Expe)

  NodeSet$Library$category=0
  
  merge_node_list = rbind(NodeSet$Expe,NodeSet$Library )

  return(merge_node_list)
}

## Read_rule_table - for biotransformation rule and artifacts ####
Read_rule_table = function(rule_table_file = "biotransform.csv"){
  library(enviPat)
  data("isotopes")
  biotransform = read.csv(rule_table_file,stringsAsFactors = F)
  biotransform$Formula = check_chemform(isotopes,biotransform$Formula)$new_formula
  biotransform = biotransform[with(biotransform, order(mass)),]
  return(biotransform)
}
## Edge_list for biotransformation ####
Edge_biotransform = function(Mset, mass_abs = 0.001, mass_ppm = 5/10^6, read_from_csv=F)
{
  
  
  if(!read_from_csv){
    
    
  
  merge_node_list = Mset$NodeSet[with(Mset$NodeSet, order(mz)),]
  merge_nrow = nrow(merge_node_list)
  
  timer=Sys.time()
  {
    edge_ls = list()
    temp_mz_list=merge_node_list$mz
    
    for (k in 1:nrow(Mset$Biotransform)){
      temp_fg=Mset$Biotransform$mass[k]
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
      print(paste("Biotransform", Mset$Biotransform$category[k], nrow(temp_edge_list),"found."))
    }
  }
  
  print(Sys.time()-timer)
  edge_list = bind_rows(edge_ls)
  
  edge_list$linktype=Mset$Biotransform$Formula[edge_list$linktype]
  
  edge_list_sub = subset(edge_list, 
                         (edge_list$node1<=nrow(Mset$Data)|
                            edge_list$node2<=nrow(Mset$Data))&
                           edge_list$node1!=edge_list$node2
  )
  
  edge_list_sub["category"]=1
  
  write_csv(edge_list_sub,"edge_list_sub.txt")
  
  } else{
    
    edge_list_sub = read.csv("edge_list_sub.txt", na="NA", stringsAsFactors = F)
    
  }
  
  
  return(edge_list_sub)
}
### Scoring edge based on mass accuracy ####
Edge_score = function(Biotransform){
  if(nrow(Biotransform)>10000){
    edge_mzdif_FIT <- fitdist(as.numeric(Biotransform$mass_dif[base::sample(nrow(Biotransform),10000)]), "norm")    
  } else {
    edge_mzdif_FIT <- fitdist(as.numeric(Biotransform$mass_dif), "norm")    
  }
  
  #hist(Biotransform$mass_dif)
  
  plot(edge_mzdif_FIT)  
  Biotransform["edge_massdif_score"]=dnorm(Biotransform$mass_dif, edge_mzdif_FIT$estimate[1], edge_mzdif_FIT$estimate[2])
  Biotransform["edge_massdif_score"]=Biotransform["edge_massdif_score"]/max(Biotransform["edge_massdif_score"])
  return(Biotransform)
}
## Variance between peaks ####
Peak_variance = function(Mset, 
                         time_cutoff=0.1,
                         mass_cutoff=0,
                         correlation_cutoff = 0.2)
{
  df_raw = Mset$Data[,c("ID","medMz","medRt",Mset$Cohort$sample_names)]
  df_raw["mean_inten"]=rowMeans(df_raw[,Mset$Cohort$sample_names])
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
      
      temp_x = t(temp_df_raw[1,Mset$Cohort$sample_names])
      temp_y = t(temp_df_raw[,Mset$Cohort$sample_names])
      temp_df_raw$correlation = as.numeric(t(cor((temp_x), (temp_y))))
      
      temp_df_raw["node1"]=temp_df_raw$ID[1]
      temp_df_raw["node2"]=temp_df_raw$ID
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
      temp_edge_ls=temp_edge_ls[temp_edge_ls$correlation>correlation_cutoff,]
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
  print("High correlation peaks identified.")
  return(edge_ls)
}

### Edge_list for artifacts ####
Artifact_prediction = function(Mset, Peak_inten_correlation, search_ms_cutoff=0.001,read_from_csv=F)
{
  if(read_from_csv == F){
    #edge_ls_highcor = EdgeSet$Peak_inten_correlation
    edge_ls_highcor=Peak_inten_correlation
    edge_ls_highcor = edge_ls_highcor[with(edge_ls_highcor, order(mz_dif)),]
    junk_df = Mset$Artifacts
    
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
      ppm=5/10^6
      temp_df_oligo = temp_df[0,]
      
      df_raw = Mset$Data
      
      temp_edge_ls = data.frame(ratio=edge_ls_highcor$mz_dif/edge_ls_highcor$mz_node1)
      temp_edge_ls["rounding"]=round(temp_edge_ls[,1],digit=0)
      temp_edge_ls["dif"]=temp_edge_ls["rounding"]-temp_edge_ls[,1]
      temp_edge_ls = temp_edge_ls[(temp_edge_ls$rounding!=0)
                                  &(abs(temp_edge_ls$dif)<0.05),]
      
      for(i in 1:nrow(temp_edge_ls)){
        temp_data = edge_ls_highcor[rownames(temp_edge_ls)[i],]
        temp_mz1 = df_raw$medMz[which(df_raw$ID==temp_data$node1)]
        temp_mz2 = temp_mz1 + temp_data$mz_dif
        for(j in 2:10){
          #if charge data
          #if(abs(temp_mz1*j-(H_mass-e_mass)*mode*(j-1)-temp_mz2)<temp_mz2*ppm){
          #if neutral data
          if(abs(temp_mz1*j-temp_mz2)<search_ms_cutoff){
            temp_df_oligo[(nrow(temp_df_oligo)+1),]=c(temp_data,"oligomer",paste("x",j,sep=""), (temp_mz1*j-temp_mz2)/temp_mz1*10^6)
          }
        }
      }
    }
    rm(temp_edge_ls)
    temp_df_oligo
    nrow(temp_df_oligo)
    
    oligo_summary=table(temp_df_oligo$linktype)
    print(oligo_summary)
    
    #data$parent[data$ID==6]
    test_time = Sys.time()-test_time
    
    edge_ls_annotate=rbind(temp_df_isotope,temp_df_oligo)
    
    edge_ls_annotate_network = edge_ls_annotate[,c("node1","node2","linktype","mass_dif","category")]
    
    write_csv(edge_ls_annotate_network,"artifact_edge_list.txt")
  } else{
    
    edge_ls_annotate_network = read.csv("artifact_edge_list.txt",stringsAsFactors = F)
  }
  
  
  return(edge_ls_annotate_network)
}

## Merge_edgeset ####
Merge_edgeset = function(EdgeSet){
  
  edge_merge = rbind(EdgeSet$Artifacts,EdgeSet$Biotransform)
  Mdata = Mset$Data$log10_inten
  node1 = edge_merge$node1
  node2 = edge_merge$node2
  
  
  c1 = list()
  c2 = list()
  for(i in 1:nrow(edge_merge)){
    c1[[length(c1)+1]] = Mdata[node1[i]]
    c2[[length(c2)+1]] = Mdata[node2[i]]
  }
  
  edge_merge["node1_log10_inten"] = unlist(c1)
  edge_merge["node2_log10_inten"] = unlist(c2)
  edge_merge["edge_id"]=1:nrow(edge_merge)
  edge_merge["msr_inten_dif"] = edge_merge["node2_log10_inten"]-edge_merge["node1_log10_inten"]
  
  return(edge_merge)
}


## Network_prediction used to connect nodes to library and predict formula ####
Network_prediction = function(Mset, edge_list_sub, 
                              top_formula_n=2,
                              read_from_csv = F
                              )
{
  if(!read_from_csv){
  mnl=Mset$NodeSet
  # top_formula_n=3
  # edge_list_sub = EdgeSet$Merge
  #Initialize predict_formula from HMDB known formula
  {
    data_str = data.frame(id=as.numeric(),
                          formula=as.character(), 
                          steps=as.numeric(), 
                          parent=as.numeric(), 
                          is_metabolite = as.logical(),
                          score=as.numeric(), stringsAsFactors = F)
    #sf stands for summary_formula
    
    
    
    sf = lapply(1:nrow(mnl),function(i)(data_str))
    
    for(i in (1+nrow(Mset$Data)):nrow(mnl)){
      Initial_formula =sf[[i]]
      Initial_formula[1,]= list(i,mnl$MF[i],0,0,T,1)
      sf[[i]]=Initial_formula
      
    }
  }
  
  #while loop to Predict formula based on known formula and edgelist 
  
  nrow_experiment = nrow(Mset$Data)
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
      n=3
      for (n in 1: nrow(new_nodes_df)){
        #if(n%%1000==0){print(paste("Head_n =",n))}
        
        head = new_nodes_df$id[n]
        head_formula = new_nodes_df$formula[n]
        head_is_metabolite = new_nodes_df$is_metabolite[n]
        
        
        temp_edge_list=subset(edge_list_node1, edge_list_node1$node1==head)
        
        
        if(head <= nrow_experiment){
          #If head signal is < defined cutoff, then prevent it from propagating out, but it can still get formula from others.
          if(Mset$Data$mean_inten[head]< 2e4){next}
          #If head is an isotopic peak, then only look for isotopic peaks
          if(grepl("\\[",head_formula)){
            temp_edge_list = temp_edge_list[grepl("\\[",temp_edge_list$category),]
          }
          #If head is a metal or Boron or silicon adduct, then only look for adducts or skip
          if(!head_is_metabolite){
            #next
            temp_edge_list = temp_edge_list[temp_edge_list$category!=1, ]
          }
        }
        
        if(nrow(temp_edge_list)==0){next}
        
        head_score = sf[[head]]$score[match(head_formula,sf[[head]]$formula)]
        if(head_score==0){next}
        
        i=1
        for(i in 1:nrow(temp_edge_list)){
          tail=temp_edge_list$node2[i]
          if(tail>nrow_experiment){next}
          
          
          temp_fg = temp_edge_list$linktype[i]
          if(temp_fg==""){
            temp_formula = head_formula
          }else if (grepl("x",temp_fg)){
            fold = as.numeric(gsub("x","",temp_fg))
            temp_formula=my_calculate_formula(head_formula,head_formula,fold-1,Is_valid = T)
          }else {
            
            temp_formula=my_calculate_formula(head_formula,temp_fg,1,Is_valid = T)
          }
          
          #function return false if not valid formula
          if(is.logical(temp_formula))
          { temp_score=0
          temp_formula=paste(head_formula,temp_fg,sep="+")
          }else{
            temp_score = head_score*temp_edge_list$edge_massdif_score[i]
          }
          
          temp_parent = head
          temp_steps = step+1
          temp_is_metabolite = all(head_is_metabolite, temp_edge_list$category[i]==1)
          #Criteria to enter new entry into formula list
          #1. new formula
          temp = sf[[tail]]
          temp_subset=subset(temp, temp$formula==temp_formula)
          if(nrow(temp_subset)!=0){
            #2. If not a new metabolite status entry, then next
            if(any(temp_is_metabolite == temp_subset$is_metabolite)){
              next
            }
            #3. if not much higher scores, then next
            if(temp_score<=(1.2*max(temp_subset$score))){
              next
            }
          }
          

          #Enter new entry
          temp[nrow(temp)+1, ] = list(tail, temp_formula, temp_steps, temp_parent,temp_is_metabolite, temp_score)
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
        #if(n%%1000==0){print(paste("Tail_n =",n))}
        tail = new_nodes_df$id[n]
        tail_formula = new_nodes_df$formula[n]
        tail_is_metabolite = new_nodes_df$is_metabolite[n]
        
        if(tail <= nrow_experiment){
          #If tail signal is < defined cutoff, then prevent it from propagating out, but it can still get formula from others.
          if(Mset$Data$mean_inten[tail]<2e4){next}
          #If tail is an isotopic peak, then do not propagate
          if(grepl("\\[",tail_formula)){next}
          if(!tail_is_metabolite){
            next
            #temp_edge_list = temp_edge_list[temp_edge_list$category!=1, ]
          }
        }
        
        
        temp_edge_list=subset(edge_list_node2, edge_list_node2$node2==tail)
        
        if(nrow(temp_edge_list)==0){next}
        
        tail_score = sf[[tail]]$score[match(tail_formula,sf[[tail]]$formula)]
        
        if(tail_score==0){next}
        i=1
        for(i in 1:nrow(temp_edge_list)){
          head=temp_edge_list$node1[i]
          if(head>nrow_experiment){next}
          
          temp_fg = temp_edge_list$linktype[i]
          if(temp_fg==""){
            temp_formula = tail_formula
          }else if (grepl("x",temp_fg)){
            fold = as.numeric(gsub("x","",temp_fg))
            temp_formula=my_calculate_formula(head_formula,head_formula,-(fold-1)/fold,Is_valid = T)
            if(grepl(".", temp_formula)){temp_formula=F}
          }else {
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
          temp_is_metabolite = all(tail_is_metabolite, temp_edge_list$category[i]==1)
          #Criteria to enter new entry into formula list
          #1. new formula
          temp = sf[[head]]
          temp_subset=subset(temp, temp$formula==temp_formula)
          if(nrow(temp_subset)!=0){
            #2. If not a new metabolite status entry, then next
            if(any(temp_is_metabolite == temp_subset$is_metabolite)){
              next
            }
            #3. if not much higher scores, then next
            if(temp_score<=(1.2*max(temp_subset$score))){
              next
            }
          }
          #Enter new entry
          temp[nrow(temp)+1, ] = list(head, temp_formula, temp_steps, temp_parent, temp_is_metabolite, temp_score)
          temp = temp[with(temp, order(-score)),]
          sf[[head]]=temp
        }
      }
    }
    
    step=step+1
    
    
  }
  
  # Pruning formula  #
  {
    pred_formula = bind_rows(sf)
    pred_formula = pred_formula[!grepl("-",pred_formula$formula),]
    pred_formula = pred_formula[!duplicated(pred_formula[,c("id","formula")]),]
    
    pred_formula["mz"] = formula_mz(pred_formula$formula)
    pred_formula["measured_mz"]=NA
    
    
    for(i in 1:nrow(pred_formula)){
      pred_formula$measured_mz[i] = mnl$mz[mnl$ID==pred_formula$id[i]]
    }
    pred_formula["abs_dif"] = pred_formula["mz"]-pred_formula["measured_mz"]
    pred_formula["ppm_dif"] = pred_formula["abs_dif"]/pred_formula$measured_mz*1E6
    pred_formula_1 = pred_formula[abs(pred_formula$abs_dif)<0.001|abs(pred_formula$ppm_dif)<5,]
    
    
    pred_formula_1["rdbe"] = NA
    pred_formula_1$rdbe = formula_rdbe(pred_formula_1$formula)
    
    merge_formula = pred_formula_1[,1:6]
    sf = list()
    for(n in 1: max(merge_formula$id)){
      sf[[n]]=merge_formula[merge_formula$id==n,]
    }
    
    write_csv(merge_formula,"All_formula_predict.txt")
  }
  

  
  } else {
    
    merge_formula = read.csv("All_formula_predict.txt",stringsAsFactors = F)
    sf = list()
    for(n in 1: max(merge_formula$id)){
      sf[[n]]=merge_formula[merge_formula$id==n,]
    }
  }

  
  return(sf)
}



# Function for CPLEX ####
## Prepare_CPLEX parameter ####
Prepare_CPLEX = function(Mset, EdgeSet, read_from_csv = F){
  All_formula_predict=bind_rows(Mset[["NodeSet_network"]])
  raw_pred_formula = All_formula_predict
  raw_node_list = Mset$NodeSet
  raw_edge_list = EdgeSet$Merge
  
  #Clean up
  {
    pred_formula = raw_pred_formula[!grepl("-",raw_pred_formula$formula),]
    pred_formula = pred_formula[!duplicated(pred_formula[,c("id","formula")]),]
    
    
    lib_nodes = raw_node_list[raw_node_list$category==0,]
    lib_nodes_cutoff = nrow(raw_node_list)-nrow(lib_nodes)
    unknown_nodes = raw_node_list[raw_node_list$category!=0,]
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
    for(n in 1: max(merge_formula$id)){
      pred_formula_ls[[n]]=merge_formula[merge_formula$id==n,]
    }
  }
  
  ##Core codes
  
  #Construct constraint matrix 
  
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
      if(n%%10000==0){
        print(paste("n=",n,"elapsed="))
        print(Sys.time()-timer)
      }
      
      temp_edge = edge_list[n,]
      node_1 = temp_edge$node1
      node_2 = temp_edge$node2
      formula_1 = pred_formula_ls[[node_1]]
      formula_2 = pred_formula_ls[[node_2]]
      temp_fg = temp_edge$linktype
      

      temp_formula = formula_1$formula[2]
      for(temp_formula in unique(formula_1$formula)){
        temp_score = log10(temp_edge$edge_massdif_score+1e-10)+1

        #drop off edge where an isotopic peak reaches out to other
        #if(!grepl("\\[",temp_fg ) & grepl("\\[", temp_formula )){next}
        #modify score based on isotopic abundance of formula
        if(grepl("\\[",temp_fg )){
          calc_abun = isotopic_abundance(temp_formula, temp_fg)
          abun_ratio = calc_abun/10^temp_edge$msr_inten_dif
          score_modifier = dnorm(abun_ratio,1,0.2)/dnorm(1,1,0.2)
          #log score
          temp_score = temp_score + log10(score_modifier+1e-10)+1
          
          #linear score
          #temp_score = temp_score * score_modifier + temp_score * score_modifier * (1-temp_score * score_modifier)
        }
        
        #Assuming formula in node_1 is always smaller than node_2
        if(temp_fg==""){
          temp_formula_2=temp_formula
        }else{
          temp_formula_2 = my_calculate_formula(temp_formula, temp_fg)
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
            edge_info[[temp_i]] = list(edge_id=temp_edge$edge_id,
                                       edge_score=temp_score,
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
            edge_info[[temp_i]] = list(edge_id=temp_edge$edge_id,
                                       edge_score=temp_score,
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
      #edge_info_sum = CPLEXset$data$edge_info_sum
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
      #test1$edge_ilp_id[duplicated(test1[,c("ilp_index1","formula1")])
      edge_info_sum$edge_score[test1$edge_ilp_id[duplicated(test1[,c("ilp_index1","formula1")])]]=0
      
      
    }
    
    write_csv(triplet_df,"triplet_df.txt")
    write_csv(edge_info_sum,"edge_info_sum.txt")
    mat = simple_triplet_matrix(i=triplet_df$i,
                                j=triplet_df$j,
                                v=triplet_df$v)
  }else{
    triplet_df = read.csv("triplet_df.txt")
    edge_info_sum = read.csv("edge_info_sum.txt",stringsAsFactors = F)
    mat = simple_triplet_matrix(i=triplet_df$i,
                                j=triplet_df$j,
                                v=triplet_df$v)
    
  }
  
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
  return(CPLEX = list(data = CPLEX_data,
                      para = CPLEX_para)
  )
}

### Score_formula ####
Score_formula = function(CPLEXset)
{
  unknown_formula = CPLEXset$data$unknown_formula
  
  #when measured and calculated mass differ, score based on normal distirbution with mean=0 and sd=1e-3
  unknown_formula["msr_mass"] = Mset$Data$medMz[unknown_formula$id]
  unknown_formula["cal_mass"] = formula_mz(unknown_formula$formula)
  unknown_formula["msr_cal_mass_dif"] = unknown_formula["msr_mass"]-unknown_formula["cal_mass"]
  unknown_formula["Mass_score"] = dnorm(unknown_formula$msr_cal_mass_dif, 0, 1e-3)/dnorm(0, 0, 1e-3)
  
  #when rdbe < -1, penalty to each rdbe less
  unknown_formula["rdbe"] = formula_rdbe(unknown_formula$formula)
  unknown_formula["rdbe_score"] = sapply((0.1*(unknown_formula$rdbe)), min, 0)
  
  #when step is large, the likelihood of the formula is true decrease from its network score
  #Penalty happens when step > 5 on the existing score
  unknown_formula["step_score"] = sapply(-0.1*(unknown_formula$steps-5), min, 0) 
  
  #the mass score x step score evaluate from mass perspective how likely the formula fits the peak
  #the rdbe score penalizes unsaturation below -1
  #Each node should be non-positive, to avoid node formula without edge connection
  unknown_formula["cplex_score"] = log10(unknown_formula["Mass_score"])+log10(unknown_formula["score"])+unknown_formula["step_score"]+unknown_formula["rdbe_score"]
  
  # hist(unknown_formula$cplex_score)
  # length(unknown_formula$cplex_score[unknown_formula$cplex_score<1])
  
  return(unknown_formula)
}
### Score_edge_cplex ####
Score_edge_cplex = function(CPLEXset, edge_penalty = -0.5)
{
  edge_info_sum = CPLEXset$data$edge_info_sum
  edge_info_sum$edge_score = edge_info_sum$edge_score + edge_penalty
  unknown_formula = CPLEXset$data$unknown_formula
  
  test3 = edge_info_sum[edge_info_sum$ilp_index2<=nrow(unknown_formula)&
                          edge_info_sum$ilp_index1<=nrow(unknown_formula),]
  
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
  
  return(temp_edge_info_sum)
  
  
}
## Run_CPLEX ####
Run_CPLEX = function(CPLEXset, obj){

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
  tictoc::tic()
  return_code = mipoptCPLEX(env, prob)
  result_solution=solutionCPLEX(env, prob)
  
  print(paste(return_codeCPLEX(return_code),"-",
              status_codeCPLEX(env, getStatCPLEX(env, prob)),
              " - OBJ_value =", result_solution$objval))
  
  tictoc::toc()
  
  writeProbCPLEX(env, prob, "prob.lp")
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
  # 
  # CPLEX_x = result_solution$x
  # CPLEX_slack = result_solution$slack
  # 
  # 
  # if(write_to_csv){
  #   write_csv(as.data.frame(result_solution$x),"CPLEX_x.txt")
  #   write_csv(as.data.frame(result_solution$slack),"CPLEX_slack.txt")
  # }
  # 
  # } else{
  #   
  #   CPLEX_x = read.csv("CPLEX_x.txt")
  #   CPLEX_slack = read.csv("CPLEX_slack.txt")
  #   if(CPLEXset$CPLEX_para$nc!=nrow(CPLEX_x)){ print("CPLEX_x row number is incosistent with data!")}
  # }
  
  return(list(obj = obj, result_solution = result_solution))
}
## CPLEX_permutation ####
CPLEX_permutation = function(CPLEXset, n_pmt = 5, sd_rel_max = 0.5){
  
  edge_info_sum = CPLEXset$data$edge_info_sum
  unknown_formula = CPLEXset$data$unknown_formula
  obj = CPLEXset$Init_solution$obj
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

## CPLEX_screen_edge ####
CPLEX_screen_edge = function(CPLEXset, edge_penalty_range = seq(-.6, -0.9, by=-0.1)){

  solution_ls = list()
  #result_solution =1
  for(edge_penalty in edge_penalty_range){
    edge_info_sum = Score_edge_cplex(CPLEXset, edge_penalty = edge_penalty)
    temp_obj = c(CPLEXset$data$unknown_formula$cplex_score,
                            edge_info_sum$edge_score)
    result_solution = Run_CPLEX(CPLEXset, obj = temp_obj)
    solution_ls[[length(solution_ls)+1]] = result_solution
  }
  return(solution_ls)
  
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
# Function for graph ####
## Analysis of Subnetwork  ####
Subnetwork_analysis = function(g_sub, member_lb = 1, member_ub = 10)
{
  clu=components(g_sub)
  #subnetwork criteria 
  subnetwork = igraph::groups(clu)[table(clu$membership)<= member_ub & table(clu$membership)>= member_lb]
  
  g_subnetwork_list = lapply(subnetwork, make_ego_graph, graph=g_sub, order=diameter(g_sub), mode="all")
  for (i in 1:length(subnetwork)){
    #if(!any(merge_node_list$category[as.numeric(subnetwork[[i]])]>2)){next}
    plot(g_subnetwork_list[[i]][[1]],
         #vertex.color = 'white',
         vertex.label = vertex.attributes(g_subnetwork_list[[i]][[1]])$formula,
         #vertex.label = vertex.attributes(g_subnetwork_list[[i]][[1]])$mz,
         vertex.label.color = "red",
         vertex.label.cex = 1,
         #edge.color = 'black',
         edge.label = edge.attributes(g_subnetwork_list[[i]][[1]])$linktype,
         vertex.size = 10,
         edge.arrow.size = 2/length(vertex.attributes(g_subnetwork_list[[i]][[1]])$formula),
         main = paste("Subnetwork",names(subnetwork)[[i]])
    )
  }
}
## Analysis of Specific node ####
subgraph_specific_node = function(interested_node, g, step = 2)
{
  # interested_node = "2115"
  # g=g_sub
  interested_node = interested_node
  g.degree <- degree(g, mode = c("all"))
  g_intrest <- make_ego_graph(g, 
                              step, 
                              #1,
                              nodes = interested_node, 
                              mode = c("all"))[[1]]
  dists = distances(g_intrest, interested_node)
  colors <- c("black", "red", "orange", "blue", "dodgerblue", "cyan")
  # V(g_intrest)$color <- colors[dists+1]
  # png(filename=paste("Subnetwork of node ", interested_node,".png",sep=""),
  #     width = 2400, height=2400,
  #     res=300)

  plot(g_intrest,
       #vertex.color = 'white',
       vertex.label = vertex.attributes(g_intrest)$formula,
       #vertex.label = vertex.attributes(g_intrest)$medRt,
       vertex.label.color = "blue",
       vertex.label.cex = 1,
       #vertex.label.dist = 2,
       #vertex.size = vertex.attributes(g_intrest)$log10_inten,
       #edge.width = edge.attributes(g_intrest)$Confidence*2-2,
       # edge.color = color_palette[edge.attributes(g_intrest)$color],
       #edge.label = edge.attributes(g_intrest)$mz_dif,
       edge.label = edge.attributes(g_intrest)$linktype,
       edge.label.color = "red",
       edge.label.cex = 1,
       edge.arrow.size = 0.5,
       edge.arrow.width = 1,
       main = paste("Subnetwork of node", interested_node),
       layout = layout_nicely(g_intrest)
       #layout = layout_as_tree(g_intrest)
  )
  # dev.off()
}
# Debugging tools ####
## Trace formula history ####
Trace_step = function(query_id, unknown_node_CPLEX)
{
  #query_id=13
  df = unknown_node_CPLEX[0,]
  while(query_id <= max(unknown_node_CPLEX$ID)){
    df[nrow(df)+1,] = unknown_node_CPLEX[unknown_node_CPLEX$ID==query_id,]
    temp_parent = unknown_node_CPLEX$parent[unknown_node_CPLEX$ID==query_id]
    query_id = temp_parent
    
  }
  return(df)
}
#————————————————————————#####
# Main Codes ####
## Read files ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  Mset = list()
  Mset[["Library"]] = read_library("hmdb_unique.csv")
  Mset[["Biotransform"]]=Read_rule_table(rule_table_file = "biotransform.csv")
  Mset[["Artifacts"]]=Read_rule_table(rule_table_file = "artifacts.csv")
  
  datapath = ("C:/Users/Li Chen/Dropbox/temp/Melanie data/Brain")
  setwd(datapath)
  
  filename = c("brain.csv")
  
  Mset[["Raw_data"]] <- read_csv(filename)
  #Mset[["Raw_data"]] = Mset$Raw_data[base::sample(nrow(Mset$Raw_data),8000),]
  

  #write.csv(Mset[["Library"]], "HMDB_detected_nodes_clean.csv")
}

## Initialise ####
{
  Mset[["Global_parameter"]]=  list(mode = 1,
                                    normalized_to_col_median = F)
  Mset[["Cohort"]]=Cohort_Info(Mset)

  
  #Clean-up duplicate peaks 
  Mset[["Data"]] = Peak_cleanup(Mset,
                                ms_dif_ppm=1/10^6, 
                                rt_dif_min=0.01,
                                detection_limit=500)
  #View(Mset$Data)
  Mset[["ID"]]=Mset$Data$ID
}

## Feature generation ####
{
  #Identify peaks with high blanks
  Mset[["High_blanks"]]=High_blank(Mset, fold_cutoff = 2)
  #View(Mset$High_blanks)
  
  #library_match
  
  Mset[["library_match"]] = library_match(Mset, ppm=5/10^6)

  #Metaboanalyst_Statistic
  Mset[["Metaboanalyst_Statistic"]]=Metaboanalyst_Statistic(Mset)
  

  # output assigned formula
  Mset[["Summary"]] = Summary_Mset(Mset)
  write_csv(Mset$Summary, paste("Mdata",filename,sep="_"))

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
  EdgeSet[["Biotransform"]] = Edge_score(EdgeSet$Biotransform)
  
  
  EdgeSet[["Peak_inten_correlation"]] = Peak_variance(Mset,
                                                      time_cutoff=0.1,
                                                      mass_cutoff = 2e4,
                                                      correlation_cutoff = 0)
  EdgeSet[["Artifacts"]] = Artifact_prediction(Mset, 
                                               EdgeSet$Peak_inten_correlation, 
                                               search_ms_cutoff=0.001,
                                               read_from_csv = read_from_csv)
  EdgeSet[["Artifacts"]] = Edge_score(EdgeSet$Artifacts)
  
  
  EdgeSet[["Merge"]] = Merge_edgeset(EdgeSet)
  

  Mset[["NodeSet_network"]] = Network_prediction(Mset, 
                                                 EdgeSet$Merge, 
                                                 top_formula_n = 3,
                                                 read_from_csv = F)
  
  CPLEXset = Prepare_CPLEX(Mset, EdgeSet, read_from_csv = F)
  CPLEXset$data$unknown_formula = Score_formula(CPLEXset)
  
}

# Run CPLEX ####
{
  edge_info_sum = Score_edge_cplex(CPLEXset, edge_penalty = -log10(.5)-1)
  obj_cplex = c(CPLEXset$data$unknown_formula$cplex_score, edge_info_sum$edge_score)

  CPLEXset[["Init_solution"]] = list(Run_CPLEX(CPLEXset, obj_cplex))
  CPLEXset[["Init_solution2"]] = list(Run_CPLEX(CPLEXset, obj_cplex))

  #CPLEXset[["Screen_solution"]] = CPLEX_screen_edge(CPLEXset, edge_penalty_range = seq(-.6, -0.9, by=-0.1))
  #CPLEXset[["Pmt_solution"]] = CPLEX_permutation(CPLEXset, n_pmt = 5, sd_rel_max = 0.3)
}

# Read CPLEX result ####
{

  CPLEX_all_x = Read_CPLEX_result(CPLEXset$Init_solution)
  
  CPLEX_x = rowMeans(CPLEX_all_x,na.rm=T)
  
  unknown_nodes = CPLEXset$data$unknown_nodes
  unknown_formula = CPLEXset$data$unknown_formula
  
  unknown_formula["ILP_result"] = CPLEX_x[1:nrow(unknown_formula)]
  unknown_formula_CPLEX = unknown_formula[unknown_formula$ILP_result !=0,]
  print(paste("pred formula num =", nrow(unknown_formula_CPLEX)))
  
  unknown_node_CPLEX = merge(unknown_nodes,unknown_formula_CPLEX,by.x = "ID", by.y = "id",all=T)
  
  #edge_info_sum = CPLEXset$data$edge_info_sum
  
  edge_info_sum["ILP_result"] = CPLEX_x[(nrow(unknown_formula)+1):length(CPLEX_x)]
  edge_info_CPLEX = edge_info_sum[edge_info_sum$ILP_result!=0,]
  

  # output assigned formula
  Mdata = Mset$Summary
  formula = unknown_node_CPLEX[,c("ID","formula","is_metabolite")]
  Mdata2 = merge(formula, Mdata, all = T)
  write.csv(Mdata2, paste("Mdata",filename,sep="_"),row.names = F)
  
}


  # Helper function
{
  id = 233
  unknown_formula_id = unknown_formula[unknown_formula$id==id,]
  edge_list_id = EdgeSet$Merge[EdgeSet$Merge$node1==id | EdgeSet$Merge$node2==id,]
  edge_info_sum_id = edge_info_sum[edge_info_sum$edge_id %in% edge_list_id$edge_id,]
  
  Mset$NodeSet_network[[id]]
  
  sf = bind_rows(Mset$NodeSet_network)
  test = sf[duplicated(sf[,c("id","formula")],]
}

# Graphic analysis ####
{
  Graphset = list()
  merge_edge_list = EdgeSet$Merge[edge_info_CPLEX$edge_id,]
  merge_node_list = merge(Mset$NodeSet,unknown_node_CPLEX,all=T)
  for(i in 1:nrow(merge_node_list)){
    if(is.na(merge_node_list$formula[i]) & !is.na(merge_node_list$MF[i])){
      merge_node_list$formula[i] = merge_node_list$MF[i]
    }
  }
  
  merge_node_list["log10_inten"] = NA
  merge_node_list$log10_inten[Mset$Data$ID] = Mset$Data$log10_inten
  
  merge_node_list_1e5 = merge_node_list[merge_node_list$log10_inten>5 &
                                          !is.na(merge_node_list$log10_inten), ]
  nrow(merge_node_list_1e5[is.na(merge_node_list_1e5$formula),])
  
  colors <- c("white", "red", "orange", "yellow", "green")
  merge_node_list["color"] = colors[merge_node_list$category+1]
  
  g <- graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list, directed = T)
  #E(g)$color = colors[merge_edge_list$category+1]
  merge_node_list["degree"]=degree(g, mode = c("all"))
  
  g_sub = graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list[merge_node_list$ID %in% c(merge_edge_list$node1, merge_edge_list$node2),], directed = T)
  #E(g_sub)$color = colors[merge_edge_list$category+1]
  
  
  #Basic graph characteristics, distance, degree, betweeness
  {
    #distance
    farthest_vertices(g_sub) 
    #degree
    g.d = degree(g_sub, mode = c("all"))
    #Closeness
    #g.c = closeness(g_sub)
    #Betweenness
    g.b <- betweenness(g_sub, directed = T)
    #which.max(g.b)
    
    
    clu=components(g_sub)
    clu_member = table(clu$membership)
    #subnetwork criteria 
    mainnetwork = igraph::groups(clu)[table(clu$membership)>3000]
    
    g_sub_node_list = merge_node_list[merge_node_list$ID %in% c(merge_edge_list$node1, merge_edge_list$node2),]
    g_sub_main_node_list = g_sub_node_list[g_sub_node_list$ID %in% mainnetwork[[1]], ]
    
    nrow(g_sub_main_node_list[!is.na(g_sub_main_node_list$MF),])
    
    
    # plot(g_sub,
    #      vertex.label = NA,
    #      #edge.color = 'black',
    #      vertex.size = sqrt(degree(g_sub, mode = c("all"))),
    #      edge.arrow.size = 0.05,
    #      layout = layout_nicely(g_sub))
  }
}
# 
# png(filename=paste("full dataset plot",".png",sep=""),
#     width = 2400, height=2400,
#     res=600)
# plot(g_sub,
#      vertex.label = NA,
#      edge.color = 'red',
#      vertex.color = 'white',
#      vertex.size = sqrt(degree(g_sub, mode = c("all"))),
#      edge.arrow.size = 0.01,
#      
#      layout = layout_nicely(g_sub))
# dev.off()



{
  
  
  subgraph_specific_node("2115", g_sub,1)
  Subnetwork_analysis(g_sub, member_lb = 4, member_ub = 10)
  
}

  
  #Analyze the network/subgraph of specific node
  mainDir = dirname(rstudioapi::getSourceEditorContext()$path)
  subDir = "subgraph of specific node"
  dir.create(file.path(mainDir, subDir),showWarnings=F)
  setwd(file.path(mainDir, subDir))
  
  
  
  clu=components(g_sub)
  #subnetwork criteria 
  subnetwork = igraph::groups(clu)[table(clu$membership)<10000]
  
  g_subnetwork_list = lapply(subnetwork, make_ego_graph, graph=g_sub, order=diameter(g_sub), mode="all")
  
  
  find_node_in_subgraph=function(node_name, subnetwork){
    for(i in 1:length(subnetwork)){
      if(node_name%in% subnetwork[[i]]){return(i)}
    }
  }
  
  find_node_in_subgraph(331, subnetwork)
  
  for(i in 1:5){
    #Analyze the network/subgraph of specific node
    {
      temp = merge_node_list[i,]
      #temp = merge_node_list[612,]
      interested_node = paste(temp$ID)
      target_mz = round(temp$mz,digits=4)
      target_rt = round(temp$RT,digits=2)
      step = temp$steps+1
      
      target_subgraph = find_node_in_subgraph(as.character(i),subnetwork)
      #target_subgraph = find_node_in_subgraph(as.character(612),subnetwork)
      if(target_subgraph!=1){
        g_interest = make_ego_graph(g_sub, min(3,diameter(g_sub)), nodes = interested_node, mode = c("all"))[[1]]
      }else
      {g_interest <- make_ego_graph(g_sub,2, nodes = interested_node, mode = c("all"))[[1]]}
      
      #V(g_interest)$color <- colors[dists+1]
      
      plot(g_interest,
           vertex.label = vertex.attributes(g_interest)$formula,
           #vertex.label = vertex.attributes(g_interest)$mz,
           #vertex.label = vertex.attributes(g_interest)$ID,
           
           vertex.label.color = "black",
           vertex.label.cex = 1,
           #edge.color = 'black',
           edge.label = edge.attributes(g_interest)$linktype,
           vertex.size = 10,
           edge.arrow.size = .25,
           main = paste("mz=", target_mz," RT=",target_rt," formula=", temp$formula, sep="")
      )
    }
    
    png(filename=paste("mz=", target_mz," RT=",target_rt,".png",sep=""),
        width = 2400, height=2400,
        res=300)
    plot(g_interest,
         vertex.label = vertex.attributes(g_interest)$formula,
         #vertex.label = vertex.attributes(g_interest)$mz,
         #vertex.label = vertex.attributes(g_interest)$ID,
         vertex.label.color = "black",
         vertex.label.cex = 1,
         #edge.color = 'black',
         edge.label = edge.attributes(g_interest)$linktype,
         vertex.size = 10,
         edge.arrow.size = .5,
         main = paste("mz=", target_mz," RT=",target_rt," formula=", temp$formula, sep="")
    )
    dev.off()
    output_network_csv=merge_node_list[vertex.attributes(g_interest)$name,]
    write.csv(output_network_csv[,-5], paste("mz=", target_mz," RT=",target_rt,".csv",sep=""), row.names=F)
  }
  

  
  
  



  
  
## legacy code ####
# # select first few formula in the NodeSet_network to simplify 
#   {
#     all_formula = bind_rows(Mset[["NodeSet_network"]])
#     sf = list()
#     for(n in 1: max(all_formula$id)){
#       sf[[n]]=head(all_formula[all_formula$id==n,],3)
#     }
#     
#     All_formula_predict = bind_rows(sf)
#     
#     edge_merge = EdgeSet$Merge
#     hist(edge_merge$edge_massdif_score)
#     Mset[["NodeSet_network"]] = sf
#     
#     write_csv(all_formula,"All_formula_predict.txt")
#   }
# 
#   # HMDB #
#   {
# 
#     
#     unknown_node_CPLEX_HMDB = merge(unknown_node_CPLEX, df, all=T)
#     unknown_node_CPLEX_HMDB_dif = unknown_node_CPLEX_HMDB[unknown_node_CPLEX_HMDB$formula!=
#                                                             unknown_node_CPLEX_HMDB$Library_match_formula &
#                                                             (!is.na(unknown_node_CPLEX_HMDB$Library_match_formula)),]
#   }
#   #calculate natural abundance expectation
#   {
#     edge_info_isotope = EdgeSet$Merge
#     edge_info_isotope = edge_info_isotope[grepl("\\[",edge_info_isotope$linktype),]
#     edge_info_sum_isotope = edge_info_sum[grepl("\\[",edge_info_sum$formula2),]
#     
#     edge_info_sum_isotope["msr_iso_abun"] = NA
#     edge_info_sum_isotope["calc_iso_abun"] = NA
#     i=1
#     for(i in 1:nrow(edge_info_sum_isotope)){
#       edge_info_sum_isotope$msr_iso_abun[i] = 10^EdgeSet$Merge$msr_inten_dif[edge_info_sum_isotope$edge_id[i]]
#       edge_info_sum_isotope$calc_iso_abun[i] = (isotopic_abundance(edge_info_sum_isotope$formula1[i], 
#                                                                    EdgeSet$Merge$linktype[edge_info_sum_isotope$edge_id[i]]))
#     }
#     
#     edge_info_sum_isotope["iso_abun_ratio"] = edge_info_sum_isotope$msr_iso_abun / edge_info_sum_isotope$calc_iso_abun
#     edge_mzdif_FIT <- fitdist(edge_info_sum_isotope$iso_abun_ratio[edge_info_sum_isotope$iso_abun_ratio<1.5 &
#                                                                      edge_info_sum_isotope$iso_abun_ratio >0.5], "norm")    
#     plot(edge_mzdif_FIT)
#     
#     edge_info_sum_isotope["iso_abun_score2"]=dnorm(edge_info_sum_isotope$iso_abun_ratio, edge_mzdif_FIT$estimate[1], edge_mzdif_FIT$estimate[2])
#     edge_info_sum_isotope["iso_abun_score2"]=edge_info_sum_isotope["iso_abun_score2"]/max(edge_info_sum_isotope["iso_abun_score2"])
#     
#     
#   }
#   
#   
#   # Evaluate Xi's data from annotation
#   {
#     {
#       setwd("C:/Users/Li Chen/Desktop/Github local/myrepo/Merge/Xi_full_hmdb_all n=2")
#       
#       CPLEX_x2 = read_csv("CPLEX_x.txt")
#       raw_pred_formula = read.csv("All_formula_predict.txt",stringsAsFactors = F)
#       raw_node_list = Mset$NodeSet
#       pred_formula = raw_pred_formula[!grepl("-",raw_pred_formula$formula),]
#       pred_formula = pred_formula[!duplicated(pred_formula[,c("id","formula")]),]
#       lib_nodes = raw_node_list[raw_node_list$category==0,]
#       lib_nodes_cutoff = nrow(raw_node_list)-nrow(lib_nodes)
#       unknown_nodes = raw_node_list[raw_node_list$category!=0,]
#       unknown_nodes = unknown_nodes[unknown_nodes$ID %in% unique(pred_formula$id),]
#       num_unknown_nodes = nrow(unknown_nodes)
#       lib_formula = pred_formula[pred_formula$id %in% lib_nodes$ID,]
#       lib_formula = lib_formula[lib_formula$steps==0,]
#       unknown_formula = pred_formula[pred_formula$id %in% unknown_nodes$ID,]
#       unknown_formula["ILP_result"] = CPLEX_x2$x[1:nrow(unknown_formula)]
#       unknown_formula_CPLEX2 = unknown_formula[unknown_formula$ILP_result !=0,]
#       
#       unknown_node_CPLEX2 = merge(unknown_nodes,unknown_formula_CPLEX2,by.x = "ID", by.y = "id",all=T)
#     }
#     
#     unknown_node_CPLEX_merge = merge(unknown_node_CPLEX,unknown_node_CPLEX2, by="ID", all=T)
#     
#     unknown_node_CPLEX_merge_dif = unknown_node_CPLEX_merge[unknown_node_CPLEX_merge$formula.x!=
#                                                               unknown_node_CPLEX_merge$formula.y,]
#     
#     unknown_node_CPLEX_Xi = merge(unknown_node_CPLEX_merge, Mset$Data, by.x = "ID", by.y = "ID", all=T)
#     unknown_node_CPLEX_Xi = unknown_node_CPLEX_Xi[,c(1:ncol(unknown_node_CPLEX_merge), 28,29)]
#     
#     unknown_node_CPLEX_Xi_met = unknown_node_CPLEX_Xi[unknown_node_CPLEX_Xi$isotopeLabel=="'Metabolite'",]
#     unknown_node_CPLEX_Xi_met_dif = unknown_node_CPLEX_Xi_met[unknown_node_CPLEX_Xi_met$formula.x!=
#                                                                 unknown_node_CPLEX_Xi_met$formula.y,]
#     
#   }
#   {
#     wl_result = read_csv("WL_data_190405.csv")
#     merge_result = merge(unknown_node_CPLEX[,c("ID","formula","is_metabolite")],wl_result, by.x="ID", by.y = "id", all =T)
#     merge_result$formula.y = check_chemform(isotopes, merge_result$formula.y)$new_formula
#   }  
#   
#   # signal > e5
#   {
#     Mdata = Mset$Data
#     e5_id = Mdata$ID[Mdata$log10_inten>=5]
#     merge_result_e5 = merge_result[merge_result$ID %in% e5_id, ]
#     e6_id = Mdata$ID[Mdata$log10_inten>=6]
#     merge_result_e6 = merge_result[merge_result$ID %in% e6_id, ]
#     e7_id = Mdata$ID[Mdata$log10_inten>=7]
#     merge_result_e7 = merge_result[merge_result$ID %in% e7_id, ]
#   }
#   
#   # Compare formula #
#   {
#     
#     merge_result_with_formula = merge_result[merge_result$formula.y!="[]",]
#     merge_result_with_formula_correct = merge_result_with_formula[merge_result_with_formula$formula.x==merge_result_with_formula$formula.y &
#                                                                     (!is.na(merge_result_with_formula$formula.x)),]
#     merge_result_with_formula_dif = merge_result_with_formula[merge_result_with_formula$formula.x!=merge_result_with_formula$formula.y |
#                                                                 (is.na(merge_result_with_formula$formula.x)),]
#     nrow(merge_result_with_formula_correct)/nrow(merge_result_with_formula)
#   }
#   
#   