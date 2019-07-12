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
  # library(cplexAPI)
  library(MetaboAnalystR)
  library(pracma)
  library(tictoc)
  # install.packages("janitor")
  library(janitor)
  #devtools::install_github("LiChenPU/Formula_manipulation")
  library(lc8)
  library(profvis)
  
  # setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

# Fucntions ####
# Function for parsing#### 
## read_library ####
read_library = function(library_file){
  data(isotopes)
  hmdb_lib = read_csv(library_file)
  hmdb_lib$MF = check_chemform(isotopes, hmdb_lib$MF)$new_formula
  hmdb_lib$MF = sapply(hmdb_lib$MF, my_calculate_formula,"C1")
  hmdb_lib$MF = sapply(hmdb_lib$MF, my_calculate_formula,"C1",-1)
  hmdb_lib$Exact_Mass = formula_mz(hmdb_lib$MF)
  hmdb_lib["rdbe"]=formula_rdbe(hmdb_lib$MF)
  return(hmdb_lib)
}
## Read_rule_table - for biotransformation rule and artifacts ####
Read_rule_table = function(rule_table_file = "biotransform.csv"){
  data("isotopes")
  biotransform = read.csv(rule_table_file,stringsAsFactors = F)
  for(i in 1: nrow(biotransform)){
    if(biotransform$Formula[i]==""){next}
    biotransform$Formula[i] = check_chemform(isotopes,biotransform$Formula[i])$new_formula
    biotransform$Formula[i] = my_calculate_formula(biotransform$Formula[i], "C1",Is_valid = F)
    biotransform$Formula[i] = my_calculate_formula(biotransform$Formula[i], "C1", -1 ,Is_valid = F)
    biotransform$mass[i] = formula_mz(biotransform$Formula[i])
  }
  
  biotransform = biotransform[with(biotransform, order(mass)),]
  return(biotransform)
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
  if(length(Mset$Cohort$sample_cohort) != length(Mset$Cohort$sample_names))
  {print("Warning! cohort number does not match sample number.")}
  
  return(list("sample_names"=sample_names,"blank_names"=blank_names, "sample_cohort"=sample_cohort))
}

## Peak_cleanup - Clean up duplicate peaks from peak picking ####
Peak_cleanup = function(Mset, 
                        ms_dif_ppm=1/10^6, 
                        rt_dif_min=0.01,
                        detection_limit=500
)
{
  raw = Mset$Raw_data
  #raw = raw[complete.cases(raw[, (1+which(colnames(raw)=="parent")):ncol(raw)]),]
  colnames(raw)[colnames(raw)=="groupId"] = "ID"
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ion_mode = Mset$Global_parameter$mode
  raw$medMz = raw$medMz*abs(ion_mode) - (H_mass-e_mass)*ion_mode
  
  
  
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
    
    
    ncol_raw = ncol(raw)
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
      s3$medMz[k_min]=mean(s3$medMz[k_min:(k_max-1)], na.rm=TRUE)
      s3$medRt[k_min]=mean(s3$medRt[k_min:(k_max-1)], na.rm=TRUE)
      # s3$goodPeakCount[k_min]=max(s3$goodPeakCount[k_min:(k_max-1)], na.rm=TRUE)
      s3$ID[k_min]=min(s3$ID[k_min:(k_max-1)], na.rm=TRUE)
      temp = s3[k_min:(k_max-1),14:ncol_raw]
      temp[1,] = apply(temp, 2, function(x){
        if(any(!is.na(x))){
          return(max(x, na.rm = T))
        } else {
          return(NA)
        }
      })
      s3[k_min, 14:ncol_raw] = temp[1,]
      # 
      # for (n in 14:ncol(raw)){
      #   if(any(!is.na(s3[k_min:(k_max-1),n]))){
      #     s3[k_min,n]=max(s3[k_min:(k_max-1),n], na.rm=TRUE)
      #   }
      # }
    }
  }
  
  #intermediate files, replace below detection number to random small number
  {
    
    s4=s3
    # s4[,4:ncol(s4)][s4[,4:ncol(s4)]<detection_limit]=sample(1:detection_limit, 
    #                                                         size=sum(s4[,4:ncol(s4)]<detection_limit), 
    #                                                         replace=T)
    
    s4["mean_inten"]=NA
    s4["mean_inten"]=rowMeans(s4[,Mset$Cohort$sample_names], na.rm = T)
    s4$flag[s4$mean_inten<detection_limit]=F
  }
  
  
  s5 = s4[s4$flag, 1:ncol(raw)]
  s5 = s5[with(s5, order(ID)),]
  s5$ID = 1:nrow(s5)
  
  s5["mean_inten"]=rowMeans(s5[,Mset$Cohort$sample_names], na.rm = T)
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
  
  
  Mdata = Mset$library_match$library_match_formula
  Mdata["Merge_ID"] = with(Mdata, paste(ID, library_match_formula, library_match_name, sep="_" ))
  MA_output = Mset$Data[,c("ID",  Mset$Cohort$sample_names)]
  MA_output$ID = Mdata$Merge_ID
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
  mSet<-ANOVA.Anal(mSet, F, 0.5, "fisher")
  # gc()
  # if(ncol(mSet$dataSet$norm) > 15000){
  #   mSet<-my_PlotSubHeatMap(mSet, paste(gsub(".csv","", filename),"top15000_"), "png", 600, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 15000, "overview", T, T, T, F)
  # } else{
  #   mSet<-PlotHeatMap(mSet, paste(gsub(".csv","", filename),"full_"), "png", 600, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview", T, T, NA, T, F)
  # }
  # gc()
  # mSet<-my_PlotSubHeatMap(mSet, paste(gsub(".csv","", filename),"top50_"), "png", 600, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 50, "overview", T, T, T, F)
  # gc()
  # 
  
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
  
  gc()
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
  if(!is.null(Mset$Metaboanalyst_Statistic)){
    ANOVA_FDR = Mset$Metaboanalyst_Statistic
    ANOVA_FDR$ID = as.numeric(gsub("_.*","", ANOVA_FDR$ID))
    Mdata = merge(ANOVA_FDR,Mdata, all=T)
  } else{
    ID = Mdata$ID
    log10_FDR = rep(NA, length(ID))
    ANOVA_FDR = cbind(ID, log10_FDR)
    Mdata = merge(ANOVA_FDR,Mdata, all=T)
  }
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
  NodeSet$Expe["rdbe"]=NA
  colnames(NodeSet$Expe) = c("ID","mz","RT","MF", "category","compound_name","rdbe")
  
  NodeSet[["Library"]] = Mset$Library
  NodeSet$Library = cbind(NodeSet$Library, RT= -1)
  colnames(NodeSet$Library) = c("ID","compound_name","MF","mz","category","rdbe","RT")
  NodeSet$Library$ID = 1:nrow(NodeSet$Library)+nrow(NodeSet$Expe)
  NodeSet$Library$category=0
  
  merge_node_list = rbind(NodeSet$Expe,NodeSet$Library)
  
  NodeSet[["Adduct"]] = Mset$Artifacts[Mset$Artifacts$category=="adduct",]
  NodeSet$Adduct["RT"] = -1
  NodeSet$Adduct["ID"] = 1:nrow(NodeSet$Adduct)+nrow(merge_node_list)
  NodeSet$Adduct$category = -1
  NodeSet$Adduct = NodeSet$Adduct[,c("ID", "mass", "RT", "Formula", "category", "Symbol", "rdbe")]
  colnames(NodeSet$Adduct) = c("ID","mz","RT","MF", "category","compound_name","rdbe")
  
  merge_node_list = rbind(merge_node_list,NodeSet$Adduct)
  
  return(merge_node_list)
}

## Edge_biotransform - generate Edge_list for biotransformation ####
Edge_biotransform = function(Mset, mass_abs = 0.001, mass_ppm = 5, read_from_csv=F)
{
  if(!read_from_csv){
    mass_ppm = mass_ppm/10^6
    merge_node_list = Mset$NodeSet[Mset$NodeSet$category!=-1,]
    merge_node_list = merge_node_list[with(merge_node_list, order(mz)),]
    merge_nrow = nrow(merge_node_list)
    
    timer=Sys.time()
    {
      edge_ls = list()
      temp_mz_list=merge_node_list$mz
      
      for (k in 1:nrow(Mset$Biotransform)){
        temp_fg=Mset$Biotransform$mass[k]
        temp_direction = Mset$Biotransform$direction[k]
        temp_rdbe = Mset$Biotransform$rdbe[k]
        i=j=1
        temp_edge_list = data.frame(node1=as.numeric(), 
                                    node2=as.numeric(), 
                                    linktype=as.numeric(), 
                                    mass_dif=as.numeric(), 
                                    direction = as.numeric(),
                                    rdbe = as.numeric())
        while(i<=merge_nrow){
          temp_ms=0
          mass_tol = max(temp_mz_list[i]*mass_ppm,0.001)
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
    # write_csv(edge_list_sub,"edge_list_sub.txt")
    
  } else{
    edge_list_sub = read.csv("edge_list_sub.txt", na="NA", stringsAsFactors = F)
  }
  return(edge_list_sub)
}
## Check_sys_measure_error - Check systematic error ####
Check_sys_measure_error = function(Biotransform, inten_threshold=1e5){
  
  Biotransform = EdgeSet$Biotransform
  Biotransform = Biotransform[Biotransform$linktype=="",]
  # Biotransform = Biotransform[Biotransform$node1>nrow(Mset$Data)|Biotransform$node2>nrow(Mset$Data),]
  high_inten_node = Mset$Data$ID[Mset$Data$mean_inten>inten_threshold]
  Biotransform = Biotransform[Biotransform$node1 %in% high_inten_node | Biotransform$node2 %in% high_inten_node, ]
  
  A_col_1 = - Mset$NodeSet$mz[Biotransform$node1]*as.numeric(Biotransform$node1<=nrow(Mset$Data))/Mset$NodeSet$mz[Biotransform$node2]
  A_col_2 = Mset$NodeSet$mz[Biotransform$node2]*as.numeric(Biotransform$node2<=nrow(Mset$Data))/Mset$NodeSet$mz[Biotransform$node2]
  A_col_3 = - as.numeric(Biotransform$node1<=nrow(Mset$Data))/Mset$NodeSet$mz[Biotransform$node2]
  A_col_4 = as.numeric(Biotransform$node2<=nrow(Mset$Data))/Mset$NodeSet$mz[Biotransform$node2]
  
  A = cbind(A_col_1,A_col_2,A_col_3,A_col_4)
  b = - Biotransform$mass_dif
  C = matrix(c(1,0,-1,0,0,1,0,-1), nrow = 2)
  d = matrix(c(0,0), nrow=2)
  x <- lsqlin(A, b, C,d)
  ppm_adjust = x[1]
  abs_adjust = x[3]/10^6
  
  
  if(abs(ppm_adjust)>0.5 | abs(abs_adjust)> 1e-4){
    print(paste("expect systematic measurement error, ppm shift =", ppm_adjust, "and abs shift = ", abs_adjust ))
    return(adjust = c(ppm_adjust=ppm_adjust, abs_adjust=abs_adjust))
  } else{
    print(paste("Do not expect systematic measurement error, ppm shift =", ppm_adjust, "and abs shift = ", abs_adjust ))
    return(adjust = c(ppm_adjust=ppm_adjust, abs_adjust=abs_adjust))
  }
  
  
  # Check if distribution is gaussian
  # shapiro.test(EdgeSet$Biotransform$mass_dif[base::sample(nrow(EdgeSet$Biotransform),5000)])
  # 
  # #install.packages("mclust")
  # library(mclust)
  # mc <- Mclust(Biotransform$mass_dif, G=3)
  # mc$parameters
  # a=mc$z
  # 
  # #install.packages("flexmix")
  # library(flexmix)
  # fl = flexmix(Biotransform_0$mass_dif~1, k=4)
  # a = fl@posterior
  # 
}

### Edge_score - Scoring edge based on mass accuracy ####
Edge_score = function(Biotransform){
  #Biotransform = EdgeSet$Biotransform
  if(nrow(Biotransform)>10000){
    edge_mzdif_FIT <- fitdist(as.numeric(Biotransform$mass_dif[base::sample(nrow(Biotransform),10000)]), "norm")    
  } else {
    edge_mzdif_FIT <- fitdist(as.numeric(Biotransform$mass_dif), "norm")    
  }
  
  # plot(edge_mzdif_FIT)  
  summary(edge_mzdif_FIT)
  
  Biotransform["edge_massdif_score"]=dnorm(Biotransform$mass_dif, 0, edge_mzdif_FIT$estimate[2])
  Biotransform["edge_massdif_score"]=Biotransform["edge_massdif_score"]/max(Biotransform["edge_massdif_score"])
  return(Biotransform)
}
## Peak_variance - Variance between peaks ####
Peak_variance = function(Mset, 
                         time_cutoff=0.1,
                         TIC_cutoff=10000,
                         correlation_cutoff = -1)
{
  df_raw = Mset$Data[,c("ID","medMz","medRt",Mset$Cohort$sample_names)]
  df_raw["mean_inten"]=rowMeans(df_raw[,Mset$Cohort$sample_names], na.rm = T)
  df_raw["log10_inten"]=log10(df_raw$mean_inten)
  
  {
    df_raw = df_raw[with(df_raw, order(medRt)),]
    df_raw_medRt = df_raw$medRt
    i = 1
    
    edge_list = list()
    timer = Sys.time()
    while(i!= nrow(df_raw)){
      if(i %% 1000 ==0 ){
        print(paste("nrow",i,"elapsed="))
        print((Sys.time()-timer))
      }
      
      if(df_raw$mean_inten[i]<TIC_cutoff){
        i = i+1
        next
      }
      
      temp_t = df_raw_medRt[i]
      
      temp_t_start = max(temp_t-time_cutoff, df_raw_medRt[1])
      temp_t_end = min(temp_t+time_cutoff, df_raw_medRt[nrow(df_raw)])
      
      i_start = i_end = i
      
      
      while(i_start >= 1 & df_raw_medRt[i_start] > temp_t_start){
        i_start = i_start-1
      }
      i_start = i_start+1
      
      while(i_start <= nrow(df_raw) & df_raw_medRt[i_end] < temp_t_end){
        i_end = i_end + 1
      }
      i_end = i_end - 1
      
      temp_df_raw = df_raw[i_start:i_end,]
      
      temp_df_raw$time_dif=temp_df_raw$medRt-df_raw_medRt[i]
      temp_df_raw$mz_dif = round(temp_df_raw$medMz-df_raw$medMz[i], digits=5)
      
      
      temp_df_raw$correlation = 1
      # Ignore co-variance between peaks
      # temp_x = t(df_raw[i,Mset$Cohort$sample_names])
      # temp_y = t(temp_df_raw[,Mset$Cohort$sample_names])
      # temp_df_raw$correlation = as.numeric(t(cor(temp_x, temp_y)))
      
      # temp_x = c(0.1,1,1.1,1.2)
      # temp_y = c(0.01,0.05,0.02,0.1)
      # cor(temp_x, temp_y)
      
      temp_df_raw["node1"]=df_raw$ID[i]
      temp_df_raw["node2"]=temp_df_raw$ID
      
      temp_edge_ls=temp_df_raw[,c("node1","node2", "correlation")]
      
      # cor(temp_x, temp_y, method = "kendall")
      # cor(temp_x, temp_y, method = "spearman")
      # cosine(tem#p_x, temp_y)
      
      temp_edge_ls=temp_edge_ls[temp_edge_ls$correlation>correlation_cutoff,]
      edge_list[[i]]=temp_edge_ls
      
      i = i+1
    }
  }
  
  edge_ls = bind_rows(edge_list)
  edge_ls=edge_ls[edge_ls$node1!=edge_ls$node2,]
  # edge_ls["mz_node1"] = Mset$Data$medMz[edge_ls$node1]
  # edge_ls["mz_node2"] = Mset$Data$medMz[edge_ls$node2]
  edge_ls["mz_dif"] = Mset$Data$medMz[edge_ls$node2]-Mset$Data$medMz[edge_ls$node1]
  
  # calculate mass difference between measured peak and adducts
  {
    adduct_set = Mset$NodeSet[Mset$NodeSet$category == -1,]
    measure_set = Mset$NodeSet[Mset$NodeSet$category == 1,]
    temp = outer(adduct_set$mz, measure_set$mz, "-")
    rownames(temp) = adduct_set$ID
    colnames(temp) = measure_set$ID
    temp_gather = gather(as.data.frame(temp), key = "node1", value="mz_dif")
    temp_gather["node1"] = as.numeric(temp_gather$node1)
    temp_gather["node2"] = adduct_set$ID
    # temp_gather["mz_node1"] = adduct_set$mz
    # temp_gather["mz_node2"] = Mset$Data$medMz[temp_gather$node2]
    temp_gather["correlation"] = 1
  }
  edge_ls = rbind(edge_ls, temp_gather)
  
  
  #Switch directionality
  {
    edge_ls = edge_ls[with(edge_ls, order(mz_dif)),]
    temp_data = edge_ls[edge_ls$mz_dif<0,]
    
    temp_node=temp_data$node1
    temp_data$node1=temp_data$node2
    temp_data$node2=temp_node
    temp_data$mz_dif = -temp_data$mz_dif
    
    edge_ls[edge_ls$mz_dif<0,] = temp_data 
    rm(temp_data)
  }
  
  
  # Remove duplicated rows
  {
    # edge_ls = edge_ls[!duplicated(edge_ls[,c("node1","node2")]),]
    edge_ls = edge_ls[with(edge_ls, order(node1, node2)),]
    keep_row = rep(T,nrow(edge_ls))
    node1 = edge_ls$node1
    node2 = edge_ls$node2
    for(i in 2:nrow(edge_ls)){
      if(node1[i] == node1[i-1]){
        if(node2[i] == node2[i-1]){
          keep_row[i]=F
        }
      }
    }
    edge_ls = edge_ls[keep_row, ]
  }
  
  edge_ls = edge_ls[with(edge_ls, order(mz_dif)),]
  edge_ls["mz_node1"] = Mset$NodeSet$mz[edge_ls$node1]
  edge_ls["mz_node2"] = Mset$NodeSet$mz[edge_ls$node2]
  # edge_ls["mz_dif"] = Mset$NodeSet$mz[edge_ls$node2] -Mset$NodeSet$mz[edge_ls$node1]
  # edge_ls["log10_inten_node1"] = Mset$Data$log10_inten[edge_ls$node1]
  # edge_ls["log10_inten_node2"] = Mset$Data$log10_inten[edge_ls$node2]
  # edge_ls["log10_inten_ratio"] = edge_ls["log10_inten_node2"]-edge_ls["log10_inten_node1"]
  # edge_ls["time_dif"] = Mset$Data$medRt[edge_ls$node2] - Mset$Data$medRt[edge_ls$node1]
  
  print("High correlation peaks identified.")
  return(edge_ls)
}

### Artifact_prediction - Edge_list for artifacts ####
Artifact_prediction = function(Mset, Peak_inten_correlation, search_ms_cutoff=0.002, search_ppm_cutoff = 10, read_from_csv=F)
{
  if(read_from_csv == F){
    # edge_ls_highcor = EdgeSet$Peak_inten_correlation
    edge_ls_highcor = Peak_inten_correlation
    edge_ls_highcor = edge_ls_highcor[with(edge_ls_highcor, order(mz_dif)),]
    junk_df = Mset$Artifacts
    
    i=j=1
    temp_ls = list()
    temp_df = edge_ls_highcor[1,]
    temp_df["category"]=NA
    temp_df["linktype"]=NA
    temp_df["direction"]=NA
    temp_df["rdbe"]=NA
    temp_df["mass_dif"]=double()
    
    temp_df = temp_df[0,]
    
    edge_ls_highcor_mz_dif = edge_ls_highcor$mz_dif
    timer = Sys.time()
    while (i <= nrow(edge_ls_highcor)){
      if(i%%1000000==0){print(paste(i, "edges screened."))
        print(Sys.time()-timer)}
      
      search_cutoff = max(search_ms_cutoff,search_ppm_cutoff*edge_ls_highcor$mz_node1[i]/10^6)
      
      
      if(edge_ls_highcor_mz_dif[i]<(junk_df$mass[j]-search_cutoff)){
        i=i+1
        next
      }
      
      if(edge_ls_highcor_mz_dif[i]<(junk_df$mass[j]+search_cutoff)){
        temp_df[(nrow(temp_df)+1),]=c(edge_ls_highcor[i,],
                                      as.character(junk_df$Symbol[j]), 
                                      junk_df$Formula[j], 
                                      junk_df$direction[j],
                                      junk_df$rdbe[j],
                                      (edge_ls_highcor_mz_dif[i]-junk_df$mass[j])/edge_ls_highcor$mz_node2[i]*10^6
        )
        
      }
      if(edge_ls_highcor_mz_dif[i]>(junk_df$mass[j]+10*search_ms_cutoff)){
        temp_ls[[j]]=temp_df
        temp_df = temp_df[0,]
        j=j+1
        if(j>nrow(junk_df)){break}
        #search_cutoff = 10 * initial_FIT$estimate[2]/1000 
        while(edge_ls_highcor_mz_dif[i]>(junk_df$mass[j]-10*search_ms_cutoff)){i=i-1}
        next
      }
      i=i+1
    }
    
    
    
    temp_df_isotope = bind_rows(temp_ls)
    
    junk_summary=table(temp_df_isotope$category)
    print(junk_summary)
    
    #Oligomers. Note that it is indistinguishable between 2-charge-parent pair and parent-dimer pair
    
    test_time = Sys.time()
    {
      
      temp_df_oligo = temp_df[0,]
      
      nodeset = Mset$NodeSet[Mset$NodeSet$category!=0,]
      
      temp_edge_ls = data.frame(ratio=edge_ls_highcor_mz_dif/edge_ls_highcor$mz_node1)
      temp_edge_ls["rounding"]=round(temp_edge_ls[,1],digit=0)
      temp_edge_ls["dif"]=temp_edge_ls["rounding"]-temp_edge_ls[,1]
      temp_edge_ls = temp_edge_ls[(temp_edge_ls$rounding!=0)
                                  &(abs(temp_edge_ls$dif)<0.05),]
      
      for(i in 1:nrow(temp_edge_ls)){
        temp_data = edge_ls_highcor[rownames(temp_edge_ls)[i],]
        temp_mz1 = nodeset$mz[which(nodeset$ID==temp_data$node1)]
        temp_mz2 = temp_mz1 + temp_data$mz_dif
        search_cutoff = max(search_ms_cutoff,search_ppm_cutoff*temp_mz2/10^6)
        for(j in 2:10){
          if(abs(temp_mz1*j-temp_mz2)<search_cutoff){
            temp_df_oligo[(nrow(temp_df_oligo)+1),]=c(temp_data,
                                                      "oligomer",
                                                      paste("x",j,sep=""), 
                                                      0,
                                                      paste("x",j,sep=""),
                                                      (temp_mz1*j-temp_mz2)/temp_mz2*10^6)
          }
        }
      }
    }
    rm(temp_edge_ls)
    
    nrow(temp_df_oligo)
    
    oligo_summary=table(temp_df_oligo$linktype)
    print(oligo_summary)
    
    #data$parent[data$ID==6]
    test_time = Sys.time()-test_time
    
    edge_ls_annotate=rbind(temp_df_isotope,temp_df_oligo)
    
    
    
    edge_ls_annotate_network = edge_ls_annotate[,c("node1","node2","linktype","mass_dif","category","direction","rdbe")]
    
    # write_csv(edge_ls_annotate_network,"artifact_edge_list.txt")
  } 
  else{
    edge_ls_annotate_network = read.csv("artifact_edge_list.txt",stringsAsFactors = F)
  }
  return(edge_ls_annotate_network)
}

### Hetero_dimer - Edge_list for hetero_dimer ####
Hetero_dimer = function(Peak_inten_correlation)
{
  # e = EdgeSet$Peak_inten_correlation
  e = Peak_inten_correlation
  e2 = e[e$node1 %in% Mset$Data$ID & e$node2 %in% Mset$Data$ID, ]
  e3 = e2[e2$mz_dif > min(e2$mz_node1) & e2$mz_dif < max(e2$mz_node1),]
  e3_list = split(e3, e3$node2)
  hetero_dimer_ls = list()
  for(i in 1: length(e3_list)){
    temp_e = e3_list[[i]]
    temp_matrix = outer(temp_e$mz_node1, temp_e$mz_node1, FUN = "+")  # mz_node1 and mz_dif are the same.
    temp_matrix = (temp_matrix - temp_e$mz_node2[1])/temp_e$mz_node2[1] * 10^6
    temp_index = which(abs(temp_matrix) < 5, arr.ind = T)
    if(length(temp_index)>0){
      temp_ppm = temp_matrix[temp_index]
      temp_node_1 = temp_e$node1[temp_index[,1]]
      linktype = temp_e$node1[temp_index[,2]]
      temp_df = data.frame(node1 = temp_node_1, linktype = linktype, node2 = temp_e$node2[1], mass_dif = temp_ppm)
      hetero_dimer_ls[[length(hetero_dimer_ls)+1]] = temp_df
    }
  }
  hetero_dimer_df = bind_rows(hetero_dimer_ls)
  hetero_dimer_df["category"]="Heterodimer"
  hetero_dimer_df["direction"]=1
  hetero_dimer_df["rdbe"]=0
  # hetero_dimer_df_duplicate = hetero_dimer_df[duplicated(hetero_dimer_df[,c("node2", "linktype")]) | duplicated(hetero_dimer_df[,c("node2", "linktype")], fromLast = T),]
  
  print("Potentail hetero dimer identified.")
  return(hetero_dimer_df)
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


## Network_prediction - used to connect nodes to library and predict formula ####
Network_prediction = function(Mset, 
                              edge_biotransform, 
                              edge_artifact,
                              biotransform_step = 5,
                              artifact_step = 5,
                              propagation_score_threshold = 0.5,
                              top_n = 50,
                              read_from_csv = F
)
{
  gc()
  if(!read_from_csv){
    mnl=Mset$NodeSet
    edge_biotransform = EdgeSet$Biotransform
    edge_artifact = EdgeSet$Artifacts
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
      for(i in mnl$ID[mnl$category==-1]){
        Initial_formula =sf[[i]]
        Initial_formula[1,]= list(i,mnl$MF[i],0,0,F,1,mnl$rdbe[i])
        sf[[i]]=Initial_formula
      }
    }
    
    #while loop to Predict formula based on known formula and edgelist 
    nrow_experiment = nrow(Mset$Data)
    step=0
    
    timer=Sys.time()
    while(step <= biotransform_step){
      
      all_nodes_df = bind_rows(sf)
      
      # Handle artifacts
      sub_step = 0
      while(sub_step <= 0.01 * artifact_step){
        print(paste("sub_step",sub_step,"elapsed="))
        print((Sys.time()-timer))
        all_nodes_df = bind_rows(sf)
        new_nodes_df = all_nodes_df[all_nodes_df$steps==(step + sub_step)
                                    &all_nodes_df$score>propagation_score_threshold,]
        sub_step = sub_step+0.01
        if(nrow(new_nodes_df)==0){break}
        edge_artifact_sub = edge_artifact[edge_artifact$node1 %in% new_nodes_df$id | 
                                            edge_artifact$node2 %in% new_nodes_df$id,]
        
        for(n in 1:nrow(new_nodes_df)){
          temp_new_node = new_nodes_df[n,]
          flag_id = temp_new_node$id
          flag_formula = temp_new_node$formula
          flag_is_metabolite = temp_new_node$is_metabolite
          flag_score = temp_new_node$score
          flag_rdbe = temp_new_node$rdbe
          
          
          
          # temp_edge_list=subset(edge_artifact_sub, edge_artifact_sub$node1==flag_id |
          #                         edge_artifact_sub$node2==flag_id)
          flag_id_in_node1_or_node2 = edge_artifact_sub$node1==flag_id | edge_artifact_sub$node2==flag_id
          temp_edge_list = edge_artifact_sub[flag_id_in_node1_or_node2,]
          
          
          #If head signal is < defined cutoff, then prevent it from propagating out, but it can still get formula from others.
          if(flag_id <= nrow_experiment){
            if(Mset$Data$mean_inten[flag_id]< 2e4){next}
          }
          
          #If flag is an isotopic peak, then only look for isotopic peaks
          if(grepl("\\[",flag_formula)){
            temp_edge_list = temp_edge_list[grepl("\\[",temp_edge_list$category),]
            if(nrow(temp_edge_list)==0){next}
          }
          
          #Filter edge direction
          temp_edge_list = temp_edge_list[(temp_edge_list$node1 == flag_id & temp_edge_list$direction != -1) 
                                          |(temp_edge_list$node2 == flag_id & temp_edge_list$direction != 1),]
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
            partner_steps = step+sub_step
            partner_is_metabolite = F
            
            temp = sf[[partner_id]]
            
            # temp_subset=subset(temp, temp$formula==partner_formula & temp$is_metabolite==partner_is_metabolite)
            formula_metabolite_status_matched = temp$formula==partner_formula & temp$is_metabolite==partner_is_metabolite
            temp_subset = temp[formula_metabolite_status_matched,]
            
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
      
      # Handle biotransform
      {
        all_nodes_df = bind_rows(sf)
        new_nodes_df = all_nodes_df[all_nodes_df$steps==step
                                    &all_nodes_df$score>propagation_score_threshold,]
        print(paste("nrow",nrow(all_nodes_df),"in step",step,"elapsed="))
        print((Sys.time()-timer))
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
          
          #If head signal is < defined cutoff, then prevent it from propagating out, but it can still get formula from others.
          if(flag_id <= nrow_experiment){
            if(Mset$Data$mean_inten[flag_id]< 2e4){next}
          }
          
          #If flag is an isotopic peak, then only look for isotopic peaks
          if(grepl("\\[",flag_formula)){
            temp_edge_list = temp_edge_list[grepl("\\[",temp_edge_list$category),]
            if(nrow(temp_edge_list)==0){next}
          }
          
          #Filter edge direction
          temp_edge_list = temp_edge_list[(temp_edge_list$node1 == flag_id & temp_edge_list$direction != -1) 
                                          |(temp_edge_list$node2 == flag_id & temp_edge_list$direction != 1),]
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
    gc()
    temp_i=temp_j=1
    triplet_edge_ls_edge=triplet_edge_ls_node=list()
    edge_info = list()
    timer=Sys.time()
    n=1
    for(n in 1:nrow(edge_list)){
      # for(n in 1:20000){
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
      edge_info_sum = rbind(edge_info_sum,edge_info_sum)
      
      
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
  print("finish CPLEXset.")
  return(CPLEX = list(data = CPLEX_data,
                      para = CPLEX_para)
  )
}

### Score_formula ####
Score_formula = function(CPLEXset, rdbe=T, step_score=T, iso_penalty_score=F)
{
  unknown_formula = CPLEXset$data$unknown_formula
  
  #when measured and calculated mass differ, score based on normal distirbution with mean=0 and sd=1e-3
  unknown_formula["msr_mass"] = Mset$Data$medMz[unknown_formula$id]
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
  
  
  unknown_formula["intensity"] = Mset$Data$mean_inten[unknown_formula$id] 
  if(iso_penalty_score==TRUE){
    No_expect_isotope_penalty = function(unknown_formula, linktype){
      isotope_peaks = EdgeSet$Artifacts[grepl("\\[",EdgeSet$Artifacts$linktype),]
      
      unknown_formula["No_isotope_peak"] = !unknown_formula$id %in% isotope_peaks$node1[isotope_peaks$category==linktype]
      if(linktype == "[10]B"){
        unknown_formula["No_isotope_peak"] = !unknown_formula$id %in% isotope_peaks$node2[isotope_peaks$category==linktype]
      }
      unknown_formula["iso_expect_ratio"] = 0
      unknown_formula$iso_expect_ratio[unknown_formula$No_isotope_peak] = unlist(lapply (unknown_formula$formula[unknown_formula$No_isotope_peak], isotopic_abundance,linktype))
      unknown_formula["iso_intensity"] = unknown_formula["intensity"]*unknown_formula["iso_expect_ratio"]
      unknown_formula["iso_penalty_score"] = log10(dnorm(0,unknown_formula$iso_intensity, 2e4) / 
                                                     dnorm(unknown_formula$iso_intensity,unknown_formula$iso_intensity, 2e4)+1e-10)
      return(unknown_formula["iso_penalty_score"])
    }
    unknown_formula["Cl_iso_penalty_score"] = No_expect_isotope_penalty(unknown_formula, linktype="[37]Cl")
    unknown_formula["S_iso_penalty_score"] = No_expect_isotope_penalty(unknown_formula, linktype="[34]S")
    unknown_formula["K_iso_penalty_score"] = No_expect_isotope_penalty(unknown_formula, linktype="[41]K")
    unknown_formula["B_iso_penalty_score"] = No_expect_isotope_penalty(unknown_formula, linktype="[10]B")
    
    unknown_formula["sum_iso_penalty_score"] = base::rowSums(unknown_formula[,grepl("iso_penalty_score", colnames(unknown_formula))])
  } else {
    unknown_formula["sum_iso_penalty_score"] = 0
  }
  
  # the mass score x step score evaluate from mass perspective how likely the formula fits the peak
  # the rdbe score penalizes unsaturation below -1
  #Each node should be non-positive, to avoid node formula without edge connection
  unknown_formula["cplex_score"] = log10(unknown_formula["Mass_score"]) +
    log10(unknown_formula["score"]) +
    unknown_formula["step_score"] + 
    unknown_formula["rdbe_score"] +
    unknown_formula["sum_iso_penalty_score"]
  
  
  # hist(unknown_formula$cplex_score)
  # length(unknown_formula$cplex_score[unknown_formula$cplex_score<1])
  print("Finish scoring formula.")
  return(unknown_formula)
}
### Score_edge_cplex ####
Score_edge_cplex = function(CPLEXset, edge_bonus = -log10(0.5))
{
  edge_info_sum = CPLEXset$data$edge_info_sum
  edge_info_sum = edge_info_sum[with(edge_info_sum, order(edge_ilp_id)),]
  unknown_formula = CPLEXset$data$unknown_formula
  
  {
    edge_info_sum["isotope_score"] = NA
    edge_info_sum["category"] = EdgeSet$Merge$category[edge_info_sum$edge_id]
    edge_info_sum["msr_inten_dif"] =  EdgeSet$Merge$msr_inten_dif[edge_info_sum$edge_id]
    
    edge_info_isotope = edge_info_sum[grepl("\\[", edge_info_sum$category) & !is.na(edge_info_sum$msr_inten_dif),]
    i=8706
    for(i in 1:nrow(edge_info_isotope)){
      temp_edge = EdgeSet$Merge[edge_info_isotope$edge_id[i],]
      temp_iso = temp_edge$category
      if(temp_iso=="[10]B"){
        temp_formula = edge_info_isotope$formula2[i]
        parent_inten = temp_edge$node1_log10_inten
        temp_edge$msr_inten_dif = -temp_edge$msr_inten_dif
      }else {
        temp_formula = edge_info_isotope$formula1[i]
        parent_inten = temp_edge$node2_log10_inten
      }
      calc_abun = isotopic_abundance(temp_formula, temp_iso)
      abun_ratio = 10^temp_edge$msr_inten_dif/calc_abun
      edge_info_isotope$isotope_score[i] = log10(dnorm(abun_ratio,1,0.2+10^(4-parent_inten))/dnorm(1,1,0.2+10^(4-parent_inten))+1e-10)
      # edge_info_sum$isotope_score[edge_info_isotope$edge_ilp_id[i]] = edge_info_isotope$isotope_score[i]
    }
    edge_info_sum$isotope_score[edge_info_isotope$edge_ilp_id] = edge_info_isotope$isotope_score
    # edge_info_sum[edge_info_sum$edge_ilp_id==8716,]
  }
  
  
  edge_info_sum$edge_score = log10(edge_info_sum$edge_score) + edge_bonus
  edge_info_sum$isotope_score = edge_info_sum$isotope_score + edge_bonus*2
  #edge_info_sum$isotope_score = edge_info_sum$isotope_score + edge_bonus
  edge_info_sum$isotope_score[is.na(edge_info_sum$isotope_score)]=0
  edge_info_sum$edge_score = edge_info_sum$edge_score+edge_info_sum$isotope_score
  
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
  
  print("Finish scoring edges.")
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
  
  tictoc::tic()
  
  return_code = mipoptCPLEX(env, prob)
  result_solution=solutionCPLEX(env, prob)
  
  print(paste(return_codeCPLEX(return_code),"-",
              status_codeCPLEX(env, getStatCPLEX(env, prob)),
              " - OBJ_value =", result_solution$objval))
  
  tictoc::toc()
  
  # writeProbCPLEX(env, prob, "prob.lp")
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
  
  return(list(obj = obj, result_solution = result_solution))
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

## CPLEX_screen_edge ####
CPLEX_screen_edge = function(CPLEXset, edge_bonus_range = seq(-.6, -0.9, by=-0.1)){
  
  solution_ls = list()
  #result_solution =1
  for(edge_bonus in edge_bonus_range){
    edge_info_sum = Score_edge_cplex(CPLEXset, edge_bonus = edge_bonus)
    temp_obj = c(CPLEXset$data$unknown_formula$cplex_score,
                 edge_info_sum$edge_score)
    result_solution = Run_CPLEX(CPLEXset, obj = temp_obj)
    solution_ls[[length(solution_ls)+1]] = result_solution
  }
  return(solution_ls)
  
}


## Add_constraint_CPLEX ####
Add_constraint_CPLEX = function(CPLEXset, obj){
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
  
  
  # nnz = sum((unknown_formula$ILP_result!=0)==T)
  # matbeg = 0
  # matval = rep(1,nnz)
  # matind = which(unknown_formula$ILP_result!=0)-1
  # addRowsCPLEX(env, prob, ncols=0, nrows=1, nnz=nnz, matbeg=matbeg, matind=matind, matval=matval,
  #              rhs = base::floor(nnz*.99), sense = "L",
  #              cnames = NULL, rnames = NULL)
  # addRowsCPLEX(env, prob, ncols=0, nrows=1, nnz=1, matbeg=0, matind=3909, matval=1,
  #              rhs = 1, sense = "E",
  #              cnames = NULL, rnames = NULL)
  # delRowsCPLEX(env, prob, begin = nr, end = getNumRowsCPLEX(env, prob)-1)
  # getNumRowsCPLEX(env, prob)
  
  # addMIPstartsCPLEX(env, prob, mcnt = 1, nzcnt = nc, beg = 0, varindices = 1:nc,
  #                   values = CPLEXset$Init_solution2$CPLEX_x, effortlevel = 1, mipstartname = NULL)
  # 
  
  
  tictoc::tic()
  return_code = mipoptCPLEX(env, prob)
  result_solution=solutionCPLEX(env, prob)
  
  print(paste(return_codeCPLEX(return_code),"-",
              status_codeCPLEX(env, getStatCPLEX(env, prob)),
              " - OBJ_value =", result_solution$objval))
  
  tictoc::toc()
  
  # writeProbCPLEX(env, prob, "prob.lp")
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
Subnetwork_analysis = function(g_sub, member_lb = 3, member_ub = 10)
{
  g_sub=g
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
         vertex.label.color = "black",
         vertex.label.cex = 1,
         #edge.color = 'black',
         edge.label = edge.attributes(g_subnetwork_list[[i]][[1]])$linktype,
         vertex.size = 10,
         edge.arrow.size = 2/(length(vertex.attributes(g_subnetwork_list[[i]][[1]])$MF)+1),
         main = paste("Subnetwork",names(subnetwork)[[i]])
    )
  }
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