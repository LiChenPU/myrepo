# !diagnostics off

# Import library ####

{
  library(lc8)
  library(enviPat)
  library(dplyr)
  library(tidyr)
  # library(fitdistrplus)
  library(slam)
  # library(cplexAPI)
  library(readr)
  library(stringi)
  library(pracma)
  library(igraph)
}

# {
#   library(MetaboAnalystR)
#   library(readr)
#   library(igraph)
#   library(fitdistrplus)
#   library(tidyr)
#   library(dplyr)
#   library(enviPat)
#   library(stringi)
#   library(matrixStats)
#   library(Matrix)
#   library(slam)
#   # library(cplexAPI)
#   
#   library(pracma)
#   library(tictoc)
#   library(janitor)
#   
#   devtools::install_github("LiChenPU/Formula_manipulation")
#   library(lc8)
#   library(profvis)
#   
#   # setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# }

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
  sample_cohort=stri_replace_last_regex(sample_names,'_\\d+|-\\d+|\\.\\d+|\\d+|-a|-b|-c|_mean', '',stri_opts_regex(case_insensitive=T))
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
  
  #raw = raw[complete.cases(raw[, (1+which(colnames(raw)=="parent")):ncol(raw)]),]
  # colnames(raw)[colnames(raw)=="groupId"] = "ID"
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ion_mode = Mset$Global_parameter$mode
  raw = Mset$Raw_data %>%
    rename(ID = groupId) %>%
    mutate(medMz = medMz - (H_mass-e_mass)*ion_mode)
  
  
  ##Group MS groups
  {
    s = raw[with(raw, order(medMz, medRt)),]
    
    mzs = s$medMz
    count = 1
    MZ_group = rep(1,(length(mzs)))
    for(i in 2:length(mzs)){
      if(mzs[i]-mzs[i-1]>mzs[i-1]*ms_dif_ppm){
        count = count+1
      }
      MZ_group[i]=count
    }
    s["MZ_group"]=MZ_group
    }
  
  ##Group RT similar groups based on MS groups
  {
    s2 = s[with(s, order(MZ_group, medRt)),]
    
    rts = s2$medRt
    
    MZRT_group = rep(1,(length(rts)))
    MZ_group = s2$MZ_group
    
    count = 1
    for(i in 2:length(rts)){
      if(MZ_group[i]!=MZ_group[i-1] | rts[i]-rts[i-1]>rt_dif_min){
        count = count+1
      }
      MZRT_group[i]=count
    }
    s2["MZRT_group"] = MZRT_group
    
  }
  
  # Take median of the mz and rt for peaks with same MZRTgroup
  {
    s3 = s2 %>% 
      arrange(MZRT_group)
    ncol_raw = ncol(raw)
    MZRT_group = s3$MZRT_group
    medMz = s3$medMz
    medRt = s3$medRt
    
    k_max=k_min=1
    while (k_max <= length(MZRT_group)){
      k_min = k_max
      while (MZRT_group[k_min] == MZRT_group[k_max]){
        k_max = k_max+1
        if(k_max > length(MZRT_group)){break}
      }
      if(k_max-k_min ==1){next}
      medMz[k_min:(k_max-1)]=median(medMz[k_min:(k_max-1)], na.rm = T)
      medRt[k_min:(k_max-1)]=median(medRt[k_min:(k_max-1)], na.rm = T)
      temp = s3[k_min:(k_max-1),15:ncol_raw]
      temp[1,] = apply(temp, 2, function(x){
        if(any(!is.na(x))){
          return(max(x, na.rm = T))
        } else {
          return(NA)
        }
      })
      s3[k_min:(k_max-1), 15:ncol_raw] = temp[1,]
    }
    
    s3$medMz = medMz
    s3$medRt = medRt
  }
  
  #intermediate files, replace below detection number to random small number
  {
    s4=s3 %>%
      mutate(mean_inten = rowMeans(.[,Mset$Cohort$sample_names], na.rm=T)) %>%
      filter(mean_inten > detection_limit)
    # s4[,4:ncol(s4)][s4[,4:ncol(s4)]<detection_limit]=sample(1:detection_limit, 
    #                                                         size=sum(s4[,4:ncol(s4)]<detection_limit), 
    #                                                         replace=T)
  }
  
  s5 = s4 %>%
    distinct(MZRT_group, .keep_all=T) %>%
    arrange(ID) %>%
    mutate(ID = 1:nrow(.)) %>%
    dplyr::select(-c("MZ_group", "MZRT_group")) %>%
    # mutate(mean_inten = rowMeans(.[,Mset$Cohort$sample_names], na.rm=T)) %>%
    mutate(log10_inten = log10(mean_inten))
  
  
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
  gc()
  if(ncol(mSet$dataSet$norm) > 15000){
    mSet<-my_PlotSubHeatMap(mSet, paste(gsub(".csv","", filename),"top15000_"), "pdf", 600, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 15000, "overview", T, T, T, F)
  } else{
    mSet<-PlotHeatMap(mSet, paste(gsub(".csv","", filename),"full_"), "pdf", 600, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview", T, T, NA, T, F)
  }
  gc()
  mSet<-my_PlotSubHeatMap(mSet, paste(gsub(".csv","", filename),"top50_"), "pdf", 600, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 50, "overview", T, T, T, F)
  gc()
  
  
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
Edge_biotransform = function(Mset, mass_abs = 0.001, mass_ppm = 5)
{
  mass_ppm = mass_ppm/10^6
  merge_node_list = Mset$NodeSet %>%
    filter(category!=-1) %>% #category -1 means adduct node
    arrange(mz)
  
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
      print(paste("Biotransform", Mset$Biotransform$category[k], nrow(temp_edge_list),"found."))
    }
  }
  
  print(Sys.time()-timer)
  edge_list = bind_rows(edge_ls) %>%
    mutate(linktype = Mset$Biotransform$Formula[linktype]) %>%
    filter(node1<=nrow(Mset$Data) | node2<=nrow(Mset$Data), 
           node1!=node2) %>%
    filter(!(linktype == "" & node1<=nrow(Mset$Data) & node2<=nrow(Mset$Data))) %>% # remove data-data isomer connection
    mutate(category = "biotransform")
  
  
  # 
  # edge_list$linktype=Mset$Biotransform$Formula[edge_list$linktype]
  # 
  # edge_list_sub = subset(edge_list, 
  #                        (edge_list$node1<=nrow(Mset$Data)|
  #                           edge_list$node2<=nrow(Mset$Data))&
  #                          edge_list$node1!=edge_list$node2
  # )
  # edge_list_sub["category"]=1
  
  return(edge_list)
}
## Check_sys_measure_error - Check systematic error ####
Check_sys_measure_error = function(Biotransform, inten_threshold=1e5, mass_dif_tol = 2){
  
  Biotransform = EdgeSet$Biotransform
  
  high_inten_node = Mset$Data %>%
    filter(mean_inten > inten_threshold) %>%
    pull(ID)
  
  Biotransform = Biotransform %>%
    filter(linktype == "") %>%
    filter(node1>nrow(Mset$Data)|node2>nrow(Mset$Data)) %>%
    filter(mass_dif<mass_dif_tol) %>%
    filter(node1 %in% high_inten_node | node2 %in% high_inten_node)
  
  Biotransform2 = Biotransform %>%
    mutate(mz_std = ifelse(node1 > node2, Mset$NodeSet$mz[node1], Mset$NodeSet$mz[node2]),
           mz_msr = ifelse(node1 > node2, Mset$NodeSet$mz[node2], Mset$NodeSet$mz[node1])) %>%
    mutate(mz_std_msr_diff = mz_std - mz_msr)
  
  
  ## minimize m_msr*ppm_adjust + abs_adjust + m_msr = m_std
  A_col_1 = Biotransform2$mz_msr
  A_col_2 = rep(1, length(A_col_1))
  b = Biotransform2$mz_std_msr_diff *10^3
  A = cbind(A_col_1, A_col_2)
  lsq_result = lm(Biotransform2$mz_std_msr_diff~Biotransform2$mz_msr)
  ppm_adjust = as.numeric(lsq_result$coefficients[2] * 10^6)
  abs_adjust = as.numeric(lsq_result$coefficients[1])
  
  fitdistData = fitdistrplus::fitdist(b, "norm")
  plot(fitdistData)
  print(fitdistData)
  # shapiro.test(b)
  # 
  # 
  # # selected_ppm_error = Biotransform %>%
  # #   mutate(mass_dif = mass_dif * ifelse(node1>node2, 1, -1)) %>%
  # #   pull(mass_dif)
  # # 
  # # selected_abs_error = Biotransform %>%
  # #   mutate(abs_mass_dif = Mset$NodeSet$mz[node2] - Mset$NodeSet$mz[node1]) %>%
  # #   mutate(abs_mass_dif = abs_mass_dif * ifelse(node1>node2, 1, -1)) %>%
  # #   pull(abs_mass_dif)
  # 
  # # shapiro.test(selected_abs_error)
  # # fitdistData = fitdistrplus::fitdist(selected_ppm_error , "norm")
  # # fitdistData = fitdistrplus::fitdist(selected_abs_error * 1000, "norm")
  # 
  # 
  # 
  # ## fit a abs adjust and a ppm adjust to minimize difference between true values
  # A_col_1 = - Mset$NodeSet$mz[Biotransform$node1]*as.numeric(Biotransform$node1<=nrow(Mset$Data))/Mset$NodeSet$mz[Biotransform$node2]
  # A_col_2 = Mset$NodeSet$mz[Biotransform$node2]*as.numeric(Biotransform$node2<=nrow(Mset$Data))/Mset$NodeSet$mz[Biotransform$node2]
  # A_col_3 = - as.numeric(Biotransform$node1<=nrow(Mset$Data))/Mset$NodeSet$mz[Biotransform$node2]
  # A_col_4 = as.numeric(Biotransform$node2<=nrow(Mset$Data))/Mset$NodeSet$mz[Biotransform$node2]
  # 
  # ## Col1,2 gives ppm_adjust, 3,4 gives abs adjust
  # ## Example: [-mz1/mz2, 0, -1/mz2, 0] * [ppm_adjust, ppm_adjust, abs_adjust, abs_adjust]' = -(mz2-mz1)/mz2
  # A = cbind(A_col_1,A_col_2,A_col_3,A_col_4)
  # b = - Biotransform$mass_dif
  # # minimizes ||A*x - b|| (i.e., in the least-squares sense) subject to C*x = d.
  # # i.e. ppm_adjust or abs_adjust apply to mz1 and mz2 equally,
  # C = matrix(c(1,0,-1,0,0,1,0,-1), nrow = 2)
  # d = matrix(c(0,0), nrow=2)
  # x <- lsqlin(A, b, C,d) 
  
  
  
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
Edge_score = function(Biotransform, mass_dist_sigma,plot_graph = F){
  # Biotransform = EdgeSet$Heterodimer
  # mass_dist_sigma = .5
  library_nodes_ID = Mset$NodeSet %>% 
    filter(category!=1) %>% 
    pull(ID)
  
  Biotransform1 = Biotransform %>%
    filter(node1 %in% library_nodes_ID | node2 %in% library_nodes_ID) %>%
    mutate(edge_massdif_score = dnorm(mass_dif, 0, mass_dist_sigma) + 1e-10) %>%
    mutate(edge_massdif_score = edge_massdif_score/max(edge_massdif_score)) 
  # max(Biotransform1$edge_massdif_score)
  
  Biotransform2 = Biotransform %>%
    filter(!node1 %in% library_nodes_ID, !node2 %in% library_nodes_ID) %>%
    mutate(edge_massdif_score = dnorm(mass_dif, 0, sqrt(2) * mass_dist_sigma)+ 1e-10) %>%
    mutate(edge_massdif_score = edge_massdif_score/max(edge_massdif_score)) 
  
  Biotransform = rbind(Biotransform1, Biotransform2)
  
  
  # if(fix_distribution_sigma){
  #   temp_sigma = ppm_error
  # } else {
  #   if(nrow(Biotransform)>10000){
  #     edge_mzdif_FIT <- fitdist(as.numeric(Biotransform$mass_dif[base::sample(nrow(Biotransform),10000)]), "norm")    
  #   } else {
  #     edge_mzdif_FIT <- fitdist(as.numeric(Biotransform$mass_dif), "norm")    
  #   }
  #   
  #   if(plot_graph){
  #     plot(edge_mzdif_FIT)
  #     print(summary(edge_mzdif_FIT))
  #   }
  #   temp_sigma = edge_mzdif_FIT$estimate[2]
  # }
  # Biotransform["edge_massdif_score"]=dnorm(Biotransform$mass_dif, 0, temp_sigma)
  # Biotransform["edge_massdif_score"]=Biotransform["edge_massdif_score"]/max(Biotransform["edge_massdif_score"])
  return(Biotransform)
}
## Peak_variance - Variance between peaks ####
Peak_variance = function(Mset, 
                         time_cutoff=0.1,
                         TIC_cutoff=10000,
                         correlation_cutoff = -1)
{
  df_raw = Mset$Data %>%
    dplyr::select(c("ID","medMz","medRt",Mset$Cohort$sample_names)) %>%
    mutate(mean_inten = rowMeans(.[,Mset$Cohort$sample_names], na.rm = T)) %>%
    mutate(log10_inten = log10(mean_inten)) %>%
    arrange(medRt)
  
  # df_raw["mean_inten"]=rowMeans(df_raw[,Mset$Cohort$sample_names], na.rm = T)
  # df_raw["log10_inten"]=log10(df_raw$mean_inten)
  
  {
    # df_raw = df_raw[with(df_raw, order(medRt)),]
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
  
  # calculate mass difference between measured peak and library adducts
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
  edge_ls["mz_dif"] = Mset$NodeSet$mz[edge_ls$node2] -Mset$NodeSet$mz[edge_ls$node1]
  # edge_ls["log10_inten_node1"] = Mset$Data$log10_inten[edge_ls$node1]
  # edge_ls["log10_inten_node2"] = Mset$Data$log10_inten[edge_ls$node2]
  # edge_ls["log10_inten_ratio"] = edge_ls["log10_inten_node2"]-edge_ls["log10_inten_node1"]
  # edge_ls["time_dif"] = Mset$Data$medRt[edge_ls$node2] - Mset$Data$medRt[edge_ls$node1]
  
  print("Close RT peaks identified.")
  return(edge_ls)
}

### Artifact_prediction - Edge_list for artifacts ####
Artifact_prediction = function(Mset, Peak_inten_correlation, 
                               search_ms_cutoff=0.002, search_ppm_cutoff = 10)
{
  
  search_ppm_cutoff = search_ppm_cutoff/1e6
  # edge_ls_highcor = EdgeSet$Peak_inten_correlation %>% arrange(mz_dif)
  edge_ls_highcor = Peak_inten_correlation %>% arrange(mz_dif)
  # edge_ls_highcor = edge_ls_highcor[with(edge_ls_highcor, order(mz_dif)),]
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
    
    search_cutoff = max(search_ms_cutoff,search_ppm_cutoff*edge_ls_highcor$mz_node1[i])
    
    
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
                                                    "Oligomer",
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
  
  edge_ls_annotate= temp_df_isotope %>%
    # mutate(category = "Adduct_isotopes") %>%
    rbind(temp_df_oligo) %>%
    dplyr::select(c("node1","node2","linktype","mass_dif","category","direction","rdbe")) 
  
  # edge_ls_annotate_network = edge_ls_annotate[,c("node1","node2","linktype","mass_dif","category","direction","rdbe")]
  
  return(edge_ls_annotate)
}

### Hetero_dimer - Edge_list for hetero_dimer ####
Hetero_dimer = function(Peak_inten_correlation, ppm_tolerance = 5, inten_threshold = 5e5)
{
  # e = EdgeSet$Peak_inten_correlation
  e = Peak_inten_correlation
  e2 = e[e$node1 %in% Mset$Data$ID & e$node2 %in% Mset$Data$ID, ]
  e3 = e2[e2$mz_dif > min(e2$mz_node1) & e2$mz_dif < max(e2$mz_node1),]
  e3_list = split(e3, e3$node2)
  
  ID_inten_threshold = Mset$Data %>%
    filter(mean_inten > inten_threshold) %>%
    pull(ID)
  # e3_list = e3_list[ID_inten_threshold]
  
  hetero_dimer_ls = list()
  for(i in 1: length(e3_list)){
    temp_e = e3_list[[i]]
    temp_matrix = outer(temp_e$mz_node1, temp_e$mz_node1, FUN = "+")  # mz_node1 and mz_dif are the same.
    temp_matrix = (temp_matrix - temp_e$mz_node2[1])/temp_e$mz_node2[1] * 10^6
    temp_index = which(abs(temp_matrix) < ppm_tolerance, arr.ind = T)
    if(length(temp_index)>0){
      temp_ppm = temp_matrix[temp_index]
      temp_node_1 = temp_e$node1[temp_index[,1]]
      linktype = temp_e$node1[temp_index[,2]]
      temp_df = data.frame(node1 = temp_node_1, linktype = linktype, node2 = temp_e$node2[1], mass_dif = temp_ppm)
      hetero_dimer_ls[[length(hetero_dimer_ls)+1]] = temp_df
    }
  }
  hetero_dimer_df = bind_rows(hetero_dimer_ls) %>%
    mutate(category = "Heterodimer",
           direction = 1,
           rdbe = 0) %>%
    filter(node1 != linktype) %>% # remove homo-dimer
    filter(node1 %in% ID_inten_threshold) # retain only high intensity as node1
 
  print("Potential heterodimer identified.")
  return(hetero_dimer_df)
}

### Ring_artifact - Edge_list for Ring_artifact ####
Ring_artifact = function(Peak_inten_correlation, ppm_range_lb = 50, ppm_range_ub = 1000, ring_fold = 50, inten_threshold = 1e6)
{
  e = Peak_inten_correlation %>%
    filter(node1 %in% Mset$Data$ID, node2 %in% Mset$Data$ID) %>%
    filter(mz_dif < max(mz_node1, mz_node2)*ppm_range_ub/1e6) %>%
    mutate(inten_ratio = Mset$Data$mean_inten[node1]/Mset$Data$mean_inten[node2]) %>%
    filter(inten_ratio > ring_fold | inten_ratio < (1/ring_fold))
  
  e2 = e %>%
    filter(inten_ratio < 1) %>%
    mutate(temp_node = node1, node1 = node2, node2 = temp_node) %>%
    mutate(temp_mz = mz_node1, mz_node1 = mz_node2, mz_node2 = temp_mz) %>%
    mutate(mz_dif = -mz_dif, inten_ratio = 1/inten_ratio) %>%
    dplyr::select(-temp_node, -temp_mz)
  
  e3 = e %>%
    filter(inten_ratio > 1) %>%
    rbind(e2)
  
  e3_list = split(e3, e3$node1)
  
  ID_inten_threshold = Mset$Data %>%
    filter(mean_inten > inten_threshold) %>%
    pull(ID)
  # e3_list = e3_list[ID_inten_threshold]
  
  Mass_ring_artifact_ls = list()
  for(i in 1: length(e3_list)){
    temp_e = e3_list[[i]] 
    temp_vector= temp_e$mz_dif / temp_e$mz_node1 * 1e6
    temp_index = which(abs(temp_vector) > ppm_range_lb & abs(temp_vector) < ppm_range_ub) 
    if(length(temp_index)>0){
      temp_ppm = temp_vector[temp_index]
      temp_node_2 = temp_e$node2[temp_index]
      linktype = "Ring_artifact"
      temp_df = data.frame(node1 = temp_e$node1[1], linktype = linktype, node2 = temp_node_2, mass_dif = temp_ppm)
      Mass_ring_artifact_ls[[length(Mass_ring_artifact_ls)+1]] = temp_df
    }
  }
  Mass_ring_artifact_df = bind_rows(Mass_ring_artifact_ls) %>%
    mutate(category = "Ring_artifact",
           direction = 1,
           rdbe = 0, 
           edge_massdif_score = 1) %>%
    filter(node1 %in% ID_inten_threshold)
    
  
  print("Potential Mass_ring_artifact identified.")
  return(Mass_ring_artifact_df)
}

## Merge_edgeset ####
Merge_edgeset = function(EdgeSet, Include_Heterodimer=T, Include_Ring_artifact=T){
  edge_merge = rbind(EdgeSet$Artifacts,EdgeSet$Biotransform)
  if(Include_Heterodimer){
    edge_merge = rbind(edge_merge, EdgeSet$Heterodimer)
  }
  if(Include_Ring_artifact){
    edge_merge = rbind(edge_merge, EdgeSet$Ring_artifact)
  }
  
  edge_merge = edge_merge %>%
    # filter(log10(edge_massdif_score) > filter_log10score) %>%
    mutate(node1_log10_inten = Mset$Data$log10_inten[node1],
           node2_log10_inten = Mset$Data$log10_inten[node2]) %>%
    mutate(edge_id = 1:nrow(.)) %>%
    mutate(msr_inten_dif = node2_log10_inten - node1_log10_inten)
  
  
  return(edge_merge)
}


## Network_prediction - used to connect nodes to library and predict formula ####
Network_prediction = function(Mset, 
                              EdgeSet,
                              biotransform_step = 1,
                              artifact_step = 1,
                              propagation_score_threshold = 0.2,
                              propagation_artifact_intensity_threshold = 2e4,
                              max_formula_num = 1e6,
                              top_n = 50
)
{
  gc()
  # Initialize
  { 
    mnl=Mset$NodeSet
    # edge_artifact = rbind(edge_artifact, edge_heterodimer)
    edge_biotransform = EdgeSet$Merge %>%
      filter(category == "biotransform")
    edge_artifact = EdgeSet$Merge %>%
      filter(category != "biotransform")
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
    }}
 
  
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
      new_nodes_df = all_nodes_df %>%
        filter(steps==(step + sub_step),
               score>propagation_score_threshold) %>%
        filter(!grepl("\\.", formula)) %>%  # Do not propagate from formula with decimal point
        filter(!grepl("Ring_artifact", formula)) # Do not propagate from ring artifacts
      
      sub_step = sub_step+0.01
      if(nrow(new_nodes_df)==0){break}
      
      
      edge_artifact_sub = edge_artifact[edge_artifact$node1 %in% new_nodes_df$id | 
                                          edge_artifact$node2 %in% new_nodes_df$id,]
      
      
      n=9
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
          if(Mset$Data$mean_inten[flag_id]< propagation_artifact_intensity_threshold){next}
        }
        
        #If flag is an isotopic peak, then only look for isotopic peaks or Oligomer/multi-charge peaks or ring artifact
        if(grepl("\\[",flag_formula)){
          temp_edge_list = temp_edge_list %>%
            filter(grepl("\\[|Oligomer|Ring_artifact",category))
          if(nrow(temp_edge_list)==0){next}
        }
        
        #Filter edge direction
        temp_edge_list = temp_edge_list[(temp_edge_list$node1 == flag_id & temp_edge_list$direction != -1) 
                                        |(temp_edge_list$node2 == flag_id & temp_edge_list$direction != 1),]
        if(nrow(temp_edge_list)==0){next}
        
        i=8
        for(i in 1:nrow(temp_edge_list)){
          if(temp_edge_list$category[i] == "Ring_artifact"){
            partner_id = temp_edge_list$node2[i]
            partner_rdbe = flag_rdbe
            partner_formula = paste("Ring_artifact", flag_formula, sep="_")
          } else if(temp_edge_list$category[i] == "Heterodimer"){
            link_edge = sf[[as.numeric(temp_edge_list$linktype[i])]]
            link_edge = link_edge[!grepl("Ring_artifact", link_edge$formula),]  # prevent mass artifact form heterodimer
            if(nrow(link_edge) == 0){next}
            temp_fg = link_edge$formula[1]
            temp_rdbe = link_edge$rdbe[1]
            temp_edge_list$edge_massdif_score[i] = temp_edge_list$edge_massdif_score[i] * link_edge$score[1]
            
            partner_id = temp_edge_list$node2[i]
            partner_rdbe = flag_rdbe + as.numeric(temp_rdbe)
            partner_formula=my_calculate_formula(flag_formula,temp_fg,1,Is_valid = T)
          } else {
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
                partner_formula=my_calculate_formula(flag_formula,flag_formula,-(fold-1)/fold,Is_valid = F)
                # Remove decimal point in formula
                # if(grepl("\\.", partner_formula)){partner_formula=F}
              }else {
                partner_formula=my_calculate_formula(flag_formula,temp_fg,-1,Is_valid = T)
              }
            }
            
            #function return false if not valid formula
            if(is.logical(partner_formula)){next}
            
          }
          
          
          
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
      new_nodes_df = all_nodes_df %>%
        filter(steps==step,
               score>propagation_score_threshold, 
               is_metabolite)
      print(paste("nrow",nrow(all_nodes_df),"in step",step,"elapsed="))
      print((Sys.time()-timer))
      step = step + 1
      if(nrow(new_nodes_df)==0 | step > biotransform_step | nrow(new_nodes_df) > max_formula_num){break}
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
    pred_formula = pred_formula[!duplicated(pred_formula[,c("id","formula")]),]
    
    merge_formula = pred_formula
    sf = list()
    for(n in 1: max(merge_formula$id)){
      sf[[n]]=merge_formula[merge_formula$id==n,]
    }
    
  }
  
  return(sf)
}


# Function for CPLEX ####
### Prepare_CPLEX_formula ####
Prepare_CPLEX_formula = function(Mset, mass_dist_sigma, rdbe=F, step_score=F, iso_penalty_score=F, filter_low_score = -9)
{
  raw_pred_formula=bind_rows(Mset[["NodeSet_network"]])
  raw_node_list = Mset$NodeSet
  
  #Clean up
  {
    pred_formula = raw_pred_formula
    
    #1 is measured peaks; 0 is library; -1 is system adduct
    lib_nodes = raw_node_list[raw_node_list$category!= 1,]
    lib_nodes_cutoff = nrow(Mset$Data)
    unknown_nodes = raw_node_list[raw_node_list$category==1,]
    unknown_nodes = unknown_nodes[unknown_nodes$ID %in% unique(pred_formula$id),]
    num_unknown_nodes = nrow(unknown_nodes)
    
    lib_formula = pred_formula[pred_formula$id %in% lib_nodes$ID,]
    lib_formula = lib_formula[lib_formula$steps==0,]
    unknown_formula = pred_formula[pred_formula$id %in% unknown_nodes$ID,]
    
  }
  
  unknown_formula = unknown_formula %>%
    mutate(msr_mass = Mset$NodeSet$mz[id]) %>%
    mutate(ILP_id = 1:nrow(.))
  
  #when measured and calculated mass differ, score based on normal distirbution with mean=0 and sd=1e-3
  unknown_formula_ring_artifact = unknown_formula %>%
    filter(grepl("Ring_artifact", formula)) %>%
    mutate(cal_mass = msr_mass,
           msr_cal_mass_dif = 0,
           msr_cal_mass_dif_ppm = 0,
           Mass_score = 1)
  unknown_formula = unknown_formula %>%
    filter(!grepl("Ring_artifact", formula)) %>%
    mutate(cal_mass = formula_mz(formula)) %>%
    mutate(msr_cal_mass_dif = msr_mass - cal_mass) %>%
    mutate(msr_cal_mass_dif_ppm = msr_cal_mass_dif / msr_mass * 1e6) %>%
    mutate(Mass_score = dnorm(msr_cal_mass_dif_ppm, 0, mass_dist_sigma)/dnorm(0, 0, mass_dist_sigma)) %>%
    mutate(Mass_score = Mass_score+1e-10) %>%
    rbind(unknown_formula_ring_artifact)
  
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
    # log10(unknown_formula["score"]) +
    unknown_formula["step_score"] + 
    unknown_formula["rdbe_score"] +
    unknown_formula["sum_iso_penalty_score"]
  
  unknown_formula = unknown_formula %>%
    arrange(ILP_id) %>%
    # filter(cplex_score > filter_low_score) %>% arrange(ILP_id) %>% # reassigned ILP_id for each formula to avoid gapping in id
    mutate(cplex_score = round(cplex_score, digits = 6))
    
  merge_formula = unknown_formula %>%
    select(colnames(lib_formula)) %>%
    rbind(lib_formula) %>%
    mutate(ILP_id = 1:nrow(.))
  
  #Select edge list where it relates only to unknown nodes that have propagated formula
  
  pred_formula_ls = list()
  merge_formula_id=merge_formula$id
  for(n in 1: max(merge_formula$id)){
    pred_formula_ls[[n]]=merge_formula[merge_formula_id==n,]
  }
  
  # hist(unknown_formula$cplex_score)
  # length(unknown_formula$cplex_score[unknown_formula$cplex_score<1])
  print("Finish scoring formula.")
  
  
  CPLEX_formula = list(pred_formula_ls = pred_formula_ls,
                    unknown_nodes = unknown_nodes,
                    unknown_formula = unknown_formula
  )
  return(CPLEX_formula)
}

### Prepare_CPLEX_edge ####
Prepare_CPLEX_edge = function(EdgeSet, CPLEXset, edge_bonus, isotope_bonus, artifact_bonus)
{
  raw_edge_list = EdgeSet$Merge
  unknown_nodes = CPLEXset$formula$unknown_nodes
  unknown_formula = CPLEXset$formula$unknown_formula
  pred_formula_ls = CPLEXset$formula$pred_formula_ls
  
  edge_list = raw_edge_list[raw_edge_list$node1 %in% unknown_nodes$ID |
                              raw_edge_list$node2 %in% unknown_nodes$ID,]
  
  
  edge_info = list()
  counts_edge_info = 1
  timer=Sys.time()
  
  # edge_list = edge_list %>%
  #   sample_n(1000)
  for(n in 1:nrow(edge_list)){
    if(n%%10000==0){
      print(paste("n=",n,"elapsed="))
      print(Sys.time()-timer)
    }
    
    temp_edge = edge_list[n,]
    
    node_1 = temp_edge$node1
    node_2 = temp_edge$node2
    formula_1 = pred_formula_ls[[node_1]]
    formula_2 = pred_formula_ls[[node_2]]

    if(temp_edge$category[1] == "Heterodimer"){
      link_fg = pred_formula_ls[[as.numeric(temp_edge$linktype)]]$formula
      link_fg = link_fg[!grepl("Ring_artifact", link_fg)]  # prevent mass artifact form heterodimer
    } else {
      link_fg = temp_edge$linktype
    }
    
    temp_score = temp_edge$edge_massdif_score
    for(temp_formula in unique(formula_1$formula)){
      for(temp_fg in unique(link_fg)){
        #Assuming formula in node_1 is always smaller than node_2
        if(temp_fg==""){
          temp_formula_2 = temp_formula
        }else if (grepl("Ring_artifact",temp_fg)){
          temp_formula_2 = paste("Ring_artifact", temp_formula, sep="_")
        }else if (grepl("x",temp_fg)){
          fold = as.numeric(gsub("x","",temp_fg))
          temp_formula_2=my_calculate_formula(temp_formula,temp_formula,fold-1,Is_valid = F)
        }else {
          temp_formula_2=my_calculate_formula(temp_formula,temp_fg,1,Is_valid = T)
        }
        
        #Write edge_info for edge and corresponding 2 nodes
        if(any(temp_formula_2 == formula_2$formula)){
          temp_j1 = formula_1$ILP_id[which(formula_1$formula==temp_formula )]
          temp_j2 = formula_2$ILP_id[which(formula_2$formula==temp_formula_2 )]
          edge_info[[counts_edge_info]] = list(edge_id=temp_edge$edge_id,
                                     edge_score=temp_score,
                                     formula1 = temp_formula,
                                     formula2 = temp_formula_2,
                                     ILP_id1 = temp_j1,
                                     ILP_id2 = temp_j2)
          counts_edge_info = counts_edge_info + 1

        }
      } 
    }
  }
  
  ##Objective parameter 
  {
    edge_info_sum = bind_rows(edge_info) %>%
      mutate(edge_ilp_id = 1:nrow(.))
    
    edge_info_sum_with_library1 = edge_info_sum %>%
      filter(ILP_id1>nrow(unknown_formula))
      
    edge_info_sum_with_library2 = edge_info_sum %>%
      filter(ILP_id2>nrow(unknown_formula))
    colnames(edge_info_sum_with_library2)=sub(1,3,colnames(edge_info_sum_with_library2))
    colnames(edge_info_sum_with_library2)=sub(2,1,colnames(edge_info_sum_with_library2))
    colnames(edge_info_sum_with_library2)=sub(3,2,colnames(edge_info_sum_with_library2))
      
    library_keep_id = rbind(edge_info_sum_with_library1, edge_info_sum_with_library2) %>%
      arrange(-edge_score) %>%
      distinct(formula1, ILP_id1, .keep_all = T) %>%
      pull(edge_ilp_id)
    
    library_remove_id = rbind(edge_info_sum_with_library1, edge_info_sum_with_library2) %>%
      filter(!edge_ilp_id %in% library_keep_id) %>%
      pull(edge_ilp_id)
  
    edge_info_sum = edge_info_sum %>%
      mutate(edge_score = ifelse(edge_ilp_id %in% library_remove_id, 1e-10, edge_score))
  }
  
  
  
  edge_info_sum = edge_info_sum %>%
    arrange(edge_ilp_id) %>%
    mutate(node1 = EdgeSet$Merge$node1[edge_id],
           node2 = EdgeSet$Merge$node2[edge_id],
           direction = EdgeSet$Merge$direction[edge_id],
           linktype = EdgeSet$Merge$linktype[edge_id],
           category = EdgeSet$Merge$category[edge_id]) %>%
    mutate(mz1 = Mset$NodeSet$mz[node1], mz2 = Mset$NodeSet$mz[node2])
  
  # Give artifact bonus to all artifact connections
  {
    edge_info_sum = edge_info_sum %>%
      mutate(artifact_score = ifelse(category == "biotransform", 0, artifact_bonus))
    
  }
  
  
  # Calculate isotope scores 
  {
    edge_info_sum = edge_info_sum %>%
      mutate(isotope_score = NA) %>%
      mutate(category = EdgeSet$Merge$category[edge_id]) %>%
      mutate(msr_inten_dif = EdgeSet$Merge$msr_inten_dif[edge_id])
    # edge_info_sum["isotope_score"] = NA
    # edge_info_sum["category"] = EdgeSet$Merge$category[edge_info_sum$edge_id]
    # edge_info_sum["msr_inten_dif"] =  EdgeSet$Merge$msr_inten_dif[edge_info_sum$edge_id]
    
    edge_info_isotope = edge_info_sum %>%
      filter(grepl("\\[", category), !is.na(msr_inten_dif))
    
    i=1
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
  
  # combine isotope scores and edge_difference scores
  
  edge_info_sum = edge_info_sum %>%
    mutate(edge_score = log10(edge_score) + edge_bonus) %>%
    mutate(isotope_score = isotope_score + isotope_bonus) %>%
    mutate(isotope_score = replace_na(isotope_score, 0)) %>%
    mutate(edge_score = edge_score + isotope_score + artifact_score)

  {
    edge_info_sum2 = edge_info_sum %>%
      filter(ILP_id2<=nrow(unknown_formula), 
             ILP_id1<=nrow(unknown_formula))
    
    edge_info_same12 = edge_info_sum2 %>% filter(formula1 == formula2)
    df_same12 = table(edge_info_same12$formula1)
    
    edge_info_dif12 = edge_info_sum2 %>% filter(formula1 != formula2)
    edge_info_dif12 = edge_info_dif12[duplicated(edge_info_dif12[,c("formula1","formula2")]) | 
                                        duplicated(edge_info_dif12[,c("formula1","formula2")], fromLast=TRUE),]
    temp_merge = with(edge_info_dif12, paste0(formula1, formula2))
    df_dif12 = table(temp_merge)
    
    # Scenerio: Each duplicated mz will generate duplicated edge connectoin, exaggerating the score
    # Solution: when n duplicated mz exist, the intra mz edge number are n*(n+1)/2, but effectively they should only get n/2
    
    sol_mat = data.frame(n=1:max(df_same12, df_dif12), div = NA)
    sol_mat$div = (-1+sqrt(1+8*sol_mat$n))/2
    
    edge_info_same12 = edge_info_same12 %>%
      mutate(edge_score = edge_score / sol_mat$div[df_same12[formula1]])
    edge_info_dif12 = edge_info_dif12 %>%
      mutate(edge_score = edge_score / sol_mat$div[df_dif12[temp_merge]])
    
  }
  
  
  temp_edge_info_sum = rbind(edge_info_same12,edge_info_dif12,edge_info_sum) %>%
    distinct(edge_ilp_id, .keep_all = T) %>%
    arrange(edge_ilp_id) 
  
  print("Finish scoring edges.")
  return(temp_edge_info_sum)
  
}
## Prepare_CPLEX parameter ####
Prepare_CPLEX_para = function(Mset, EdgeSet, CPLEXset){
  unknown_nodes = CPLEXset$formula$unknown_nodes
  unknown_formula = CPLEXset$formula$unknown_formula
  pred_formula_ls = CPLEXset$formula$pred_formula_ls
  lib_nodes_cutoff = nrow(Mset$Data)
  
  edge_info_sum = CPLEXset$edge
  edge_list = edge_info_sum %>%
    filter(edge_score>0)
    
  
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
  count_isotope = 1
  triplet_edge_ls_edge = triplet_edge_ls_node = list()
  triplet_isotope_ls_edge = triplet_isotope_ls_node = list()
  edge_info = list()
  timer=Sys.time()
  
  n=37
  # edge_list = edge_list %>%
  #   sample_n(1000)
  for(n in 1:nrow(edge_list)){
    # for(n in 1:20000){
    if(n%%10000==0){
      # print(paste("n=",n,"elapsed="))
      # print(Sys.time()-timer)
    }
    
    temp_edge = edge_list[n,]
    # Select edges that are included in CPLEX optimization
    # if(temp_edge$edge_score<0){next}
    # if(temp_edge$category == "Heterodimer"){next}
    # if(temp_edge$category == "Ring_artifact"){next}
    
    node_1 = temp_edge$node1
    node_2 = temp_edge$node2
    
    #Write triplet for edge and corresponding 2 nodes
  
    temp_j1 = temp_edge$ILP_id1
    temp_j2 = temp_edge$ILP_id2
    
    #if one node is library node,
    if(node_1>lib_nodes_cutoff){
      triplet_edge_ls_edge[[temp_i]] = list(i=temp_i,
                                            j=temp_j,
                                            v=1)
      triplet_edge_ls_node[[temp_i]] = list(i=temp_i,
                                            j=temp_j2,
                                            v=-1)
    }
    if(node_2>lib_nodes_cutoff){
      triplet_edge_ls_edge[[temp_i]] = list(i=temp_i,
                                            j=temp_j,
                                            v=1)
      triplet_edge_ls_node[[temp_i]] = list(i=temp_i,
                                            j=temp_j1,
                                            v=-1)
    }
    
    #if both nodes are unknown nodes
    if(node_1<=lib_nodes_cutoff&node_2<=lib_nodes_cutoff){
      triplet_edge_ls_edge[[temp_i]] = list(i=temp_i,
                                            j=temp_j,
                                            v=2)
      triplet_edge_ls_node[[temp_i]] = list(i=c(temp_i,temp_i),
                                            j=c(temp_j1,temp_j2),
                                            v=c(-1,-1))
      

      # Force a isotope formula comes with an edge connection
      if(grepl("\\[", temp_edge$category)){
        if(temp_edge$category != "[10]B"){
          triplet_isotope_ls_edge[[count_isotope]] = list(i=count_isotope,
                                                          j=temp_j,
                                                          v=1)
          triplet_isotope_ls_node[[count_isotope]] = list(i=count_isotope,
                                                          j=temp_j2,
                                                          v=-1)
        } else{
          triplet_isotope_ls_edge[[count_isotope]] = list(i=count_isotope,
                                                          j=temp_j,
                                                          v=1)
          triplet_isotope_ls_node[[count_isotope]] = list(i=count_isotope,
                                                          j=temp_j1,
                                                          v=-1)
          
        }

        count_isotope = count_isotope + 1
      }
      
    }
    
    temp_i = temp_i+1
    temp_j = temp_j+1
  }
  

  ## Because triplet_unknown_node take nrow(unknown_nodes) rows, and nrow(unknown_formula) columns
  triplet_edge_ls_edge_sum = bind_rows(triplet_edge_ls_edge) %>%
    mutate(i = i + nrow(unknown_nodes),
           j = j + nrow(unknown_formula))
  triplet_edge_ls_node_sum = bind_rows(triplet_edge_ls_node) %>%
    mutate(i = i + nrow(unknown_nodes))
  
  triplet_isotope_ls_edge_sum = bind_rows(triplet_isotope_ls_edge) %>%
    mutate(i = i + max(triplet_edge_ls_node_sum$i),
           j = j + nrow(unknown_formula))
  triplet_isotope_ls_node_sum = bind_rows(triplet_isotope_ls_node) %>%
    mutate(i = i + max(triplet_edge_ls_node_sum$i))

  
  #Generate sparse matrix on left hand side
  triplet_df = rbind(
    triplet_unknown_node, 
    triplet_edge_ls_edge_sum,
    triplet_edge_ls_node_sum,
    triplet_isotope_ls_edge_sum,
    triplet_isotope_ls_node_sum
  )
  
  mat = simple_triplet_matrix(i=triplet_df$i,
                              j=triplet_df$j,
                              v=triplet_df$v)
  
  
  #CPLEX solver parameter
  {
    nc <- max(mat$j)
    obj <- c(CPLEXset$formula$unknown_formula$cplex_score, 
                           CPLEXset$edge$edge_score[CPLEXset$edge$edge_score>=0])
    lb <- rep(0, nc)
    ub <- rep(1, nc)
    ctype <- rep("B",nc)
    
    nr <- max(mat$i)
    rhs = c(rep(1,nrow(unknown_nodes)),rep(0,nrow(edge_list)),rep(0, length(triplet_isotope_ls_edge)))
    sense <- c(rep("L",nrow(unknown_nodes)), rep("L", nrow(edge_list)), rep("E", length(triplet_isotope_ls_edge)))
    
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
  print("finish CPLEXset parameter.")
  return(CPLEX_para)
  
}

## Run_CPLEX ####
Run_CPLEX = function(CPLEXset, obj_cplex){
  
  # obj_cplex = CPLEXset$para$obj
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
  
  
  
  copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj = obj_cplex, rhs, sense,
                    beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
  
  
  copyColTypeCPLEX(env, prob, ctype)
  
  # Conserve memory true
  setIntParmCPLEX(env, CPX_PARAM_MEMORYEMPHASIS, CPX_ON)
  setIntParmCPLEX(env, CPX_PARAM_PROBE, 3)
  # setIntParmCPLEX(env, CPX_PARAM_INTSOLLIM, 2)
  # setDefaultParmCPLEX(env)
  # getChgParmCPLEX(env)
  
  # Assess parameters
  # getParmNameCPLEX(env, 1082)
  
  
  # Access Relative Objective Gap for a MIP Optimization Description
  # getMIPrelGapCPLEX(env, prob)
  
  tictoc::tic()
  # test = basicPresolveCPLEX(env, prob)
  return_code = mipoptCPLEX(env, prob)
  result_solution=solutionCPLEX(env, prob)
  # result_solution_info = solnInfoCPLEX(env, prob)
  
  print(paste(return_codeCPLEX(return_code),"-",
              status_codeCPLEX(env, getStatCPLEX(env, prob)),
              " - OBJ_value =", result_solution$objval))
  tictoc::toc()
  
  # writeProbCPLEX(env, prob, "prob.lp")
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
  
  return(list(obj = obj_cplex, result_solution = result_solution))
}
## Test_para_CPLEX ####
Test_para_CPLEX = function(CPLEXset, obj_cplex,  test_para = -1:4){
  
  # for(test_para_CPX_PARAM_PROBE in test_para1){
    for(temp_para in test_para){

    # print(test_para_CPX_PARAM_PROBE)
      print(temp_para)
    
    # obj_cplex = CPLEXset$para$obj
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
    
    
    
    copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj = obj_cplex, rhs, sense,
                      beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
    
    
    copyColTypeCPLEX(env, prob, ctype)
    
    # Conserve memory true
    setIntParmCPLEX(env, CPX_PARAM_MEMORYEMPHASIS, CPX_ON)
    setIntParmCPLEX(env, CPX_PARAM_PROBE, 3)
    
    # Set time is Dbl not Int
    # setDblParmCPLEX(env, CPX_PARAM_TILIM, 1000) # total run time
    # setIntParmCPLEX(env, CPX_PARAM_TUNINGTILIM, 200) # run time for each tuning (each optimizatoin run will test serveral tuning)
    
    
    
    # setIntParmCPLEX(env, CPX_PARAM_BBINTERVAL, temp_para)
    # setIntParmCPLEX(env, CPX_PARAM_NODESEL, temp_para) # 0:3 No effect
    setIntParmCPLEX(env, CPX_PARAM_CLIQUES, temp_para) # 0:3 No effect
    # 
    
    # setIntParmCPLEX(env, CPX_PARAM_CLIQUES, temp_para)
    
    # setIntParmCPLEX(env, CPX_PARAM_INTSOLLIM, 2)
    # setIntParmCPLEX(env, CPX_PARAM_PROBE, 2)
    # setDefaultParmCPLEX(env)
    # getChgParmCPLEX(env)
    
    # Assess parameters
    # getParmNameCPLEX(env, 1082)
    
    
    # Access Relative Objective Gap for a MIP Optimization Description
    # getMIPrelGapCPLEX(env, prob)
    
    tictoc::tic()
    # test = basicPresolveCPLEX(env, prob)
    return_code = mipoptCPLEX(env, prob)
    result_solution=solutionCPLEX(env, prob)
    # result_solution_info = solnInfoCPLEX(env, prob)
    
    print(paste(return_codeCPLEX(return_code),"-",
                status_codeCPLEX(env, getStatCPLEX(env, prob)),
                " - OBJ_value =", result_solution$objval))
    tictoc::toc()
    
    # writeProbCPLEX(env, prob, "prob.lp")
    delProbCPLEX(env, prob)
    closeEnvCPLEX(env)
  # }
  }
  
  return(0)
}
## CPLEX_permutation ####
CPLEX_permutation = function(CPLEXset, n_pmt = 5, sd_rel_max = 0.5){
  unknown_formula = CPLEXset$formula$unknown_formula
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
## determine_is_metabolite ####
determine_is_metabolite = function(){
  ## Prepare formula list 
  {
    # formula_list = merge(Mset$NodeSet, unknown_formula,by.x = "ID", by.y = "id",all=T)
    
    pred_formula_df = bind_rows(CPLEXset$formula$pred_formula_ls)
    nodeset_df = Mset$NodeSet %>% dplyr::select(-MF, -rdbe)
    
    formula_list = merge(unknown_formula, pred_formula_df, all = T) %>%
      merge(nodeset_df, by.x = "id", by.y = "ID", all = T) %>%
      select(ILP_id, mz, RT, everything()) %>%
      arrange(ILP_id)
# 
#     formula_list$formula[formula_list$ID>nrow(Mset$Data)] = formula_list$MF[formula_list$ID>nrow(Mset$Data)]
#     formula_list$rdbe.y[formula_list$ID>nrow(Mset$Data)] = formula_list$rdbe.x[formula_list$ID>nrow(Mset$Data)]
#     formula_list=formula_list[,!(colnames(formula_list) %in% c("MF", "rdbe.x", "is_metabolite"))]

    formula_list = formula_list %>%
      select(ILP_id, everything())
    
    ilp_id_non0 = formula_list %>%
      filter(ILP_result!=0) %>%
      pull(ILP_id)
  }
  ## prepare relation list
  {
    # Select all connections with both sides are confirmed formula
    relation_list = edge_info_sum %>%
      filter((ILP_id1 %in% ilp_id_non0 | node1 > nrow(Mset$Data)),
             (ILP_id2 %in% ilp_id_non0 | node2 > nrow(Mset$Data))) %>%
      select(ILP_id1, ILP_id2, everything())
  }
  ## Define if formula is artifact
  {
    artifact_edgeset = relation_list %>%
      filter(category != "biotransform") 
    
    Heterodimer_ILP_id = artifact_edgeset %>% 
      filter(category == "Heterodimer") %>%
      pull(ILP_id2)
    
    formula_list_step = formula_list %>%
      dplyr::select(ILP_id, steps) %>%
      drop_na()
    # formula_list_step$steps[7]
    # formula_list_step$steps[15856]
    
    Oligomer_ILP_id = artifact_edgeset %>% 
      filter(category == "Oligomer") %>%
      filter(formula_list_step$step[ILP_id1] %% 1 <= formula_list_step$step[ILP_id2] %% 1) %>%
      pull(ILP_id2)
    
    Double_charge_ILP_id = artifact_edgeset %>% 
      filter(category == "Oligomer") %>%
      filter(formula_list$steps[ILP_id1] %% 1 > formula_list$steps[ILP_id2] %% 1) %>%
      filter(grepl("\\.", formula_list$formula[ILP_id1])) %>%
      pull(ILP_id1)
    
    Isotope_ILP_id = artifact_edgeset %>% 
      filter(grepl("\\[", category)) %>%
      filter(direction == 1) %>%
      pull(ILP_id2) 
    
    Isotope_ILP_id = artifact_edgeset %>% 
      filter(grepl("\\[", category)) %>%
      filter(direction != 1) %>%
      pull(ILP_id1) %>%
      c(Isotope_ILP_id)

    Ring_artifact_ILP_id = artifact_edgeset %>% 
      filter(category == "Ring_artifact") %>%
      pull(ILP_id2)
    
    Fragment_ILP_id = artifact_edgeset %>% 
      filter(direction != 1, 
             !grepl("\\[", category), 
             category != "Oligomer") %>%
      pull(ILP_id1)
    
    Adduct_ILP_id_fragment_link = artifact_edgeset %>%
      filter(direction != 1, 
             !grepl("\\[", category), 
             category != "Oligomer") %>%
      pull(ILP_id2)
    
    Adduct_ILP_id = artifact_edgeset %>%
      filter(category != "Heterodimer",
             category != "Oligomer",
             !grepl("\\[", category),
             category != "Ring_artifact",
             direction == 1) %>%
      pull(ILP_id2) %>%
      c(Adduct_ILP_id_fragment_link)
    
    
    
    
    # length((c(Fragment_ILP_id, Ring_artifact_ILP_id, Isotope_ILP_id, Oligomer_ILP_id, Heterodimer_ILP_id)))
    # length(Adduct_ILP_id)
    
    Artifact_assignment = rep("", max(formula_list$ILP_id, na.rm = T))
    
    Artifact_assignment[Adduct_ILP_id] = paste0(Artifact_assignment[Adduct_ILP_id],"|","Adduct")
    Artifact_assignment[Fragment_ILP_id] = paste0(Artifact_assignment[Fragment_ILP_id],"|","Fragment")
    Artifact_assignment[Heterodimer_ILP_id] = paste0(Artifact_assignment[Heterodimer_ILP_id],"|","Heterodimer")
    Artifact_assignment[Oligomer_ILP_id] = paste0(Artifact_assignment[Oligomer_ILP_id],"|","Oligomer")
    Artifact_assignment[Isotope_ILP_id] = paste0(Artifact_assignment[Isotope_ILP_id],"|","Isotope")
    Artifact_assignment[Ring_artifact_ILP_id] = paste0(Artifact_assignment[Ring_artifact_ILP_id],"|","Ring_artifact")
    Artifact_assignment[Double_charge_ILP_id] = paste0(Artifact_assignment[Double_charge_ILP_id],"|","Double_charge")
    
    
    formula_list["Artifact_assignment"] = Artifact_assignment[formula_list$ILP_id]
  }
  ## Define if formula is metabolite
  {
    biotransform_edgeset = relation_list %>%
      filter(category == "biotransform") 
    
    g_bio = graph_from_data_frame(d = biotransform_edgeset, 
                                  vertices =  formula_list %>%
                                    filter(ILP_id %in% c(biotransform_edgeset$ILP_id1, biotransform_edgeset$ILP_id2)),
                                  directed = F)
    clu=components(g_bio)
    # Bio_subnetwork uses all biotransformation to make network
    # a ILP_id is defined metabolite, if the network contains library metabolites
    library_ILP_id = formula_list %>% filter(category == 0) %>% pull(ILP_id)
    g_bio_subnetwork = igraph::groups(clu)
    g_bio_subnetwork_membership = as.data.frame(clu$membership) %>%
      filter(as.numeric(row.names(.)) %in% library_ILP_id) %>%
      pull(1) %>%
      unique()
    
    biotranform_nodes_id = c()
    for(i in g_bio_subnetwork_membership){
      biotranform_nodes_id = c(biotranform_nodes_id, as.numeric(g_bio_subnetwork[[i]]))
    }
    
    #biotranform_nodes = unique(c(biotranform_edgeset$node1, biotranform_edgeset$node2))
    formula_list = formula_list %>%
      mutate(Biotransform = ifelse(is.na(formula_list$formula), NA, ILP_id %in% biotranform_nodes_id)) %>%
      mutate(is_metabolite = case_when(
        Artifact_assignment == "" & Biotransform ~ "Yes",
        Artifact_assignment != "" & Biotransform ~ "Maybe",
        Artifact_assignment != "" & !Biotransform ~ "No",
        Artifact_assignment == "" & !Biotransform ~ "NA"
      ))
  }
  return(list(formula_list = formula_list,
              relation_list = relation_list))
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
