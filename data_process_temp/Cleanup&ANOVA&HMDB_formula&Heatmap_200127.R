## Library ####
{
  library(readr)
  library(tidyr)
  library(ggplot2)
  #install.packages("stringi")
  library(stringi)
  #install.packages("matrixStats")
  library(matrixStats)
  library(dplyr)
}

## Parameter setup ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  
  hmdb_df = read_csv("hmdb_unique.csv")
  known_mz = read_csv("known_mz.csv")
  hmdb_full_df = read_csv("hmdb_full.csv")
  
  setwd("C:/Users/Li Chen/Desktop/Matt")
  filename = "1.csv"
  # filename = "allM.csv"
  
  raw <- read_csv(filename) %>%
    dplyr::rename(ID = groupId)
  
  mode = -1
  full_hmdb = T
  
  ms_dif_ppm=5
  ms_dif_ppm = ms_dif_ppm/10^6
  
  rt_dif_min=0.2
  
  detection_limit=2500
  replace_random = F
  

}

## Read cohorts ####
{
  first_sample_col_num = 15
  raw = raw %>% 
    filter(complete.cases(.[, first_sample_col_num:ncol(.)]))
  # Data quality control
  # ggplot(raw, aes(x=goodPeakCount, y=maxQuality))+
  #   geom_point()
  # hist(raw$goodPeakCount)
  # hist(raw$maxQuality)
  
  all_names=colnames(raw)[first_sample_col_num:ncol(raw)]
  
  if(length(grep("blank|blk", all_names, ignore.case = T))!=0){
    sample_names=all_names[-grep("blank|blk", all_names, ignore.case = T)]
  } else {
    sample_names=all_names
  }
  blank_names=all_names[grep("blank|blk", all_names, ignore.case = T)]
  sample_cohort=stri_replace_last_regex(sample_names,'[:punct:]?[:alnum:]+', '')
  
  print(sample_names)
  print(sample_cohort)
  print(blank_names)
}


###### Don't need to change from here #####
## Flag similar m/z and RT ####
##Group MS groups
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

##Group RT similar groups based on MS groups
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

##Flag groups for deletion & combine signal
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
    for (n in 14:ncol(raw)){
      s3[k_min,n]=max(s3[k_min:(k_max-1),n])
    }
  }
}

## Flag high blank samples ####
{
  s4 = s3 %>%
    mutate(high_blank = F)
  if(length(blank_names)>0){
    s4["high_blank"]= rowMeans(s4[,sample_names]) < 2*rowMeans(s4[,blank_names])
  }
}

## Annontate base on HMDB mz ####
{
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  
  hmdb_df = hmdb_df %>%
    mutate(adjustmz = Exact_Mass + (H_mass-e_mass)*mode) %>%
    arrange(adjustmz)
  
  s5 = s4 %>%
    mutate(Formula = "",
           Metabolite = "") %>%
    arrange(medMz)
  
  i_min=1
  j=1
  ppm=ms_dif_ppm
  while(j<=nrow(s5)){
    while(hmdb_df$adjustmz[i_min+1] < (s5$medMz[j]*(1-ppm))){
      i_min=i_min+1
    }
    temp_metabolite=as.character()
    temp_formula=as.character()
    k=1
    while(hmdb_df$adjustmz[i_min+k]<(s5$medMz[j]*(1+ppm))){
      if(length(temp_metabolite)==0){
        temp_metabolite=hmdb_df$Name[i_min+k]
        temp_formula=hmdb_df$MF[i_min+k]
      }
      else{
        temp_metabolite=paste(temp_metabolite,";",hmdb_df$Name[i_min+k])
        temp_formula=paste(temp_formula,";",hmdb_df$MF[i_min+k])
      }
      k=k+1
    }
    if(length(temp_metabolite)!=0){
      s5$Metabolite[j]=temp_metabolite
      s5$Formula[j]=temp_formula
    }
    j=j+1
  }
}

## Output files for MetaboAnalyst ####

if(min(table(sample_cohort)) >= 3 & length(table(sample_cohort)) > 1){
  {
    MA_output = s5 %>%
      dplyr::select(ID, sample_names)
    MA_output = rbind(c("Cohort", sample_cohort), MA_output)
    MetaboAnalyst_filename = paste("MA_", filename, sep="")
    write.csv(MA_output, file=MetaboAnalyst_filename, row.names=F)
  }
  
  {
    library(MetaboAnalystR)
    mSet <- InitDataObjects("pktable", "stat", FALSE)
    mSet<-Read.TextData(mSet, MetaboAnalyst_filename, "colu", "disc");
    mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
    mSet<-SanityCheckData(mSet)
    mSet<-ReplaceMin(mSet)
    mSet<-PreparePrenormData(mSet)
    mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
    if(length(table(sample_cohort)) == 2){
      mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, FALSE)
      FDR_file = "t_test.csv"
    } else {
      mSet<-ANOVA.Anal(mSet, F, 0.05, "fisher")
      FDR_file = "anova_posthoc.csv"
    }

    FDR_raw <- read_csv(FDR_file) %>%
      dplyr::select(c("X1","FDR")) %>%
      dplyr::rename(ID = X1) %>%
      mutate(FDR = -log10(FDR))
    
  }
  
  {
    fn <-MetaboAnalyst_filename
    if (file.exists(fn)) file.remove(fn)
    fn <-FDR_file
    if (file.exists(fn)) file.remove(fn)
  }
}

## Output ####
{
  hmdb_match_output=s5 %>%
    arrange(ID) %>%
    filter(flag) %>% 
    filter(!high_blank) %>%
    filter(T) %>%
    dplyr::select(-flag)
  
  if(min(table(sample_cohort)) >= 3 & length(table(sample_cohort)) > 1){
    hmdb_match_output = hmdb_match_output %>%
      merge(FDR_raw, by = "ID") %>%
      mutate(FDR = ifelse(is.na(FDR), 0, FDR))
  } else {
    hmdb_match_output = hmdb_match_output %>%
      mutate(FDR = 0)
  }
  
  refcols = c("ID","medMz","medRt","Formula","Metabolite","FDR","high_blank")
  hmdb_match_output = hmdb_match_output %>%
    dplyr::select(refcols, all_names)
  
  write.csv(hmdb_match_output, file=paste("hmdb_",filename,sep=""), row.names=F)
}

# ## Include alternative names and SMILE ####
# if(full_hmdb){
#   
#   
#   hmdb_SMILE = hmdb_full_df[,c("MF","Name","HMDB_ID","SMILES") ]
#   
#   s6 = s5[!is.na(s5$Formula),]
#   
#   ls = list()
#   for (i in 1:nrow(s6)){
#     temp_df = hmdb_SMILE[hmdb_SMILE$MF==s6$Formula[i],]
#     temp_data = s6[rep(i,nrow(temp_df)),1:(ncol(s6)-2)]
#     ls[[i]] = cbind(temp_df, temp_data)
#   }
#   
#   s7 = bind_rows(ls)
#   
#   
#   write.csv(s7, file=paste("full_hmdb_SMILE_",filename,sep=""), row.names=F)
# }


## Data pre-clean ####
{
  mdata_pre_clean = hmdb_match_output %>%
    filter(Metabolite != "") %>%
    mutate(row_label = paste(stri_sub(Metabolite,1,30), medRt, sep="_")) %>%
    arrange(-FDR) 
  
  row_labels = mdata_pre_clean$row_label
  row_labels = NULL
  
  mdata_pre_clean = mdata_pre_clean %>%
    dplyr::select(sample_names)
  
  mdata_clean = mdata_pre_clean %>% 
    data_impute(impute_method = "threshold", random = F) %>%  # "min", "percentile", "threshold"
    data_normalize_row(nor_method = "") %>% # row_median, row_mean, sample, cohort
    data_normalize_col(nor_method = "col_median") %>% # col_median, col_sum
    data_transform(transform_method = "log2") %>%
    data_scale(scale_method = "mean_center") # mean_center, auto(mean_center + divide sd)
  
  # method "ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"
  # distance "correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  # correlation is default as pearson correlation
  row_cluster = cluster_mat(mdata_clean, method = "complete", distance = "correlation")
  col_cluster = cluster_mat(t(mdata_clean), method = "complete", distance = "correlation")
  
  # col_cluster = F
}


## Heat map ####
{
  my_plot_heatmap(raw_data = mdata_clean,
                  cohort = sample_cohort,
                  row_labels = row_labels, 
                  imgName = sub(".csv","", filename), 
                  format = "pdf", # pdf
                  dpi = 72,
                  width = NA, # define output graph width
                  palette = "RdBu",  # RdBu, gbr, heat, topo
                  viewOpt = "overview", # Detail
                  rowV = row_cluster, # cluster by row, "F" if not clutser
                  colV = col_cluster, # cluster by column, "F" if not clutser
                  border = T, # border for each pixel
                  grp.ave = F, # group average
                  scale_ub = 3, scale_lb = -3 # heatmap scale, auto if ub = lb
  )
}


