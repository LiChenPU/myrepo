## Library ####
{
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringi)
  library(matrixStats)
  library(dplyr)
  library(rstudioapi)
  
  # install.packages("stringi")
}

## Parameter setup ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  hmdb_df = read_csv("hmdb_unique.csv")
  # put the peak table file directory here
  setwd("./sample files") 
  setwd("C:/Users/Li Chen/Desktop/Github local/myrepo/NetID_1204/matt")

  # Filename
  filename = "raw_data.csv"

  # specify sample cohorts, if auto, then the code will try to find patterns
  # Otherwise, specify sample cohorts
  sample_cohort = "auto"
  # sample_cohort = c(rep("MS_pos_mp",30),rep("MS_neg", 32), rep("MS_neg_mp", 32))
  
  # Other parameters
  mode = 1 # ionization mode
  ms_dif_ppm=5 # Define mass error tolerance for HMDB matching
  rt_dif_min=0.2 # When two peaks have same mass and  RT < 0.2 min, the smaller one will be deleted
  High_blank_ratio = 2 # filter peaks that are 2x less than blanks

}

###### Don't need to change from here #####
## Read cohorts ####
{
  
  raw <- read_csv(filename) %>%
    dplyr::rename(ID = groupId)
  ms_dif_ppm = ms_dif_ppm/10^6
  
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
  if(sample_cohort[[1]] == "auto"){
    sample_cohort=stri_replace_last_regex(sample_names,'[:punct:]?[:alnum:]+', '')
  }
  
  print(data.frame(sample_names, sample_cohort))
  print(blank_names)
}

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
    mutate(High_blank = F)
  if(length(blank_names)>0){
    s4["High_blank"]= rowMeans(s4[,sample_names]) < 2*rowMeans(s4[,blank_names])
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
           Metabolite = "",
           HMDB_ID = "") %>%
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
        temp_ID = hmdb_df$HMDB_ID[i_min+k]
      }
      else{
        temp_metabolite=paste(temp_metabolite,";",hmdb_df$Name[i_min+k])
        temp_formula=paste(temp_formula,";",hmdb_df$MF[i_min+k])
        temp_ID=paste(temp_ID,";",hmdb_df$HMDB_ID[i_min+k])
      }
      k=k+1
    }
    if(length(temp_metabolite)!=0){
      s5$Metabolite[j]=temp_metabolite
      s5$Formula[j]=temp_formula
      s5$HMDB_ID[j]=temp_ID
    }
    j=j+1
  }
}

## QC filtering ####
# Handle QC #
# {
#   rel.sd = list()
#   for(i in unique(sample_cohort)){
#     temp = s5 %>%
#       dplyr::select(sample_names[which(sample_cohort == i)]) %>%
#       mutate(rel.sd = rowSds(as.matrix(.))/rowMeans(.))
#     rel.sd[[i]] = temp$rel.sd
#   }
#   s6 = bind_cols(s5, rel.sd) %>%
#     mutate(QC_flag = rowMeans(.[unique(sample_cohort)]) < 0.3) %>%
#     filter(QC_flag)
# }
# Handle Sample based on QC # 
# {
#   QC = read.csv("hmdb_QC.csv")
#   s6 = s5 %>%
#     mutate(QC_filter = HMDB_ID %in% QC$HMDB_ID) %>%
#     filter(QC_filter)
# }

## Filter ####
{
  s6 = s5 %>%
    arrange(ID) %>%
    filter(flag) %>% 
    filter(!High_blank) %>%
    filter(T) %>%
    dplyr::select(-flag)
}

## Output files for MetaboAnalyst ####
if(min(table(sample_cohort)) >= 3 & length(table(sample_cohort)) > 1){
  {
    MA_output = s6 %>%
      dplyr::select(ID, sample_names)
    MA_output = rbind(c("Cohort", sample_cohort), MA_output)
    MetaboAnalyst_filename = paste("MA_", filename, sep="")
    write.csv(MA_output, file=MetaboAnalyst_filename, row.names=F)
  }
}

## Output files for MetaboAnalyst ####
if(min(table(sample_cohort)) >= 3 & length(table(sample_cohort)) > 1 & require(MetaboAnalystR)){
  {
    
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

  }
}


## Output formating ####
{
  hmdb_match_output=s6 
  
  if(min(table(sample_cohort)) >= 3 & length(table(sample_cohort)) > 1 & exists("FDR_file")){
    FDR_raw <- read_csv(FDR_file) %>%
      dplyr::select(c("X1","FDR")) %>%
      dplyr::rename(ID = X1) %>%
      mutate(FDR = -log10(FDR))
    
    hmdb_match_output = hmdb_match_output %>%
      merge(FDR_raw, by = "ID", all.x = T) %>%
      mutate(FDR = ifelse(is.na(FDR), 0, FDR))
  } else {
    hmdb_match_output = hmdb_match_output %>%
      mutate(FDR = 0)
  }
  
  refcols = c("ID","medMz","medRt","HMDB_ID", "Formula","Metabolite","FDR","High_blank")
  hmdb_match_output = hmdb_match_output %>%
    dplyr::select(refcols, all_names) %>%
    arrange(-FDR, ID) %>%
    dplyr::rename(`-log(FDR)` = FDR)
  
  write.csv(hmdb_match_output, file=paste("HMDB_",filename,sep=""), row.names=F)
}

# Clean up
{
  # fn <-MetaboAnalyst_filename
  # if (file.exists(MetaboAnalyst_filename)) file.remove(fn)
  if(exists("FDR_file")){
    if(file.exists(FDR_file)) file.remove(FDR_file)
  }
  
}


# 
# ## Data pre-clean ####
# {
#   mdata_pre_clean = hmdb_match_output %>%
#     filter(Metabolite != "") %>%
#     mutate(row_label = paste(stri_sub(Metabolite,1,30), medRt, sep="_")) %>%
#     arrange(-FDR)
# 
#   row_labels = mdata_pre_clean$row_label
#   row_labels = NULL
# 
#   mdata_pre_clean = mdata_pre_clean %>%
#     dplyr::select(sample_names)
# 
#   mdata_clean = mdata_pre_clean %>%
#     data_impute(impute_method = "threshold", random = F) %>%  # "min", "percentile", "threshold"
#     data_normalize_row(nor_method = "row_median") %>% # row_median, row_mean, sample, cohort
#     data_normalize_col(nor_method = "") %>% # col_median, col_sum
#     data_transform(transform_method = "log2") %>%
#     data_scale(scale_method = "mean_center") # mean_center, auto(mean_center + divide sd)
# 
#   # method "ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"
#   # distance "correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
#   # correlation is default as pearson correlation
#   row_cluster = cluster_mat(mdata_clean, method = "complete", distance = "correlation")
#   col_cluster = cluster_mat(t(mdata_clean), method = "complete", distance = "correlation")
# 
#   # col_cluster = F
# }
# 
# 
# ## Heat map ####
# {
#   my_plot_heatmap(raw_data = mdata_clean,
#                   cohort = sample_cohort,
#                   row_labels = row_labels, 
#                   imgName = sub(".csv","", filename), 
#                   format = "pdf", # pdf
#                   dpi = 72,
#                   width = NA, # define output graph width
#                   palette = "RdBu",  # RdBu, gbr, heat, topo
#                   viewOpt = "overview", # detail, overview
#                   rowV = row_cluster, # cluster by row, "F" if not clutser
#                   colV = col_cluster, # cluster by column, "F" if not clutser
#                   border = T, # border for each pixel
#                   grp.ave = F, # group average
#                   scale_ub = 3, scale_lb = -3 # heatmap scale, auto if ub = lb
#   )
# }
# 

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

## Temp code for merging NetID workflow####
# {
#   mdata_pre_clean = read.csv("merge.csv") %>%
#     filter(category == "Metabolite") %>%
#     filter(FDR > 2.5) %>%
#     filter(T)
#   
#   row_labels = paste(mdata_pre_clean$ID, mdata_pre_clean$formula, mdata_pre_clean$medRt)
#   row_labels = paste(mdata_pre_clean$ID, round(mdata_pre_clean$medMz,4), mdata_pre_clean$medRt, sep="_")
#   # row_labels = NULL
#   
#   mdata_clean = mdata_pre_clean %>%    
#     dplyr::select(sample_names) %>%
#     data_impute(impute_method = "threshold", random = F) %>%  # "min", "percentile", "threshold"
#     data_normalize_row(nor_method = "") %>% # row_median, row_mean, sample, cohort
#     data_normalize_col(nor_method = "") %>% # col_median, col_sum
#     data_transform(transform_method = "log2") %>%
#     data_scale(scale_method = "mean_center") # mean_center, auto(mean_center + divide sd)
#   
#   
#   row_cluster = cluster_mat(mdata_clean, method = "complete", distance = "euclidean")
#   col_cluster = cluster_mat(t(mdata_clean), method = "complete", distance = "euclidean")
#   
# }