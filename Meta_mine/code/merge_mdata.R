
library(readr)
library(tidyr)
library(dplyr)
library("stringr")

# functions ####
## Cluster dataset to see if LC condition are similar ####
DatasetDist = function(dt1, dt2, 
                       log10_inten_cutoff = 4,  # only compare peaks with intensity above the threshold
                       merge_group_ppm_tol = 3 # mz within ppm_tol will merge into one mz group
                       )
{
  # dt1 = raw_ls[[1]] # liver in QE+
  # dt2 = raw_ls[[2]] # liver QE2
  # dt3 = raw_ls[[12]] # quad QE2
  
  if(identical(dt1, dt2)){return(0)}
  
  
  dt1_clean = dt1 %>%
    dplyr::select(ID, medMz, medRt, log10_inten) %>%
    filter(!duplicated(ID)) %>%
    filter(log10_inten > log10_inten_cutoff) %>%
    arrange(medMz)
  
  
  dt2_clean = dt2 %>%
    dplyr::select(ID, medMz, medRt, log10_inten) %>%
    filter(!duplicated(ID)) %>%
    filter(log10_inten > log10_inten_cutoff) %>%
    arrange(medMz)
  
  
  dt1_clean_medMz = dt1_clean$medMz
  nrow_dt1_clean = nrow(dt1_clean)
  
  
  i = 1
  count = 1
  dt1_clean_ls = list()
  while(i<=nrow_dt1_clean){
    temp_mz = dt1_clean_medMz[i]
    i_max = i
    while(i_max < nrow_dt1_clean){
      if(dt1_clean_medMz[i_max+1] - dt1_clean_medMz[i_max] <= dt1_clean_medMz[i_max] * merge_group_ppm_tol/1E6){
        i_max = i_max+1
      } else {
        break
      }
    }
    temp_dt1 = dt1_clean[i:i_max,]
    temp_df2_filter = dt2_clean$medMz>dt1_clean_medMz[i]*(1-merge_group_ppm_tol/1E6) & 
      dt2_clean$medMz<dt1_clean_medMz[i_max]*(1+merge_group_ppm_tol/1E6)
    temp_dt2 = dt2_clean[temp_df2_filter,]
    dt1_clean_ls[[count]] = list(dt1 = temp_dt1, dt2 = temp_dt2)
    i = i_max+1
    count = count+1
  }
  
  library(rdist)
  library(clue)
  test = unlist(lapply(dt1_clean_ls, function(dt_clean){
    # temp_dt1 = as.matrix(dt_clean$dt1[,4:5])
    # temp_dt2 = as.matrix(dt_clean$dt2[,4:5])
    temp_dt1 = as.matrix(dt_clean$dt1$medRt)
    temp_dt2 = as.matrix(dt_clean$dt2$medRt)
    if(nrow(temp_dt2)==0){return(NA)}
    
    x=cdist(temp_dt1, temp_dt2)
    if(ncol(x)<nrow(x)){x = t(x)}
    y=solve_LSAP(x)
    all_dist = x[cbind(seq_along(y), y)]
    # all_dist = 1
    return(median(all_dist))
  }))
  
  return(median(test, na.rm = T))
  
}
## Clean up preliminary merge of different mdata file #### 
cleanUp = function(raw)
{
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
  
  # Remove peaks duplicated peaks in one group from the same run
  {
    s3 = s2
    s4 = s3[with(s3, order(MZRT_group, mean_inten)),]
    # test = s4[duplicated(s4[,c("formula","filename", "MZRT_group")]) |
    #             duplicated(s4[,c("formula","filename", "MZRT_group")], fromLast = T), ]
    s4 = s4[!duplicated(s4[,c("formula","filename", "MZRT_group")]), ]
  }
  
  # Take median of the mz and rt for peaks with same MZRTgroup
  {
    MZRT_group = s4$MZRT_group
    medMz = s4$medMz
    medRt = s4$medRt
    
    k_max=k_min=1
    while (k_max <= length(MZRT_group)){
      k_min = k_max
      while (MZRT_group[k_min] == MZRT_group[k_max]){
        k_max = k_max+1
        if(k_max > length(MZRT_group)){break}
      }
      medMz[k_min:(k_max-1)]=median(medMz[k_min:(k_max-1)], na.rm = T)
      medRt[k_min:(k_max-1)]=median(medRt[k_min:(k_max-1)], na.rm = T)
    }
    
    s4$medMz = medMz
    s4$medRt = medRt
  }
  return(s4)
}

# Main ####
## load files and set parameters ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("../library/pig_blood_pos")
  
  filenames = list.files(recursive = T, pattern = "mdata.csv")
  filenames = filenames[grepl("\\/", filenames)]
  
  
  #Consider two groups have the same m/z if difference is below threshold
  ms_dif_ppm= 5
  ms_dif_ppm = ms_dif_ppm/10^6
  #Consider two groups have the same RT if difference is below threshold
  rt_dif_min=0.2
  
  num_of_files = length(filenames)
  raw_ls = list()
  for(i in 1:num_of_files){
    filename=filenames[i]
    df_temp = read_csv(filename)
    df_temp["filename"]= dirname(filenames[i])
    
    ## Clean up data specifically
    {
      
      
      colnames(df_temp) = gsub("(\\d+$)", paste("_\\1",dirname(filenames[i]),sep="_"), colnames(df_temp))
      colnames(df_temp) = gsub("_+", "_", colnames(df_temp))
      
    }
    
    raw_ls[[i]]= df_temp
  }
  
  raw = bind_rows(raw_ls)
  rm(df_temp)
}

## Cluster dataset to see if LC condition are similar ####
{
  # dataset_dist_mtrx = matrix(0,
  #                            nrow = length(filenames),
  #                            ncol = length(filenames)
  # )
  # rownames(dataset_dist_mtrx) = colnames(dataset_dist_mtrx) = dirname(filenames)
  # 
  # for(i in 1:length(filenames)){
  #   print(Sys.time())
  #   for(j in 1:length(filenames)){
  #     dataset_dist_mtrx[i,j] = DatasetDist(raw_ls[[i]], raw_ls[[j]],
  #                                          log10_inten_cutoff = 5,  # only compare peaks with intensity above the threshold
  #                                          merge_group_ppm_tol = 3 # mz within ppm_tol will merge into one mz group
  #     )
  #   }
  # }
  # dataset_dist_mtrx_sym = 0.5 *(dataset_dist_mtrx + t(dataset_dist_mtrx))
  # pdf("median_RT_shift_between_dataset.pdf")
  # pheatmap::pheatmap(dataset_dist_mtrx_sym)
  # dev.off()
}

###### Don't need to change from here #####
#Reorganize peak table
{
  s4 = cleanUp(raw)
  del = c("MZ_group|goodPeakCount|maxQuality|parent")
  s5 = s4[, !grepl(del, colnames(s4), ignore.case = T)]
}

# statTable & intenTable
{
  # statTable
  summaryColNames = c("ID","formula","ILP_result","is_metabolite","_log10_FDR","library_match_formula",
                   "library_match_name","high_blank","medMz","medRt","filename","MZRT_group", 
                   "mean_inten","log10_inten")
  summaryTable = s5[,summaryColNames]
  
  summaryTable_ls = base::split(summaryTable, summaryTable$MZRT_group)
  
  statTable = summaryTable %>%
    dplyr::select(MZRT_group, medMz, medRt, `_log10_FDR`, filename, mean_inten) %>%
    arrange(mean_inten) %>%
    distinct(MZRT_group, filename, .keep_all = T) %>%
    dplyr::select(-mean_inten) %>%
    mutate(filename = paste(filename,"log10", sep="_")) %>%
    spread(filename, `_log10_FDR`)
  
  statTable["Sum_logP"] = apply(statTable[,grepl("log10", colnames(statTable))], 1, sum, na.rm = T)
  
  # formulaTable
  formulaTable = summaryTable %>%
    dplyr::select(MZRT_group, medMz, medRt, formula, filename, mean_inten) %>%
    arrange(mean_inten) %>%
    distinct(MZRT_group, filename, .keep_all = T) %>%
    dplyr::select(-mean_inten) %>%
    mutate(filename = paste(filename,"formula", sep="_")) %>%
    spread(filename, formula)
  

  
  # intenTable
  intenTable = s5[, c("MZRT_group", colnames(s5)[!colnames(s5)%in%summaryColNames])]
  intenTable = intenTable %>%
    gather(cohorts, value = number, -MZRT_group, na.rm = T) %>%
    arrange(desc(number)) %>%
    distinct(MZRT_group, cohorts, .keep_all=T) %>%
    spread(cohorts, number)
  
  # cohortTable
  sample_names = colnames(intenTable)[-1]
  blank_names = sample_names[grepl("blank|blk", sample_names, ignore.case = T)]
  
  test = str_split(sample_names, "_")
  n.obs <- sapply(test, length)
  seq.max <- seq_len(max(n.obs))
  mat <- t(sapply(test, "[", i = seq.max))
  
  cohortTable = as.data.frame(mat, stringsAsFactors =F)
  rownames(cohortTable) = sample_names
  cohortTable[grepl("blank|blk", rownames(cohortTable), ignore.case = T),] = "blank"
}

#Output merge_mdata
{
  result_all = merge(statTable,intenTable)
  result_all = merge(formulaTable,result_all)
  result_all = result_all[with(result_all, order(-Sum_logP)),]
  merge_mdata = result_all
  write.csv(merge_mdata, "merge_mdata.csv", row.names = F)
}

save(summaryTable, merge_mdata, intenTable, file = "merge_mdata.RData")




