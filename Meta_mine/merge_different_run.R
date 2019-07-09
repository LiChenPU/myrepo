#显示中文 
#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
library(readr)
library(tidyr)
library(dplyr)
library("stringr")

#Set path here
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("./2liver_quad_QE+pos")


# foldernames = list.files()[!grepl("\\.", list.files())]
filenames = list.files()[grepl(".csv", list.files())&grepl("Mdata", list.files())]


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
  df_temp["id"]= str_extract(filename, "[:alnum:]+(?=.csv?)")
  raw_ls[[i]]= df_temp
}
rm(df_temp)



###### Don't need to change from here #####
#Read data
{
raw = bind_rows(raw_ls)
raw_cn = ncol(raw)
raw_rn = nrow(raw)
}

##Group MS groups
{
s = raw[with(raw, order(medMz, medRt)),]
s[["merge_group"]]=NA
mgMS_count = 1 
i_max=i_min=1
while (i_min <= raw_rn){
  while(s$medMz[i_max]-s$medMz[i_min]< (s$medMz[i_min]*ms_dif_ppm)){
    i_max = i_max+1
    if(i_max>raw_rn){
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
  s2[["RTmerge_group"]]=NA
  
  
  mgRT_count=length(which(s2$merge_group == 0))
  
  s2$RTmerge_group[1:mgRT_count] = 1:mgRT_count
  
  mgRT_count = mgRT_count+1
  summary_RTgroup = list()
  for(i in 1:max(s2$merge_group)){
    loc_i = which(s2$merge_group == i)
    temp_data = s2[loc_i,]
    temp_RT = c(mgRT_count)
    for(j in 2:length(loc_i)){
      if(temp_data$medRt[j] - temp_data$medRt[j-1] >= rt_dif_min){
        mgRT_count = mgRT_count+1
      } 
      temp_RT = c(temp_RT, mgRT_count)
    }
    summary_RTgroup[[i]] = temp_RT
    mgRT_count = mgRT_count+1
  }
  s2$RTmerge_group[(length(which(s2$merge_group == 0))+1): length(s2$RTmerge_group)] = unlist(summary_RTgroup)
}

#Remove peaks occurs only in one run
{
  s3 = s2
#s3 = s2[s2$RTmerge_group!=0,]
}

# s3 = s3[s3$medMz < 397 & s3$medMz >396, ]

#Remove peaks duplicated peaks in one group from the same run
{
#colnames(s3)[colnames(s3)=="[-log10(p)_FDR_corrected]"] = ".log_p_FDR_corrected"
s4 = s3[with(s3, order(RTmerge_group, `_log10_FDR`)),]
  test = s4[duplicated(s4[,c("id", "RTmerge_group")]) |
              duplicated(s4[,c("id", "RTmerge_group")], fromLast = T), ]
s4 = s4[!duplicated(s4[,c("id", "RTmerge_group")]), ]

k_max=k_min=1
while (k_max <= nrow(s4)){
  k_min = k_max
  while (s4$RTmerge_group[k_min] == s4$RTmerge_group[k_max]){
    k_max = k_max+1
    if(k_max > nrow(s4)){break}
  }
  s4$medMz[k_min:(k_max-1)]=median(s4$medMz[k_min:(k_max-1)], na.rm = T)
  s4$medRt[k_min:(k_max-1)]=median(s4$medRt[k_min:(k_max-1)], na.rm = T)
}
}

#Reorganize peak table
{
  del = c("library_match_formula|high_blank|goodPeakCount|maxQuality|library_match_name|blank|inten|parent")
  s5 = s4[, !grepl(del, colnames(s4), ignore.case = T)]
  s5_gather = gather(s5, tissue, value = number, 5:(ncol(s5)), -id, -`_log10_FDR`, -merge_group, -RTmerge_group)
  s5_select = s5_gather[, c(3,4,8,9)]
  s5_select = s5_select[complete.cases(s5_select),]
  s5_spread = spread(s5_select, tissue, number, drop=T )

s5_pvalue = s5[,c("medMz","medRt", "id", '_log10_FDR')]
s5_pvalue_spread = spread(s5_pvalue, id, `_log10_FDR`, drop=T )
s5_pvalue_spread[["sum_log_p"]]=rowSums(s5_pvalue_spread[,3:(2+num_of_files)], na.rm = T)
s5_pvalue_spread[["mean_log_p"]]=rowMeans(s5_pvalue_spread[,3:(2+num_of_files)], na.rm = T)

s5_formula = s4[,c("medMz","medRt", "id", 'library_match_formula')]

s5_formula_spread = spread(s5_formula, id, library_match_formula, drop=T )


object.size(s5_gather)

result_all=merge(s5_pvalue_spread,s5_spread)
result_all = merge(s5_formula_spread,result_all, by = c("medMz", "medRt"))

result_all = result_all[with(result_all, order(-sum_log_p)),]
result_all <- result_all[,colSums(is.na(result_all))<nrow(result_all)]
}


#Output result

write.csv(result_all, "result.csv", row.names = F)


