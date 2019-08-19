#显示中文 
#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
library(readr)
library(tidyr)
library(dplyr)
library("stringr")

#Set path here
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("./10tissues_pos_mode")


foldernames = list.files()[!grepl("\\.", list.files())]

#Consider two groups have the same m/z if difference is below threshold
ms_dif_ppm= 5
ms_dif_ppm = ms_dif_ppm/10^6
#Consider two groups have the same RT if difference is below threshold
rt_dif_min=0.2
i=1

num_of_files = length(foldernames)
raw_ls = list()
for(i in 1:num_of_files){
  foldername=foldernames[i]
  df_temp = read_csv(paste("./",foldername, "./", foldername,".csv", sep=""))
  df_temp = df_temp[!is.na(df_temp$medMz),]
  df_temp$label= foldername
  raw_ls[[i]]= df_temp
}
rm(df_temp)




###### Don't need to change from here #####
#Read data
{
raw = bind_rows(raw_ls)
ncol_raw = ncol(raw)
nrow_raw = nrow(raw)
}

##Group MS groups
{
s = raw[with(raw, order(medMz, medRt)),]
s[["merge_group"]]=NA
mgMS_count = 1 
i_max=i_min=1
while (i_min <= nrow_raw){
  while(s$medMz[i_max]-s$medMz[i_min]< (s$medMz[i_min]*ms_dif_ppm)){
    i_max = i_max+1
    if(i_max>nrow_raw){
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
# profvis::profvis(for(i in 1:1){})
{
  k_max=k_min=1
  while (k_max <= nrow(s3)){
    k_min = k_max
    while (s3$RTmerge_group[k_min] == s3$RTmerge_group[k_max]){
      k_max = k_max+1
      if(k_max > nrow(s3)){break}
    }
    s3$medMz[k_min:(k_max-1)]=median(s3$medMz[k_min:(k_max-1)], na.rm = T)
    s3$medRt[k_min:(k_max-1)]=median(s3$medRt[k_min:(k_max-1)], na.rm = T)
    temp = s3[k_min:(k_max-1),14:ncol_raw]
    temp[1,] = apply(temp, 2, function(x){
      if(any(!is.na(x))){
        return(max(x, na.rm = T))
      } else {
        return(NA)
      }
    })
    s3[k_min:(k_max-1), 14:ncol_raw] = temp[1,]
  }
  s4 = s3[!duplicated(s3[,"RTmerge_group"]), ]
}

#Reorganize peak table
{
 result_all = s4[,1:ncol_raw]
}


#Output result

write.csv(result_all, "raw_data.csv", row.names = F)


