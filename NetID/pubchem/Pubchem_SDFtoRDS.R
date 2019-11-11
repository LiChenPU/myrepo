# BiocManager::install("ChemmineR")
library(ChemmineR)
library(dplyr)

filenames = list.files(pattern = ".sdf.gz")

sink("log.txt")
print(timestamp())
for(i in 1:length(filenames)){
  if(file.exists(sub("sdf.gz","rds", filenames[i]))){next} # skip when rds file exists
  print(paste("Processing", filenames[i]))
  dsfset = read.SDFset(filenames[i])
  valid <- validSDF(dsfset) # Identifies invalid SDFs in SDFset objects 
  dsfset <- dsfset[valid] # Removes invalid SDFs, if there are any 
  
  temp_result_ls = datablock(dsfset)
  temp_result = rbind_list(temp_result_ls)
  
  saveRDS(temp_result, sub("sdf.gz","rds", filenames[i]))
}

sink()