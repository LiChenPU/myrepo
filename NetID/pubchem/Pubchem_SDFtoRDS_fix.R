# BiocManager::install("ChemmineR")
library(ChemmineR)
library(dplyr)

filenames = list.files(pattern = ".sdf.gz")

sink("log.txt")
print(timestamp())
for(i in 1:length(filenames)){
  temp_rds = readRDS(sub("sdf.gz","rds", filenames[i]))
  if(colnames(temp_rds)[1] != "CMP1"){next} # skip when rds file is good
  print(paste("Processing", filenames[i]))
  dsfset = read.SDFset(filenames[i])
  valid <- validSDF(dsfset) # Identifies invalid SDFs in SDFset objects 
  dsfset <- dsfset[valid] # Removes invalid SDFs, if there are any 
  
  temp_result_ls = datablock(dsfset)
  
  temp_result = bind_rows(!!!temp_result_ls) # Lists are treated as data frames but can be spliced explicitly with !!!
  
  saveRDS(temp_result, sub("sdf.gz","rds", filenames[i]))
}

sink()