## Library ####
{
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringi)
  library(matrixStats)
  library(dplyr)
  library(rstudioapi)
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}


# Read files
{
  setwd("C:/Users/Li Chen/Desktop/Github local/myrepo/NetID_1204/matt")
  
  
  NetID_file = "matt.RData"
  load(NetID_file)
  stat_filename = "HMDB_raw_data.csv"
  stat_df = read.csv(stat_filename)
  netid_df = ilp_nodes_annotation %>%
    filter(ilp_result > 0.1)
  
  netid_annotation = netid_df %>%
    # select(Input_id.y, everything()) %>%
    select(Input_id, log10_inten, formula, class, path, category) %>%
    filter(T)
  
  result = merge(stat_df, netid_annotation, by.x = "ID", by.y = "Input_id") %>%
    select(1:3, log10_inten, formula, class, category, path, everything() )
  
  write.csv(result, "merge_NetID.csv", row.names = F, na="")
}

# 
{
  
  
}

