# Main Codes ####
## Read files ####
print(work_dir)
print(ion_mode)
{
  Mset = list()
  Mset[["Library"]] = read.csv("./code/dependent/HMDB_clean.csv", stringsAsFactors = F)
  Mset[["Biotransform"]]=Read_rule_table(rule_table_file = "./code/dependent/biotransform.csv")
  Mset[["Artifacts"]]=Read_rule_table(rule_table_file = "./code/dependent/artifacts.csv")
  
  setwd(work_dir)
  filename = c("raw_data.csv")
  Mset[["Raw_data"]] <- read_csv(filename)
}


save.image("final.RData")