

# This code batch run NetID_run_script code 
# cd to the main folder, i.e. Meta_mine

main_dir = getwd()
while(basename(main_dir)!="Meta_mine"){
  if(!grepl("Meta_mine", main_dir)){
    setwd("./Meta_mine")
  } else {
    setwd("..")
  }
  main_dir = getwd()
}

foldernames = list.dirs()

batch_run_script = c()
for(i in 1:length(foldernames)){
  if(!any(list.files(foldernames[i]) == "final.RData") & any(list.files(foldernames[i]) == "raw_data.csv")){
    print(foldernames[i])
    batch_run_script = c(batch_run_script, foldernames[i])
  }
  setwd(main_dir)
}

select_run_script = batch_run_script

for(i in 1:length(select_run_script)){
  work_dir = select_run_script[i]

  source("./code/NetID_function.R")
  if(grepl("pos",select_run_script[i])) {
    ion_mode = 1
  } else if(grepl("neg",select_run_script[i])) {
    ion_mode = -1
  }
  print(select_run_script[i])
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
  
  ## Initialise ####
  {
    Mset[["Global_parameter"]]=  list(mode = ion_mode,
                                      normalized_to_col_median = F)
    Mset[["Cohort"]]=Cohort_Info(Mset)
    print(Mset$Cohort)
    
    #Clean-up duplicate peaks 
    Mset[["Data"]] = Peak_cleanup(Mset,
                                  ms_dif_ppm=5/10^6, 
                                  rt_dif_min=0.1,
                                  detection_limit=2500)
    
    Mset[["ID"]] = Mset$Data$ID
  }
  
  setwd(main_dir)
}




