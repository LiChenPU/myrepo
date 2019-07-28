

# This code batch run NetID_run_script code 
# cd to the main folder, i.e. Meta_mine

main_dir = getwd()

foldernames = list.dirs()

for(i in 1:length(foldernames)){
  if(!any(list.files(foldernames[i]) == "final.RData") & any(list.files(foldernames[i]) == "raw_data.csv")){
    print(foldernames[i])
    work_dir = foldernames[i]
    
    source("./code/NetID_function.R")
    if(grepl("pos",foldernames[i])) {
      ion_mode = 1
    } else if(grepl("neg",foldernames[i])) {
      ion_mode = -1
    }
    source("./code/NetID_run_script.R")
  }
  setwd(main_dir)
}

