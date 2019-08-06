

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
  source("./code/NetID_run_script.R", local = T)
  setwd(main_dir)
}




