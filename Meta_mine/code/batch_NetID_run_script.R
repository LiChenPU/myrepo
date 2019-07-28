

main_dir = getwd()
# foldernames = list.files(recursive=T)[!grepl("\\.", list.files(recursive=T))]
foldernames = list.dirs()
select = c()
i=55
for(i in 1:length(foldernames)){
  # setwd(paste("./",foldernames[i],sep=""))
  
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

