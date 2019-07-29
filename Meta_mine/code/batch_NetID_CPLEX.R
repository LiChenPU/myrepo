
# This code batch run NetID_CPLEX code 
# cd to the main folder, i.e. Meta_mine

main_dir = getwd()
while(basename(main_dir)!="Meta_mine"){
  setwd("..")
  main_dir = getwd()
}

foldernames_batch_CPLEX = list.dirs()

for(i in 1:length(foldernames_batch_CPLEX)){
  setwd(foldernames_batch_CPLEX[i])
  if(!any(list.files() == "mdata.csv") & any(list.files() == "final.RData")){
    print(getwd())
    load("final.RData")
    source(paste(main_dir, "/code/NetID_CPLEX.R", sep=""))
  }
  setwd(main_dir)
}

