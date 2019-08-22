
# This code batch run NetID_CPLEX code 
# cd to the main folder, i.e. Meta_mine
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
main_dir = getwd()
while(basename(main_dir)!="Meta_mine"){
  if(!grepl("Meta_mine", main_dir)){
    setwd("./Meta_mine")
  } else {
    setwd("..")
  }
  main_dir = getwd()
}

foldernames_batch_CPLEX = list.dirs()

batch_CPLEX = c()
for(i in 1:length(foldernames_batch_CPLEX)){
  setwd(foldernames_batch_CPLEX[i])
  if(!any(list.files() == "mdata.csv") & any(list.files() == "final.RData")){
    batch_CPLEX = c(batch_CPLEX, foldernames_batch_CPLEX[i])
  }
  setwd(main_dir)
}

print(batch_CPLEX)
select_batch_CPLEX = batch_CPLEX
print(select_batch_CPLEX)

for(i in 1:length(select_batch_CPLEX)){
  setwd(select_batch_CPLEX[i])
  print(select_batch_CPLEX[i])
  load("final.RData")
  main_dir = getwd()
  while(basename(main_dir)!="Meta_mine"){
    main_dir = dirname(main_dir)
  }
  source(paste(main_dir, "/code/NetID_CPLEX.R", sep=""))
  setwd(main_dir)
}
