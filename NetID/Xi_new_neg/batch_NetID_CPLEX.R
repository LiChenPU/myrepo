
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


for(i in 1:length(NetID_files)){
  
  if(file.exists(sub(".RData", ".csv", NetID_files[i]))){
    next
  } else {
    print(NetID_files[i])
    load(NetID_files[i])
    source("../NetID_CPLEX.R")
  }
    
    
}

library(dplyr)
NetID_files = list.files(pattern = " .csv")

formula_list_ls = list()
for(i in 1:length(NetID_files)){
  formula_list_ls[[length(formula_list_ls)+1]] = read.csv( NetID_files[i], stringsAsFactors = F) %>%
    filter(ILP_result!=0) %>%
    dplyr::select(ID, formula)
  
}



CPLEX_summary = formula_list_ls[[1]]
for(i in 2:length(formula_list_ls)){
  CPLEX_summary = merge(CPLEX_summary,formula_list_ls[[i]], all = T, by.x = "ID", by.y = "ID")
}
colnames(CPLEX_summary) = c("ID", paste("run",1:6,sep=""))

sapply(formula_list_ls, nrow)

CPLEX_summary = CPLEX_summary %>%
  mutate(all_same = run1==run2 & run1==run3 & run1==run4 & run1==run5 & run1==run6)

CPLEX_summary_inconsistent = CPLEX_summary %>%
  filter(!all_same)
  
  
  
  
  

