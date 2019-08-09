

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("C:/Users/lc8/Documents/GitHub/myrepo/Meta_mine/library/mice_NDUFS4_neg")
work_dir = getwd()

folder = list.files()[!grepl("\\.", list.files())]

i =2
for(i in 1:length(folder)){
  setwd(paste("./",folder[i],sep=""))
  long = read.csv("long.csv", stringsAsFactors = F)
  short = read.csv("short.csv", stringsAsFactors = F)
  
  raw_data = merge(long, short,all=T)
  
  col_names = colnames(raw_data[,15:ncol(raw_data)])
  col_names = col_names[order(col_names)]
  raw_data = cbind(raw_data[,1:14], raw_data[,col_names])
  write.csv(raw_data,"raw_data.csv", row.names = F)
  print(paste(nrow(long), nrow(short), nrow(raw_data)))
  setwd(work_dir)
}
