

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("C:/Users/lc8/Desktop/lifeng")
work_dir = getwd()

folder = list.files()[!grepl("\\.", list.files())]

i =2
for(i in 9:length(folder)){
  setwd(paste("./",folder[i],sep=""))
  long = read.csv("long.csv", stringsAsFactors = F)
  short = read.csv("short.csv", stringsAsFactors = F)
  raw_data = merge(long, short,all=T)
  
  col_names = colnames(raw_data[,15:ncol(raw_data)])
  col_names = col_names[order(col_names)]
  raw_data = cbind(raw_data[,1:14], raw_data[,col_names])
  write.csv(raw_data,"raw_data.csv", row.names = F)
  setwd(work_dir)
}
