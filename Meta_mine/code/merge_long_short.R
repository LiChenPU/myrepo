

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("C:/Users/Li Chen/Desktop/matt")

work_dir = getwd()

# folder = list.files()[!grepl("\\.", list.files())]
folder = list.files()

# i =2
# for(i in 1:length(folder)){
  # setwd(paste("./",folder[i],sep=""))
  long = read.csv("pos_2.csv", stringsAsFactors = F)
  short = read.csv("pos_1.csv", stringsAsFactors = F)
  
  raw_data = merge(long, short,all=T)
  
  col_names = colnames(raw_data[,15:ncol(raw_data)])
  col_names = col_names[order(col_names)]
  raw_data = cbind(raw_data[,1:14], raw_data[,col_names])
  raw_data = raw_data %>%
    mutate(groupId = 1:nrow(.))
  write.csv(raw_data,"raw_data.csv", row.names = F, na="")
  print(paste(nrow(long), nrow(short), nrow(raw_data)))
  # setwd(work_dir)
# }
