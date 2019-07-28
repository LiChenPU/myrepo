
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
work_dir = getwd()

filenames = list.files(patter = ".csv")
foldernames = sub("\\.csv", "", filenames)


for(i in 1:length(filenames)){
  setwd(work_dir)
  raw_data = read.csv(filenames[i])
  new_folder = paste("./", foldernames[i],sep="")
  dir.create(new_folder)
  setwd(new_folder)
  write.csv(raw_data, "raw_data.csv", row.names = F)
}
