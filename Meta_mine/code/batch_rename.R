
# This code creates folder for each csv file using the csv name
# copy the csv to the new folder
# Rename the csv to raw_data for NetID processing

setwd(work_dir)
filenames = list.files(pattern = ".csv")
foldernames = sub("\\.csv", "", filenames)

for(i in 1:length(filenames)){
  setwd(work_dir)
  raw_data = read.csv(filenames[i])
  new_folder = paste("./", foldernames[i],sep="")
  dir.create(new_folder)
  setwd(new_folder)
  write.csv(raw_data, "raw_data.csv", row.names = F)
}
