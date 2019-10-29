
# This code creates folder for each csv file using the csv name
# copy the csv to the new folder
# Rename the csv to raw_data for NetID processing

# setwd(work_dir)
setwd("C:/Users/lc8/Documents/myrepo/NetID")
filenames = list.files(pattern = ".csv")
foldernames = sub("\\.csv", "", filenames)

for(i in 1:length(filenames)){
  setwd("C:/Users/lc8/Documents/myrepo/NetID")
  raw_data = read.csv(filenames[i])
  new_folder = paste("./", foldernames[i],sep="")
  dir.create(new_folder)
  setwd(new_folder)
  
  
  ## Specific to connect PAVE output
  {
    # colnames(raw_data)[27:28]=c("Blank1", "Blank2")
    raw_data = raw_data[,c(1:17, 27:28)]
    
  }
  
  
  write.csv(raw_data, "raw_data.csv", row.names = F)
}
