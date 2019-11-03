
# This code creates folder for each csv file using the csv name
# copy the csv to the new folder
# Rename the csv to raw_data for NetID processing

# setwd(work_dir)
setwd("C:/Users/lc8/Documents/myrepo/NetID/Wenyun_Yeast_neg")
setwd("C:/Users/lc8/Documents/myrepo/NetID/Wenyun_Yeast_pos")
setwd("C:/Users/lc8/Documents/myrepo/NetID/Lin_Yeast_neg")
# setwd("C:/Users/lc8/Documents/myrepo/NetID/Lin_Yeast_pos")


filenames = list.files(pattern = "out.csv")
foldernames = sub("\\.csv", "", filenames)

i=1
# for(i in 1:length(filenames)){
  # setwd("C:/Users/lc8/Documents/myrepo/NetID")
  raw_data = read.csv(filenames[i])
  # new_folder = paste("./", foldernames[i],sep="")
  # dir.create(new_folder)
  # setwd(new_folder)
  colnames(raw_data)
  nrow(raw_data)
  
  
  ## Specific to connect PAVE output
  {
    # colnames(raw_data)[27:28]=c("Blank1", "Blank2")
    # raw_data = raw_data[,c(1:17, 36:37)]
    # colnames(raw_data) = gsub("control","blank", colnames(raw_data))
    raw_data = raw_data[,c(1:17, 27:28)]
    colnames(raw_data) = gsub("posi|pos","neg", colnames(raw_data))
    
  }
  
  
  write.csv(raw_data, "raw_data.csv", row.names = F)
# }
