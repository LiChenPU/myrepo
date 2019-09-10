

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("C:/study/data/exactive/Yeast_unknown/190906 Yeast unknowns + spike in")
filenames = list.files(pattern="13C_d7NH4OAc", recursive = T)


for(i in filenames){
  newname = i
  newname = gsub("13C_d7NH4OAc(.*)","12C_d7NH4OAc_M\\1",i)
  
  file.rename(i, newname)
}
