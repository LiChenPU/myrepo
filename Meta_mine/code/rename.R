

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("./Full Scans")
filenames = list.files(pattern=".mzXML")


for(i in filenames){
  newname = i
  newname = sub("336", "KO_336_M", newname)
  newname = sub("339", "KO_339_F", newname)
  newname = sub("340", "KO_340_F", newname)
  newname = sub("337", "Het_337_M", newname)
  newname = sub("341", "Het_341_F", newname)
  newname = sub("342", "Het_342_F", newname)
  newname = sub("347", "Het_347_M", newname)
  newname = sub("338", "WT_338_M", newname)
  newname = sub("345", "WT_345_F", newname)
  newname = sub("348", "WT_348_M", newname)
  
  
  file.rename(i, newname)
}
