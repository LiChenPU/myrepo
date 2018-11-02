library(rjson)
library(dplyr)

#test2

target=906.9958
massRange = 5*906.9958/10^6
time = Sys.time()
mfRange="C0-100H0-100N0-12O0-20S0-5P0-5S0-5"

query <- paste("http://www.chemcalc.org/chemcalc/em?monoisotopicMass=",target,
               "&mfRange=", mfRange,
               "&massRange=", massRange,sep=""
               )
resJSON <- fromJSON(paste(readLines(query), collapse=""))



temp = data.frame(mf="mf", em = 0, error = 0, stringsAsFactors = F)
temp_ls = list()
for(i in 1:1){
for (mf.index in 1:length(resJSON$results)) {
  
  temp[1,] = list( "mf"=resJSON$results[[mf.index]]$mf , 
                         "em"=resJSON$results[[mf.index]]$em, 
                         "error"=resJSON$results[[mf.index]]$error )
  temp_ls[[mf.index]]=temp
}
}
formula_ls2 = bind_rows(temp_ls)
time = Sys.time()-time
time
  

