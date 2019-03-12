library(rjson)
library(dplyr)

#test2

target=906.9958

massRange = 5*906.9958/10^6
time = Sys.time()
mfRange="C0-99 H0-100 N0-12 O0-20 P0-3 S0-3 Si0-2 K0-1 Cl0-5 Na0-1"

query <- paste("http://www.chemcalc.org/chemcalc/em?monoisotopicMass=",target,
               "&mfRange=", mfRange,
               "&massRange=", massRange,sep=""
               )
resJSON <- fromJSON(paste(readLines(query), collapse=""))



temp = data.frame(mf="mf", em = 0, error = 0, stringsAsFactors = F)
temp_ls = list()

for (mf.index in 1:length(resJSON$results)) {

  temp[1,] = list( "mf"=resJSON$results[[mf.index]]$mf ,
                         "em"=resJSON$results[[mf.index]]$em,
                         "error"=resJSON$results[[mf.index]]$error )
  temp_ls[[mf.index]]=temp
}
formula_ls2 = bind_rows(temp_ls)
time = Sys.time()-time
time
