library(rjson)
library(dplyr)

#test2

target=373.1214



charge = "-"
massRange = 5*target/10^6
time = Sys.time()
mfRange="C10-99 H20-100 N4-12 O6-20 P0-3 S0-3 B0-1 Si0-2 K0-3 Cl0-5 Na0-3 Cu0-1 Cr0-1 Ca0-1 Ni0-1"
# mfRange="C6-20 H10-100 N4-7 O2-20 P0-3 S0-3 Si1-2 K0-4 Cl0-5 Na0-3 Cr0-1 Ca0-1 Ni0-1 B0-1"
# mfRange="C0-20 H0-100 N0-7 O0-20 P0-3 S0-3 Si0-2 K0-4 Cl0-5 Na0-3 Cu0-1 Cr0-1 Ca0-1 Ni0-1 B0-1 Li0-1 Fe0-1 Mg0-1 Be0-1 F0-6 Al0-1"
# mfRange="C0-20 H0-100 N0-7 O0-20 P0-3 S0-3 Si0-2 K0-4 Cl0-5 Na0-3 Cu0-1 Cr0-1 Ca0-1 Ni0-1 B0-1"
mfRange = paste(mfRange, charge, sep="")

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
