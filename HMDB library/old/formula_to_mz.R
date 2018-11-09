
#显示中文 
#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
library(readr)
library(tidyr)
library(rcdk)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#df=read_csv("hmdb_unique.csv")
df=read_csv("known.csv")
df=complete.cases(df)
df["mz"]=as.numeric()


for(i in 1: nrow(df)){
  if(!is.na(df$MF[i])){
    df$mz[i]=get.formula(df$MF[i])@mass
  }
}

#df=df[,-5]
#df["test"]=df$Exact_Mass-df$mz


#write.csv(df, file="hmdb_unique_mz.csv", row.names=F)
write.csv(df, file="known_mz.csv", row.names=F)

