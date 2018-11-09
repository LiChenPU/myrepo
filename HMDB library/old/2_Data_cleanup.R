
library(readr)
#install.packages("rcdk")
library(rcdk)
library(dplyr)
library(tidyr)

time = data.frame(Sys.time())
colnames(time)="read_data"

setwd("C:/Users/Li Chen/Dropbox/HRMass_ID/180629 Triwizard Tournament/HMDB library")

#data1 records select data from raw table.
df1 <- read_csv("hmdb_summary.csv")
#df1 = df1[1:1000,]

#Remove unwanted elements
{
  unwanted_elements=c(#"F", "Cl", "Br", "As", "Se", "I", 
                      "Mg", "Na", "Al", "Si","Gd", "Te", "Fe",
                      "Mo", "Zn", "K", "Ga", "Ag", "Bi", "Sn", "Pt", "Ca", "Au","Co", "Hg", "Cu","B")
  unwanted_symbols=c(")")
  unwanted = c(unwanted_elements,unwanted_symbols)
  elements_list = sapply(unwanted, grep, df1$MF)
  elements_unlist = unlist(elements_list)
  data_known_no_elements = df1[-elements_unlist,]
  df2 = data_known_no_elements
}

#Functions to parse chemical formula
formula_table = function(X){
  Xformula = get.formula(X)
  temp=Xformula@isotopes
  df = data.frame(as.numeric(temp[,2]))
  rownames(df)=temp[,1]
  colnames(df)=Xformula@string
  return(df)
}
table_to_formula = function(result_formula){
  r = as.character()
  
  for(i in 1:nrow(result_formula)){
    if(result_formula$result[i]!=0){
      if(result_formula$result[i]==1){r=paste(r, result_formula$Row.names[i], sep="")}
      else{r=paste(r, result_formula$Row.names[i],  result_formula$result[i], sep="")}
    }
  }
  if(sum(result_formula$result<0)!=0){
    r=paste(r, "Illeagal_formula")
  }
  return (r)
}
formula_manipulate = function(F1, F2, plusminus) {
  temp_1 = formula_table(F1)
  temp_2 = formula_table(F2)
  result_formula = merge(temp_1,temp_2, all=T,  by="row.names")
  result_formula[is.na(result_formula)]=0
  result_formula["result"]=result_formula[,2]+result_formula[,3]*plusminus
  return (table_to_formula(result_formula))
}

#Convert formula to neutral molecules
{
  #df2=df1
  df_charge=df2
  for (i in 1:nrow(df_charge)){
    if(df_charge$Charge[i]!=0){
      df_charge$MF[i] = formula_manipulate(df_charge$MF[i], "H", -df_charge$Charge[i])
      df_charge$Exact_Mass[i] = get.formula(df_charge$MF[i])@mass
    }
  }
  df3 = df_charge
}

# write.table(df3, file = "hmdb_full.tsv", sep = "\t",
#             col.names = T, row.names = F, quote = F)

#Remove duplicated formula entries
{
  df_unique = df3
  df_unique = df_unique[!duplicated(df_unique$MF),]
  df4 = df_unique
}

# #Validate legit formula
# {
#   df_valid = df4
#   df_valid["valid_formula"]=FALSE
#   for (i in 1:nrow(df_valid)){
#     df_valid$valid_formula = isvalid.formula(get.formula(df_valid$MF[i]))
#   }
#   df5 = df4[df_valid$valid_formula,]
# }

# # Validate Exact mass
# {
#   df_pred_mass=df4[df4$Charge==0,]
#   df_pred_mass["pred_mass"]=0
#   df_pred_mass["mass_dif"]=0
#   for (i in 1:nrow(df_pred_mass)){
#     df_pred_mass$pred_mass[i]=get.formula(df_pred_mass$MF[i])@mass
#     df_pred_mass$mass_dif[i]=df_pred_mass$pred_mass[i]-df_pred_mass$Exact_Mass[i]
#   }
# }

write.table(df4, file = "hmdb_unique.tsv", sep = "\t",
            col.names = T, row.names = F, quote = F)

#df4 <- read_csv("hmdb_clean.csv")

# Convert formula into element table
{
  full_formula = formula_table("CHNOPS")
  full_formula["CHNOPS"]=row.names(full_formula)
  for (i in 1:nrow(df4)){
    temp_formula = formula_table(df4$MF[i])
    temp_formula["CHNOPS"]=row.names(temp_formula)
    full_formula=merge(temp_formula, full_formula, by="CHNOPS", all = T)
  }
  full_formula[is.na(full_formula)]=0
  full_formula_t = t(full_formula)
}

write.table(full_formula_t, file = "full_formula.csv", sep = ",",
            col.names = F, row.names = T, quote = F)

full_formula <- read_csv("full_formula.csv")

#Merge data table, remove NA
{
  colnames(full_formula)[colnames(full_formula)=="CHNOPS"]="MF"
  df5 = merge(full_formula, df4, by="MF", all=T)
  df6 = df5[,c(8:9,1:7, 13:14, 20,26,21:25)]
  df6 = df6[complete.cases(df6),]
  df6 = df6[with(df6, order(HMDB_ID)),]
}

#Reorganized features
{
  df_reorganize=df6
  df_reorganize["Unsaturation_degree"]=df_reorganize$C+1-(df_reorganize$H-df_reorganize$N-df_reorganize$P)/2
  df_reorganize["HPO3"] = df_reorganize$P
  df_reorganize["O_nonHPO3"] = df_reorganize$O-df_reorganize$P*3
  df_reorganize["H_nonHPO3"] = df_reorganize$H-df_reorganize$P
  df_reorganize=df_reorganize[df_reorganize$C!=0,]
  df7=df_reorganize
}

write.table(df7, file = "input.tsv", sep = "\t",
            col.names = T, row.names = F, quote = F)


time["end"]=Sys.time()
time_used=as.numeric()
for (i in 2:ncol(time)){
  time_used[i-1] = round(time[1,i]-time[1,i-1],2)
}
names(time_used)=colnames(time)[1:(ncol(time)-1)]
time_used
time
