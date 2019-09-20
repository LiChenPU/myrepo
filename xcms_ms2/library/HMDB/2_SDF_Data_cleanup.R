
library(readr)
#install.packages("rcdk")
library(rcdk)
library(dplyr)
library(tidyr)
library(lc8)
library(enviPat)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#data1 records select data from raw table.
df1 <- read_tsv("HMDB_structure_sdf.tsv")
#df1 = df1[1:1000,]

#Remove unwanted elements
{
  unwanted_elements=c("F", "Cl", "Br", "As", "Se", "I", 
                      "Mg", "Na", "Al", "Si","Gd", "Te", "Fe", "Mn", 
                      "Mo", "Zn", "K", "Ga", "Ag", "Bi", "Sn", "Pt", "Ca", "Au","Co", "Hg", "Cu","B")
  unwanted_symbols=c(")")
  unwanted = c(unwanted_elements,unwanted_symbols)
  elements_list = sapply(unwanted, grep, df1$MF)
  elements_unlist = unlist(elements_list)
  data_known_no_elements = df1[-elements_unlist,]
  df2 = data_known_no_elements
}



#Convert formula to neutral molecules
{
  df_charge=df2
  data("isotopes")
  df_charge$MF = check_chemform(isotopes, df_charge$MF)$new_formula
  
  counter = 0
  for (i in 1:nrow(df_charge)){
    
    if(df_charge$Charge[i]!=0){
      counter = counter+1
      df_charge$MF[i] = my_calculate_formula(df_charge$MF[i], "H1", -df_charge$Charge[i])
      df_charge$Exact_Mass[i] = formula_mz(df_charge$MF[i])
    }
  }
  df3 = df_charge
}

#Remove duplicated formula entries
{
  df_unique = df3
  df_unique = df_unique[!duplicated(df_unique$MF),]
  df4 = df_unique
}

# identify invalid formula (rdbe and Nrule)
{
  df_valid = df4
  valid_log = rep(T, nrow(df_valid))
  
  i=1
  for(i in 1:nrow(df_valid)){
    temp_formula = get.formula(df_valid$MF[i])
    valid_log[i] = isvalid.formula(temp_formula, rule = c("nitrogen"))
  }
  
  df5 = df_valid[valid_log,]
  # valid_log_false = which(valid_log == F)
  # df_valid_false = df_valid[valid_log_false,]
}


write.table(df5, file = "hmdb_structure_sdf_unique_neutral_formula.tsv", sep = "\t",
            col.names = T, row.names = F, quote = F)



