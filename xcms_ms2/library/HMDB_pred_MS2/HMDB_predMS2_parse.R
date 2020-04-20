library(stringr)
library(dplyr)
library(enviPat)
library(lc8)
library(ChemmineOB)
library(ChemmineR)

data("isotopes")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

df_ls = readRDS("HMDB_pred_MS2_Data.rds")
hmdb_database = read.csv("../HMDB/hmdb_database.csv", stringsAsFactors = F)
HMDB_structure_sdf = readRDS("../HMDB/HMDB_structure_sdf.rds")

metaData = readRDS("HMDB_pred_MS2_metaData.rds")
metaData = metaData[,c(1,3,4)]
colnames(metaData) = c("HMDBID", "specID1", "specID2")
metaData$HMDBID = as.character(metaData$HMDBID)

test = metaData[!metaData$HMDBID %in% hmdb_database$accession,]

test = metaData[!metaData$HMDBID %in% HMDB_structure_sdf$HMDB_ID,]

test2 =HMDB_structure_sdf[HMDB_structure_sdf$HMDB_ID %in% hmdb_database$accession,]

table(table((metaData$HMDBID)))

# MoNA_MS2_neg = readRDS("../MoNA_MS2_Negative/MoNA_MS2_neg.rds")
# 
# sample_file = MoNA_MS2_neg[[1]]

spec_ls = list()
i=1
for(i in 1:nrow(metaData)){
  # for(i in 1:1000){
  if(i%%1000 ==0){print(i)
    print(Sys.time())}
  spectrum = df_ls[[i]]
  colnames(spectrum) = c("mz", "inten")
  spectrum$inten = spectrum$inten/max(spectrum$inten)
  
  external_id = metaData$HMDBID[i]
  
  formula = hmdb_database$chemical_formula[hmdb_database$accession==external_id]
  formula = check_chemform(isotopes, formula)$new_formula
  
  SMILES = hmdb_database$smiles[hmdb_database$accession==external_id]
  
  prcursor_mz = max(spectrum$mz)
  
  prcursor_type = NA
  
  neutral_mz = formula_mz(formula)
  polarity = ifelse(prcursor_mz > neutral_mz, "positive", "negative")
  
  if(abs(prcursor_mz - neutral_mz) > 1.008){
    prcursor_mz = NA
    polarity = "positive or negative"
  }

  
  spec_ls[[i]]= list(spectrum = spectrum,
                     formula = formula,
                     external_id = external_id,
                     SMILES = SMILES,
                     prcursor_mz = prcursor_mz,
                     prcursor_type = prcursor_type,
                     polarity = polarity)
}


massdif = sapply(spec_ls, '[[', "polarity")



test_sdf = smiles2sdf(SMILES)
plotStruc(test_sdf[[1]], regenCoords=T)


# Parsing the raw txt files
# setwd("C:/Users/lc8/Downloads/hmdb_predicted_msms_peak_lists")
# filenames = list.files(pattern = '.txt')
# df_ls = list()
# i=1
# for(i in 1:length(filenames)){
#   temp_df = read.csv(filenames[i], header = F, sep = " ")
#   df_ls[[i]]=temp_df
# }
# saveRDS(df_ls, "HMDB_pred_MS2_Data.rds")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
MS2_library = readRDS("HMDB_pred_MS2_neg.rds")

H_mass = 1.00782503224
e_mass = 0.00054857990943
ion_mode = -1

MS2_library2 = lapply(MS2_library, function(x){
  x[["spectrum"]]$mz = x[["spectrum"]]$mz + e_mass
  return(x)
})
  
all_formula = sapply(MS2_library2, "[[", "formula")

data(periodic_table)
all_elem = periodic_table$Symbol

elem_ls = list()
for(i in 1:length(all_elem)){
  query_elem = all_elem[i]
  formula_elem = str_extract(all_formula, paste0("(?<=",query_elem,")([:digit:]|\\.)+"))
  elem_ls[[query_elem]] = as.numeric(formula_elem)
}

support_elem = c("C", "H", "O", "N", "P","S", "Cl",
                 "Na", "K", "F", "Br", "I", "Si", "B",
                 "Ca", "Cu", "Ni")
elem_df = bind_rows(elem_ls) %>%
  dplyr::select(-support_elem) %>%
  mutate(unsupport_elem = rowSums(., na.rm=T))

all_formula_valid = !grepl("\\(", all_formula) & (elem_df$unsupport_elem == 0)
MS2_library3 = list()
for(i in 1:length(MS2_library2)){
# MS2_library3 = lapply(MS2_library2, function(test){
  if(!all_formula_valid[i]){
    MS2_library3[[i]] = MS2_library2[[i]]
    next
  }
  test = MS2_library2[[i]]
  test_spec = test$spectrum
  test_spec = test_spec %>%
    # Note here: my_pred_formula will make mistake when having too many heavy elements
    mutate(formula = my_pred_formula(test$spectrum$mz, ion_mode, parent_formula = test$formula, N_rule = F, ppm = 15)) %>%
    mutate(cal.mz = formula_mz(formula) + (H_mass-e_mass)*ion_mode,
           mz.dif = cal.mz - mz) %>%
    # mutate(toremove = abs(mz.dif) < 1e-4)
    # mutate(formula = ifelse(abs(mz.dif) < 1e-4, formula, NA)) %>%
    filter(abs(mz.dif) < 1e-4) %>%
    
    dplyr::select(mz, inten, formula)
  
  if(nrow(test_spec) < nrow(test$spectrum)){warning(i)}
  
  # temp[i] = max(abs(test_spec$mz.dif))
  test$spectrum = test_spec
  
  MS2_library3[[i]] = test
  print(i)
  # return(test)
}
saveRDS(MS2_library3, "HMDB_pred_MS2_neg2.rds")


## predict formula ####
my_pred_formula=function(mz = df$mz, ion_mode,
                         parent_formula = "C99H100N15O15S3P3", N_rule = F, 
                         ppm=15, db_max=40){
  # parent_formula = "C54H56N2O21"
  
  support_elem = c("C", "H", "O", "N", "Cl", "P","S",
                   "Na", "K", "F", "Br", "I", "Si", "B",
                   "Ca", "Cu", "Ni")
  elem_num = 0
  for(i in 1:length(support_elem)){
    query_elem = support_elem[i]
    formula_elem = str_extract(parent_formula, paste0("(?<=",query_elem,")([:digit:]|\\.)+"))
    elem_num[i] = as.numeric(formula_elem)
  }
  elem_num[is.na(elem_num)]=0
 
  predict_formula = rep("", length(mz))
  # i=15
  for(i in 1:length(mz)){
    temp = mz_formula(mz[i], charge = ion_mode,  ppm=max(ppm, 1000/mz[i]), N_rule = N_rule, db_max = db_max, db_min = -3,
                      C_range = 0:elem_num[1], 
                      H_range = 0:elem_num[2], 
                      O_range = 0:elem_num[3], 
                      N_range = 0:elem_num[4], 
                      Cl_range = 0:elem_num[5], 
                      P_range = 0:elem_num[6], 
                      S_range = 0:elem_num[7], 
                      Na_range = 0:elem_num[8], 
                      K_range = 0:elem_num[9], 
                      F_range = 0:elem_num[10], 
                      Br_range = 0:elem_num[11], 
                      I_range = 0:elem_num[12], 
                      Si_range = 0:elem_num[13], 
                      B_range = 0:elem_num[14], 
                      Ca_range = 0:elem_num[15], 
                      Cu_range = 0:elem_num[16], 
                      Ni_range = 0:elem_num[17])
    
    if(!is.data.frame(temp)){
      predict_formula[i]=round(mz[i],digits = 3)
      next
    }
    temp["abs"]=abs(temp$differ)
    temp = temp[with(temp, order(abs)),]
    predict_formula[i]=temp$formula[1]
  }
  return(predict_formula)
}













