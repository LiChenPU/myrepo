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


