


library(enviPat)
library(jsonlite)
library(stringr)


# setwd("C:/Users/lc8/Downloads/MoNA-export-LC-MS-MS_Positive_Mode-json")
print(work_dir)
setwd(work_dir)
json_file = "MoNA-export-LC-MS-MS_Positive_Mode.json"
# jsonlite_data <- jsonlite::fromJSON(json_file, simplifyDataFrame = T, simplifyVector = F, simplifyMatrix = T)
jsonlite_data = readRDS("jsonlite_data.rds")

spec_ls = list()
i=1
data("isotopes")

for(i in 1:nrow(jsonlite_data)){
  # for(i in 1:19){
  if(i%%1000==0){
    print(i)
    print(Sys.time())
  }
  temp_json = jsonlite_data[i,]
  
  # mz, intensity
  spectrum = temp_json$spectrum
  spectrum_spl = str_split(spectrum, " ")[[1]]
  spectrum_spl2 = str_split(spectrum_spl, ":")
  spec_df = matrix(as.numeric(unlist(spectrum_spl2)), nrow=2)
  spec_df[2,] = spec_df[2,]/100
  spec_df = t(spec_df)
  colnames(spec_df) = c("mz", "inten")
  
  
  temp_json_compound_metaData = temp_json$compound[[1]]$metaData[[1]]
  if(length(temp_json_compound_metaData)!=0){
    formula = temp_json_compound_metaData$value[temp_json_compound_metaData$name=="molecular formula"]
    if(length(formula)!=0){
      formula = check_chemform(isotopes, chemforms = formula)$new_formula
    } else {
      formula =""
    }
    
    # ID, HMDBID or anything that can be searchable
    external_id = temp_json_compound_metaData[temp_json_compound_metaData$category=="external id",c("name", "value")]
    # structure
    SMILES = temp_json_compound_metaData$value[temp_json_compound_metaData$name=="SMILES" & temp_json_compound_metaData$computed]
    # compound_class
    compound_class = temp_json_compound_metaData$value[temp_json_compound_metaData$name=="compound class"]
  } else{
    external_id = NA
    SMILES = NA
    compound_class = NA
  }
  # formula
  
  temp_json_metaData = temp_json$metaData[[1]]
  # precusor mz
  prcursor_mz = as.numeric(temp_json_metaData$value[temp_json_metaData$name=="precursor m/z"])
  # ion state
  prcursor_type = temp_json_metaData$value[temp_json_metaData$name=="precursor type"]
  # polarity
  polarity = temp_json_metaData$value[temp_json_metaData$name=="ionization mode"]
  
  
  
  spec_ls[[i]]= list(spectrum = spec_df,
                     formula = formula,
                     external_id = external_id,
                     SMILES = SMILES,
                     prcursor_mz = prcursor_mz,
                     prcursor_type = prcursor_type,
                     polarity = polarity)
}





saveRDS(spec_ls, "MoNA_MS2_pos.rds")

spec_ls = readRDS("MoNA_MS2_pos.rds")





  