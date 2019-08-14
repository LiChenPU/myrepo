library(enviPat)
library(jsonlite)
library(stringr)
library(dplyr)

# setwd("C:/Users/lc8/Downloads/MoNA-export-LC-MS-MS_Negative_Mode-json")
setwd(work_dir)
print(getwd())
# json_file = "MoNA-export-LC-MS-MS_Negative_Mode.json"
# jsonlite_data <- jsonlite::fromJSON(json_file, simplifyDataFrame = T, simplifyVector = F, simplifyMatrix = T)
# saveRDS(jsonlite_data, "jsonlite_data.rds")
# setwd("C:/Users/lc8/Documents/GitHub/myrepo/xcms_ms2/library/MoNA_MS2_Positive")
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
  spec_df[2,] = spec_df[2,]/max(spec_df[2,])
  spec_df = t(spec_df)
  colnames(spec_df) = c("mz", "inten")
  spec_df = as.data.frame(spec_df)
  spec_df = spec_df %>%
    arrange(desc(inten)) %>%
    top_n(100, inten)
  spec_df = as.matrix(spec_df)
  
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





  