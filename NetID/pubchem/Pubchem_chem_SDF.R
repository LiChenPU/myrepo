# BiocManager::install("ChemmineR")
library(ChemmineR)
library(lc8)
library(dplyr)
library(enviPat)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

filenames = list.files(pattern = ".sdf.gz")
for(i in 1:10){
  if(file.exists(sub("sdf.gz","rds", filenames[i]))){next} # skip when rds file exists
  tictoc::tic()
  dsfset = read.SDFset(filenames[i])
  valid <- validSDF(dsfset) # Identifies invalid SDFs in SDFset objects 
  dsfset <- dsfset[valid] # Removes invalid SDFs, if there are any 
  
  temp_result_ls = datablock(dsfset)
  temp_result = rbind_list(temp_result_ls)
  
  saveRDS(temp_result, sub("sdf.gz","rds", filenames[i]))
  tictoc::toc()
}

rdsfilenames = list.files(pattern = "\\.rds")


rds_ls = list()
for(i in 1:length(rdsfilenames)){
  rds_ls[[i]] = readRDS(rdsfilenames[i])
}

summary = bind_rows(rds_ls)

## format ####
{
  numericcharacters <- function(x) {
    !any(is.na(suppressWarnings(as.numeric(x)))) & is.character(x)
  }
  f0 = summary %>%
    mutate_if(numericcharacters, as.numeric)
}



## Filter elements ####
{
  data(periodic_table)
  data("elem_table")
  all_elem = periodic_table$Symbol
  keep_elem = all_elem[c(1:36, 53)]
  filter_elem = all_elem[-c(1:36, 53)]
  filter_elem_char = paste(filter_elem, collapse = "|")
  
  f1 = f0 %>%
    filter(!grepl(filter_elem_char, PUBCHEM_MOLECULAR_FORMULA)) %>%
    distinct(PUBCHEM_MOLECULAR_FORMULA, .keep_all = T)
}

## Neutral formula ####
{
  data(isotopes)

  f2_charge = f1 %>%
    mutate_if(numericcharacters, as.numeric) %>%
    filter(PUBCHEM_TOTAL_CHARGE != 0) %>%
    mutate(neutral_formula = check_chemform(isotopes, gsub("\\+.*|-.*", "", PUBCHEM_MOLECULAR_FORMULA))$new_formula)
    
  
  temp_neutral_formula = character(nrow(f2_charge))
  for(i in 1:nrow(f2_charge)){
    temp_neutral_formula[i] = my_calculate_formula(f2_charge$neutral_formula[i], 
                                              paste0("H",abs(f2_charge$PUBCHEM_TOTAL_CHARGE[i])),
                                              -1 * sign(f2_charge$PUBCHEM_TOTAL_CHARGE[i]),
                                              Is_valid = F)
  }
  
  f2_charge = f2_charge %>% 
    mutate(neutral_formula = temp_neutral_formula)
  
  tictoc::tic()
  f2 = f1 %>%
    filter(PUBCHEM_TOTAL_CHARGE == 0) %>%
    # mutate(neutral_formula = PUBCHEM_MOLECULAR_FORMULA) %>%
    mutate(neutral_formula = check_chemform(isotopes, PUBCHEM_MOLECULAR_FORMULA)$new_formula) %>%
    bind_rows(f2_charge) %>%
    distinct(neutral_formula, .keep_all = T) %>%
    mutate(exact_mz = sapply(neutral_formula, formula_mz))
  tictoc::toc()
  
  test_f2 = f2 %>%
    mutate(mz_dif = exact_mz - PUBCHEM_MONOISOTOPIC_WEIGHT) %>%
    filter(PUBCHEM_ISOTOPIC_ATOM_COUNT == 0, PUBCHEM_TOTAL_CHARGE == 0)
}


write.table(hmdb_summary, file = "HMDB_structure_sdf.tsv", sep = "\t",
            col.names = T, row.names = F, quote = F)

# 
# 
# 
# 
# setwd("C:/Users/lc8/Desktop/hmdb_metabolites")
# 
# library(XML)
# input <- "hmdb_metabolites.xml"
# items <- NULL
# maxItems <- 50
# 
# parseItem = function (parser, node, ...) {
#   children <- xmlChildren(node)
#   items <<- rbind(items, sapply(children, xmlValue))
#   if (nrow(items) == maxItems) {
#     xmlStopParser(parser)
#   }
# }
# 
# # with XMLParserContextFunction, we get the parser as first parameter
# # so we can call xmlStopParser
# class(parseItem) = c("XMLParserContextFunction", "SAXBranchFunction")
# 
# xmlEventParse(input,
#               branches = list(item = parseItem),
#               ignoreBlanks = T
# )
# 
# items <- as.data.frame(items)



