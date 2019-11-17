# BiocManager::install("ChemmineR")
library(ChemmineR)
library(lc8)
library(dplyr)
library(enviPat)

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rdsfilenames = list.files(pattern = "\\.rds")

sink("log.txt")
timestamp()
time = Sys.time()

for(repeating in 0:1)
{
  output_name = substr(rdsfilenames[40*repeating+ 1], 1, 12)
  print(paste("Processing", output_name))
  if(file.exists(paste0(output_name,".rds"))){next}
  rds_ls = list()
  # for(i in 1:length(rdsfilenames)){
  for(i in 1:40){
    rds_ls[[i]] = readRDS(rdsfilenames[40*repeating+i])
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
    # data(periodic_table)
    # data("elem_table")
    # all_elem = periodic_table$Symbol
    # keep_elem = all_elem[c(1:36, 53)]
    # filter_elem = all_elem[-c(1:36, 53)]
    # filter_elem_char = paste(filter_elem, collapse = "|")
    
    f1 = f0 %>%
      filter(!grepl(filter_elem_char, PUBCHEM_MOLECULAR_FORMULA)) %>%
      # filter(PUBCHEM_ISOTOPIC_ATOM_COUNT == 0) %>%
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
    
    
    f2 = f1 %>%
      filter(PUBCHEM_TOTAL_CHARGE == 0) %>%
      mutate(neutral_formula = PUBCHEM_MOLECULAR_FORMULA) %>%
      mutate(neutral_formula = check_chemform(isotopes, PUBCHEM_MOLECULAR_FORMULA)$new_formula) %>%
      bind_rows(f2_charge) %>%
      distinct(neutral_formula, .keep_all = T)# %>%
    # mutate(exact_mz = sapply(neutral_formula, formula_mz))
    # 
    # test_f2 = f2 %>%
    #   mutate(mz_dif = exact_mz - PUBCHEM_MONOISOTOPIC_WEIGHT) %>%
    #   filter(PUBCHEM_ISOTOPIC_ATOM_COUNT == 0, PUBCHEM_TOTAL_CHARGE == 0)
  }
  
  saveRDS(f2, file = paste0("./Unique_formula/", output_name, ".rds"))
  
  unique_formula = f2 %>%
    mutate(n = 1:nrow(.)) %>%
    dplyr::select(n, PUBCHEM_COMPOUND_CID, neutral_formula)
  saveRDS(unique_formula, file = paste0("./Unique_formula/", output_name, "_unique_formula.rds"))
  
  print(Sys.time()- time)
}

sink()




