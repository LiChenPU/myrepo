library(lc8)
library(enviPat)
library(readr)
library(ChemmineR)
library(stringr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

hmdb_database = read_csv("hmdb_raw.csv")

time = Sys.time()

# Remove invalid formula entry ####
{
  data(isotopes)
  formula_check = check_chemform(isotopes, hmdb_database$chemical_formula)
  formula_check_warning = formula_check %>%
    filter(!warning)
  
  hmdb_database = hmdb_database %>%
    filter(!grepl("\\(", chemical_formula)) %>%
    mutate(chemical_formula = formula_check_warning$new_formula)

  
  hmdb_database = hmdb_database %>%
    filter(monisotopic_molecular_weight > 3)
}

# Count charges from smiles ####
{
  SMILES_error_entry = hmdb_database %>%
    filter(!str_detect(smiles, "\\+\\d]|\\-\\d]")) %>%
    filter(str_detect(smiles, "\\+\\d|\\-\\d"))
  
  SMILES_charges_ls = str_match_all(hmdb_database$smiles, "\\+\\d|\\-\\d|\\+|\\-")
  
  SMILES_charges = rep(0,length(SMILES_charges_ls))
  for(i in 1:length(SMILES_charges_ls)){
    if(dim(SMILES_charges_ls[[i]])[1] == 0){
      next
    }
    temp = unlist(SMILES_charges_ls[i])
    pos = sum(temp == "+")
    neg = sum(temp == "-")
    if(pos + neg == length(temp)){
      SMILES_charges[i] = pos - neg
      next
    }
    pos = pos + sum(as.numeric(str_extract(temp, "\\+\\d")), na.rm = T)
    neg = -neg + sum(as.numeric(str_extract(temp, "\\-\\d")), na.rm = T)
    SMILES_charges[i] = pos + neg
  }
  
  hmdb_database_charge = hmdb_database %>%
    mutate(charge = SMILES_charges) %>%
    filter(charge != 0) %>%
    mutate(formula = mapply(my_calculate_formula, chemical_formula, "H-1", charge, F))
  
  hmdb_database_neutral = hmdb_database_charge %>%
    bind_rows(hmdb_database %>% 
                mutate(charge = 0,
                       formula = chemical_formula)
              ) %>%
    distinct(accession, .keep_all=T)
}

# filter elements ####
{
  data(periodic_table)
  all_elem = periodic_table$Symbol
  keep_elem = all_elem[c(1:36, 53)]
  filter_elem = all_elem[-c(1:36, 53)]
  filter_elem_char = paste(filter_elem, collapse = "|")
  
  hmdb_database_element = hmdb_database_neutral %>%
    filter(!grepl(filter_elem_char, formula))
}

# calculate mass ####
{
  # mass = rep(0, nrow(hmdb_database_element))
  # for(i in 1:nrow(hmdb_database_element)){
  #   print(i)
  #   mass[i] = formula_mz(hmdb_database_element$formula[i])
  # }
  hmdb_database_mass = hmdb_database_element %>%
    mutate(mass = sapply(formula, formula_mz)) 

}

# arrangment ####
{
  hmdb_database_arrange = hmdb_database_mass %>%
    arrange(accession)
}

write.csv(hmdb_database_arrange, "hmdb_database.csv", row.names = F)

print(Sys.time()-time)