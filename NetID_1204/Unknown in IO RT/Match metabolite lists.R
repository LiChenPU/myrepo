{
  # devtools::install_github("LiChenPU/Formula_manipulation")
  library(lc8)
  library(enviPat)
  library(dplyr)
  library(tidyr)
  # library(fitdistrplus)
  library(slam)
  # library(cplexAPI)
  library(readr)
  library(stringi)
  library(pracma)
  library(igraph)
  library(cplexAPI)
  library(readxl)
  library(stringr)
  library(janitor)
  
  
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  options(scipen=999) # Turn off scientific expression
}

# Clean genome_met formula
# Remove symbolic group formulas
# Make neutral formula
{
  genome_met = read_xlsx("iIsor_20191016.xlsx", sheet = "metabolites")
  
  data("isotopes")
  formula_check = check_chemform(isotopes, genome_met$formula)
  
  genome_met_valid = genome_met %>%
    mutate(formula = formula_check$new_formula) %>%
    filter(!formula_check$warning) %>%
    mutate(formula_neutral = apply(., 1, function(x){
      my_calculate_formula(x["formula"], "H-1", sign = as.numeric(x["charge"]))
      })
    ) 
}


# Read NetID output files
{
  work_dir = "Io_pos"
  setwd(work_dir)
  load("Io_pos.RData")
  
  netid_dt = read.csv("Rt_neg.csv")
  
  HMDB_lib = Mset$Library_HMDB
  known_lib = Mset$Library_known
    
  
  result = netid_dt %>%
    mutate(Formula_in_HMDB = formula %in% HMDB_lib$formula) %>%
    mutate(Formula_in_genomodel = formula %in% genome_met_valid$formula_neutral) %>%
    mutate(Formula_in_known = formula %in% known_lib$formula)
  
  writexl::write_xlsx(result, "Rt_neg.xlsx")
  
    

}






