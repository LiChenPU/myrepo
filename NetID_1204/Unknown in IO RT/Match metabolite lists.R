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
  
  write_csv(genome_met_valid, "genome_met_valid.csv", na="")
}


# Read NetID output files
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  work_dir = "Rt_pos"
  setwd(work_dir)
  load(list.files()[1])
  
  netid_dt = read.csv("Rt_pos.csv")
  
  HMDB_lib = Mset$Library_HMDB
  known_lib = Mset$Library_known
    
  
  result = netid_dt %>%
    mutate(Formula_in_HMDB = formula %in% HMDB_lib$formula) %>%
    mutate(Formula_in_genomodel = formula %in% genome_met_valid$formula_neutral) %>%
    mutate(Formula_in_known = formula %in% known_lib$formula)
  
  tabyl(result, class, Formula_in_HMDB)
  test = result %>%
    filter(!Formula_in_HMDB,Formula_in_known)
  
  writexl::write_xlsx(result, "Rt_pos.xlsx")

}


## find metabolites in yeast ####
{
  path = list.dirs()[-1]
  parent_dir = getwd()
  
  raw_ls = list()
  for(i in unique(path)){
    setwd(parent_dir)
    setwd(i)
    file_names = list.files()
    file_name = file_names[grepl(".xlsx", file_names)]
    
    raw_ls[[length(raw_ls)+1]] = read_xlsx(file_name) %>%
      mutate(origin = sub(".xlsx", "", file_name))
    
  }
  
  raw_df = bind_rows(raw_ls)
  # write_csv(raw_df, "all.csv", na="")
  
  target_formula = "C8H17N1O8S1"
  target_formula = "C9H17N1O6S1"
  target_mass = formula_mz(target_formula)
  proton_mass = formula_mz("H1", 1)
  
  df_filter = raw_df %>%
    # filter(abs(mass - target_mass) < target_mass * 3e-6)
    mutate(search_ppm_error = abs(abs(medMz - target_mass)-proton_mass) / target_mass * 1e6) %>%
    mutate(adjust_search_ppm_error = abs(mass - target_mass) / target_mass * 1e6) %>%
    filter(search_ppm_error < 5)
  
  
  
}




