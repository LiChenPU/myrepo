# BiocManager::install("ChemmineR")
library(ChemmineR)
library(lc8)
library(dplyr)
library(enviPat)

sink("log.txt")
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rdsfilenames = list.files(path = "./Unique_formula",pattern = "unique_formula\\.rds")

print(timestamp())

time = Sys.time()
rds_ls = list()
for(i in 1:length(rdsfilenames))
{
  temp_result = readRDS(paste0("./Unique_formula/",rdsfilenames[i]))
  rds_ls[[i]] = temp_result %>%
    mutate(origin = rdsfilenames[i])
  # print(i)
}

unique_formulas = bind_rows(rds_ls) %>% 
  distinct(neutral_formula, .keep_all = T)
rm("rds_ls")

unique_formulas_ls = split(unique_formulas, unique_formulas$origin)
data(isotopes)
check_formula_ls = list()
bad = c()
for(i in 1:length(unique_formulas_ls)){
  print(paste("prosessing",i))
  print(Sys.time()-time)
  check_formula = try(check_chemform(isotopes, unique_formulas_ls[[i]]$neutral_formula))
  if(inherits(check_formula, "try-error")){
    bad = c(bad, i)
    next
  }
  check_formula_ls[[i]] = unique_formulas_ls[[i]] %>%
    mutate(invalid_formula = check_formula$warning,
           monoisotopic_mass = check_formula$monoisotopic_mass)
  
}

summary = bind_rows(check_formula_ls) %>%
  filter(!invalid_formula) %>%
  dplyr::select(-invalid_formula)

print(paste0("summary nrow=", nrow(summary)))
# saveRDS(summary, "summary.rds")
sink()


summary = readRDS("summary.rds")

test = summary %>%
  filter(invalid_formula)

## Filter elements ####
{
  data(periodic_table)
  data("elem_table")
  all_elem = periodic_table$Symbol
  keep_elem = elem_table$element 
  keep_elem = keep_elem[keep_elem %in% all_elem & keep_elem!="Dy"]
  filter_elem = all_elem[!all_elem %in% keep_elem]
  filter_elem_char = paste(c(filter_elem), collapse = "|")
  
  
  summary_filter = summary %>%
    filter(!grepl(filter_elem_char, neutral_formula))
}

test = summary_filter %>%
  filter(grepl("-", neutral_formula))
  
  
  
  

