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
  if(inherits(temp, "try-error")){
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
saveRDS(summary, "summary.rds")
sink()

# unique_formulas_mz = numeric(length = nrow(unique_formulas))
# 
# bad = c()
# for(i in 3441785:length(unique_formulas_mz)){
#   temp = try(formula_mz(unique_formulas_formulas[i]))
#   if(inherits(temp, "try-error")){
#     bad = c(bad, i)
#     next
#   }
#   unique_formulas_mz[i] = temp
#   if(i %% 100 == 0)print(i)
# }


## not enough memory to run ####
# unique_formulas_ls = split(unique_formulas, unique_formulas$origin)
# summary_ls = list()
# for(i in 1:length(unique_formulas_ls)){
#   temp_result = readRDS(paste0("./Unique_formula/",gsub("_unique_formula","", unique_formulas_ls[[i]]$origin[1])))
#   temp_result2 = temp_result[unique_formulas_ls[[i]]$n,]
#   summary_ls[[i]] = temp_result2
#   print(i)
# }
# 
# summary = bind_rows(summary_ls)




