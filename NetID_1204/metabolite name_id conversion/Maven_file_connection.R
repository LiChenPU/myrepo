library(lc8)
library(enviPat)
library(readr)
library(readxl)
library(ChemmineR)
library(stringr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# filename = "./dependent/Knowns-new-Hilic25min-QE.csv"
# 
# raw = read_csv(filename)
# raw = read_csv("raw.csv")
raw = read_xlsx("ID_Metabolites_xi.xlsx")

data(isotopes)

# check curated formula vs mass
check = check_chemform(isotopes, raw$Formula)
f1 = raw %>%
  mutate(formula = check$new_formula,
         exact_mass = check$monoisotopic_mass)
# write.csv(f1, "f1.csv")

# check SMILES formula vs curated formula
formulas = character(nrow(f1))
for(i in 1:nrow(f1)){
  # print(i)
  test = smiles2sdf(f1$SMILES[i])
  formulas[i] = MF(test, addH=T)
}
# 
formulas_check = check_chemform(isotopes, formulas)$new_formula

# count charge
SMILES_charges_ls = str_match_all(f1$SMILES, "\\+\\d|\\-\\d|\\+|\\-")

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

f2 = f1 %>%
  mutate(charge = SMILES_charges) %>%
  mutate(neutral_formula = mapply(my_calculate_formula, formulas_check, "H-1", charge)) %>%
  mutate(same_smiles = neutral_formula == formula) %>%
  mutate(mass_dif = exact_mass - mz)

# calculate rdbe 
f3 = f2 %>%
  mutate(rdbe = sapply(neutral_formula, formula_rdbe)) 

write.csv(f3, "f3.csv")

f4 = f3 %>%
  group_by(HMDB) %>%
  filter(n()>1)



# check for redundant entry
# f3 = f2 %>% 
#   group_by(HMDB) %>%
#   filter(n()==1)
# 
# f3_c = read_csv("./dependent/f3.csv")
# f4 = bind_rows(f3, f3_c)
# write.csv(f4, "f4.csv")

# raw = raw %>%
#   arrange(exact_mass, rt) %>%
#   mutate(id = 1:nrow(.))
# duplicate_formula = raw$formula[grepl("\\(\\d\\)", raw$name)]
# f5 = raw %>%
#   group_by(Match) %>%
#   filter(n()>1)

# Check for decimal rdbe (indicating incorrect formula)
# raw = raw %>%
#   filter(grepl("\\.", rdbe))

# write.csv(raw, "raw.csv")







