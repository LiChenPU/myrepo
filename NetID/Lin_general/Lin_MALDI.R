
raw_data = Mset$Raw_data 
CPLEX_result = formula_list %>%
  filter(ILP_result != 0) %>%
  mutate(msr_cal_mass_dif_ppm = round(msr_cal_mass_dif_ppm, 2)) %>%
  dplyr::select(ID, formula, is_metabolite, ILP_result, parent, msr_cal_mass_dif_ppm)# %>%
  # mutate(formula = gsub("([[:alpha:]])1([[:alpha:]]|$)", "\\1\\2", formula))
result = merge(CPLEX_result, raw_data, by.x = "ID", by.y = "groupId", all = T)

write.csv(result, "NetID.csv", row.names = F, na = "")
