library(janitor)
# install.packages("readxl")
library(readxl)

# Evaluate Xi's data from annotation ####
{
  # wl_result = read_csv("WL_190405_both.csv")
  wl_result = read_excel("WL_190405_both.xlsx")
  merge_result = merge(unknown_node_CPLEX[,c("ID","formula","is_metabolite","is_artifact", "is_biotransform","steps")],wl_result, by.x="ID", by.y = "ID", all =T)
  merge_result$formula.y = check_chemform(isotopes, merge_result$formula.y)$new_formula
}

# signal > e5
{
  Mdata = Mset$Data
  e5_id = Mdata$ID[Mdata$log10_inten>=5]
  merge_result_e5 = merge_result[merge_result$ID %in% e5_id, ]
  e6_id = Mdata$ID[Mdata$log10_inten>=6]
  merge_result_e6 = merge_result[merge_result$ID %in% e6_id, ]
  e7_id = Mdata$ID[Mdata$log10_inten>=7]
  merge_result_e7 = merge_result[merge_result$ID %in% e7_id, ]
}

# Compare formula #
{
  merge_result_with_formula = merge_result[!is.na(merge_result$`Ground truth`) & !is.na(merge_result$`Correct?`),]
  merge_result_with_formula_correct = merge_result_with_formula[merge_result_with_formula$formula.x==merge_result_with_formula$`Ground truth` &
                                                                  (!is.na(merge_result_with_formula$formula.x)),]
  merge_result_with_formula_dif = merge_result_with_formula[merge_result_with_formula$formula.x!=merge_result_with_formula$`Ground truth` |
                                                              (is.na(merge_result_with_formula$formula.x)),]
  nrow(merge_result_with_formula_correct)/nrow(merge_result_with_formula)
}


# Compare assignment #
{
  merge_result_notbg = merge_result_e5[merge_result_e5$feature!="background",]
  #merge_result_notbg = merge_result_e6[merge_result_e6$feature!="background",]
  
  PAVE_NetID= tabyl(merge_result_notbg, feature_1, is_metabolite)
  PAVE_NetID_filter = filter(merge_result_notbg, feature_1 == "Adduct", is_metabolite == "Yes")
  PAVE_NetID_filter = filter(merge_result_notbg, feature_1 == "Isotope", is_metabolite != "No")
  
  WY_NetID = tabyl(merge_result_notbg, feature, is_metabolite)
  WY_NetID_filter = filter(merge_result_notbg, feature == "buffer", is_metabolite == "Yes")
  
  fea_vs_fea1 = tabyl(wl_result, feature, feature_1 )
  artifact_vs_fea1 = tabyl(merge_result_notbg, is_artifact, feature_1,is_biotransform )
}

