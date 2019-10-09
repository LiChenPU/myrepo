library(janitor)
# install.packages("readxl")
library(readxl)
library(readr)
library(dplyr)
library(enviPat)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# Read previous CPLEX results
{
  NetID_files = list.files(pattern = "[0-9].csv")
  formula_list_ls = list()
  for(i in 1:length(NetID_files)){
    formula_list_ls[[length(formula_list_ls)+1]] = read.csv( NetID_files[i], stringsAsFactors = F) 
  }
  
  # formula_list2 = formula_list_ls[[6]]
}

evaluate_summary = list()
for(i in 1:length(NetID_files)){
  formula_list2 = formula_list_ls[[i]]
  


# Evaluate Xi's data from annotation ####
{
  # wl_result = read_csv("WL_190405_both.csv")
  
  
  data("isotopes")
  wl_result = read_excel("WL_190405_both_190924.xlsx") %>%
    mutate(formula = check_chemform(isotopes, formula)$new_formula)
  merge_result = formula_list2 %>%
    dplyr::select(ID,formula,is_metabolite, is_artifact, is_biotransform, ILP_result) %>%
    merge(wl_result, by.x="ID", by.y = "id", all.y = T)
}

# Evaluate all metabolites ####
{
  all_PAVE_metabolites = merge_result %>%
    filter(feature...11 == "Metabolite")
  all_PAVE_metabolites_assigned = all_PAVE_metabolites %>%
    arrange(-ILP_result) %>%
    distinct(ID, .keep_all = T)
    
  all_PAVE_metabolites_correct = all_PAVE_metabolites %>%
    filter(formula.x == formula.y) %>%
    filter(ILP_result!=0)
  all_PAVE_metabolites_wrong = all_PAVE_metabolites %>%
    filter(!ID %in% all_PAVE_metabolites_correct$ID) #%>%
    # filter(ILP_result!=0)
}

# Evaluate if unconnected in PAVE background ####
{
  ID_in_CPLEX = merge_result %>% filter(ILP_result!=0) %>% pull(ID)
  all_unconnected = merge_result %>%
    filter(!ID %in% ID_in_CPLEX) %>%
    distinct(ID, .keep_all = T) # %>% pull(feature...11) %>% table()
  
  # all_unconnected2 = merge_result %>%
  #   filter(!ID %in% ID_in_CPLEX) %>%
  #   distinct(ID, .keep_all = T)  %>% pull(feature...11) %>% table()
  # test = merge(data.frame(all_unconnected2),data.frame(all_unconnected), by = ".", all =T) %>%
  #   mutate(Freq.y = replace_na(Freq.y, 0)) %>% 
  #   mutate(percentage = Freq.y/Freq.x)
  # sum(test$Freq.y)/sum(test$Freq.x)
  
  
  all_unconnected_non_backgrounds = all_unconnected %>%
    filter(feature...11 != "Background")
  
  all_unconnected_non_backgrounds_3e5 = all_unconnected_non_backgrounds %>%
    filter(log10_inten>log10(3e5))
}


# Evaluate all isotope-labeled peaks
{
  all_labeled_peaks = merge_result %>%
    filter(feature...11 != "Background") %>%
    arrange(-ILP_result) 
  
  all_labeled_peaks_3e5 = all_labeled_peaks %>%
    filter(log10_inten > log10(3e5))
  
  all_labeled_peaks_3e5_distinct = all_labeled_peaks_3e5 %>%
    distinct(ID, .keep_all=T)
  
  all_labeled_peaks_3e5_correct_potential = all_labeled_peaks_3e5 %>%
    filter(formula.x == `Ground truth`)
  
  all_labeled_peaks_3e5_correct = all_labeled_peaks_3e5_correct_potential %>%
    filter(ILP_result!=0)

  all_labeled_peaks_3e5_wrong = all_labeled_peaks_3e5_distinct %>%
    filter(!ID %in% all_labeled_peaks_3e5_correct$ID)
}

  
  evaluate_summary[[i]] = list(Total_potential_3e5 = all_labeled_peaks_3e5,
                               Correct_potential = all_labeled_peaks_3e5_correct_potential,
                               Correct = all_labeled_peaks_3e5_correct,
                               Wrong = all_labeled_peaks_3e5_wrong,
                               # Unconnected = all_unconnected_non_backgrounds,
                               Unconnected_3e5 = all_unconnected_non_backgrounds_3e5
                               )
}
  
for(i in 1:length(NetID_files)){
  results = sapply(evaluate_summary[[i]], nrow)
  print(results)
}

selected_file = 12
test1 = evaluate_summary[[selected_file]]$Correct_potential %>%
  filter(!ID %in% evaluate_summary[[selected_file]]$Correct$ID) %>% 
  dplyr::select(ID) %>%
  inner_join(evaluate_summary[[selected_file]]$Total_potential_3e5)
result12 = evaluate_summary[[12]]$Total_potential_3e5 %>%
  dplyr::select(ID, ILP_result, formula.x)

result13 = evaluate_summary[[13]]$Total_potential_3e5 %>%
  dplyr::select(ID, ILP_result, formula.x)

result13_not12 = evaluate_summary[[13]]$Total_potential_3e5 %>%
  anti_join(result12) # %>%
  filter(formula.x == `Ground truth`, ILP_result ==1)

result12_not13 = evaluate_summary[[12]]$Total_potential_3e5 %>%
  anti_join(result13) %>%
  filter(formula.x == `Ground truth`, ILP_result ==1)





# Evaluate others ####
# Evaluate isotopes 
{
  all_PAVE_isotopes = merge_result %>% # dplyr::select(feature...11) %>% table
    filter(feature...11 == "Isotope")
  all_PAVE_isotopes_correct = all_PAVE_isotopes %>%
    filter(grepl("\\[",formula.x )) %>%
    filter(ILP_result!=0)
  all_PAVE_isotopes_wrong = all_PAVE_isotopes %>%
    filter(!grepl("\\[",formula.x )) %>%
    filter(ILP_result!=0)
}
# Evaluate Adducts 
{
  all_PAVE_adduct = merge_result %>%  #dplyr::select(feature...11) %>% table
    filter(feature...11 == "Adduct")
  all_PAVE_adduct_correct = all_PAVE_adduct %>%
    filter(is_artifact) %>%
    filter(ILP_result!=0)
  all_PAVE_adduct_wrong = all_PAVE_adduct %>%
    filter(!is_artifact) %>%
    filter(ILP_result!=0)
}
# Evaluate unknowns
{
  all_PAVE_unknowns = merge_result %>%
    filter(feature...11 == "[]") %>%
    mutate(proposed_C = sapply(formula.x, elem_num_query, "C"),
           proposed_N = sapply(formula.x, elem_num_query, "N")) %>%
    filter(ILP_result !=0 )
  
  sum(is.na(all_PAVE_unknowns$formula.x))
  
  all_PAVE_unknowns_CNmatch = all_PAVE_unknowns %>%
    filter(C==proposed_C, N==proposed_N)
    
  all_PAVE_unknowns_Cmatch = all_PAVE_unknowns %>%
    filter(C==proposed_C, N!=proposed_N)
  
  all_PAVE_unknowns_CN_2match = all_PAVE_unknowns %>%
    filter(abs(C - proposed_C)<=2, abs(N-proposed_N)<=2)
  
}

# Evaluate step info 
{
  CPLEX_step = formula_list2 %>% 
    filter(ILP_result == 1) %>%
    filter(category ==1) %>%  # measured data
    filter(steps >= 6) %>%
    pull(steps)
  
  table(CPLEX_step)
  
}

# Compare with simply select top metabolites ####
# {
#   all_top_metabolites = bind_rows(CPLEXset$data$pred_formula_ls) %>%
#     distinct(id, .keep_all = T)
#   merge_result2 = all_top_metabolites %>%
#     merge(wl_result, by.x="id", by.y = "id", all.y = T)
# }
# Compare assignment ####
{
  Mdata = Mset$Data
  e5_id = Mdata$ID[Mdata$log10_inten>=5]
  merge_result_e5 = merge_result[merge_result$ID %in% e5_id, ]
  e6_id = Mdata$ID[Mdata$log10_inten>=6]
  merge_result_e6 = merge_result[merge_result$ID %in% e6_id, ]
  e7_id = Mdata$ID[Mdata$log10_inten>=7]
  merge_result_e7 = merge_result[merge_result$ID %in% e7_id, ]
  merge_result_notbg = merge_result_e5[merge_result_e5$feature...11!="Background",]
  #merge_result_notbg = merge_result_e6[merge_result_e6$feature!="background",]
  
  PAVE_NetID= tabyl(merge_result_notbg, feature...11, is_metabolite)
  PAVE_NetID_filter = filter(merge_result_notbg, feature_1 == "Adduct", is_metabolite == "Yes")
  PAVE_NetID_filter = filter(merge_result_notbg, feature_1 == "Isotope", is_metabolite != "No")
  
  WY_NetID = tabyl(merge_result_notbg, feature, is_metabolite)
  WY_NetID_filter = filter(merge_result_notbg, feature == "buffer", is_metabolite == "Yes")
  
  fea_vs_fea1 = tabyl(wl_result, feature, feature_1 )
  artifact_vs_fea1 = tabyl(merge_result_notbg, is_artifact, feature_1,is_biotransform )
}



# 
# 
# # Compare formula #
# {
#   merge_result_with_formula = merge_result[!is.na(merge_result$`Ground truth`) & !is.na(merge_result$`Correct?`),]
#   merge_result_with_formula_correct = merge_result_with_formula[merge_result_with_formula$formula.x==merge_result_with_formula$`Ground truth` &
#                                                                   (!is.na(merge_result_with_formula$formula.x)),]
#   merge_result_with_formula_dif = merge_result_with_formula[merge_result_with_formula$formula.x!=merge_result_with_formula$`Ground truth` |
#                                                               (is.na(merge_result_with_formula$formula.x)),]
#   nrow(merge_result_with_formula_correct)/nrow(merge_result_with_formula)
# }




