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
  NetID_edge_files = list.files(pattern = "edge.csv")
  formula_list_ls = edge_info_sum_ls = list()
  for(i in 1:length(NetID_files)){
    formula_list_ls[[length(formula_list_ls)+1]] = read.csv(NetID_files[i], stringsAsFactors = F) 
    edge_info_sum_ls[[length(edge_info_sum_ls)+1]] = read.csv(NetID_edge_files[i], stringsAsFactors = F)
  }
  
  # formula_list2 = formula_list_ls[[6]]
}

# Evaluate Xi's data from annotation ####

evaluate_summary = list()
# for(i in 1:length(NetID_files)){
  formula_list2 = formula_list_ls[[i]]
  relation_list2 = edge_info_sum_ls[[i]]
  
{
  # wl_result = read_csv("WL_190405_both.csv")
  data("isotopes")
  wl_result = read_excel("WL_190405_1021category.xlsx") %>%
    mutate(formula = check_chemform(isotopes, formula)$new_formula)
  merge_result = formula_list2 %>%
    dplyr::select(ID,formula,is_metabolite, Artifact_assignment, Biotransform, ILP_result) %>%
    merge(wl_result, by.x="ID", by.y = "id", all.y = T) 
}

## 552 evaluation

{
  # All ILP!=0
  {
    all_filter = merge_result %>%
      filter(feature...11 != "Background") %>%
      filter(log10_inten>log10(3e5)) %>%
      filter(!category_manual %in% c("NA", "No_label"))
    
    all_ID_filter = all_filter %>%
      arrange(-ILP_result) %>%
      distinct(ID, .keep_all = T) %>%
      arrange(ID)# %>%
    # mutate(pred_C = sapply(formula.x, elem_num_query, "C"),
    #        pred_N = sapply(formula.x, elem_num_query, "N")) %>%
    # mutate(category = case_when(
    #   feature...11 == "Metabolite" ~ "Metabolite",
    #   feature...11 == "Fragment" ~ "Fragment",
    #   feature...11 == "Background" ~ "",
    #   is_metabolite == "No" ~ "Artifact",
    #   feature...11 %in% c("Dimer", "Heterodimer") & grepl("Oligomer|Heterodimer", Artifact_assignment) ~ "Artifact",
    #   feature...11 %in% c("Low_score", "[]") & pred_C== C & pred_N== N ~ "Metabolite_CNmatch", 
    #   is_metabolite == "Yes"  ~ "Maybe",
    #   is_metabolite == "Maybe" ~ "Maybe"
    # ))
    
  
    table(all_ID_filter$category_manual)
    
    
    test = all_ID_filter %>%
      # filter(feature...11 == "Odd_N")# %>%
      filter(category_manual == "Metabolite_CNnotmatch")
    tabyl(all_ID_filter, feature...11, category_manual)
    tabyl(all_ID_filter, is_metabolite, category_manual)
    
    # write.csv(all_ID,"all_ID.csv", row.names = F)
  }
  
  ## Matching to HMDB
  {
    HMDB_library = read.csv("../dependent/HMDB_CHNOPS_clean.csv")
    all_ID_filter_HMDB = all_ID_filter %>%
      filter(grepl("Metabolite", category_manual)) %>%
      mutate(In_HMDB = `Ground truth` %in% HMDB_library$MF)
    
    all_ID_filter_HMDB_false = all_ID_filter_HMDB %>%
      filter(!In_HMDB)
    all_ID_filter_HMDB_true = all_ID_filter_HMDB %>%
      filter(In_HMDB)
    # PAVE performance
    table(all_ID_filter_HMDB_true$feature...11)
    table(all_ID_filter_HMDB_false$feature...11)
    # NetID performance
    all_ID_filter_HMDB_NetID_true = all_ID_filter_HMDB_true %>%
      filter(ILP_result == 1) %>%
      filter(formula.x == `Ground truth`)
    table(all_ID_filter_HMDB_NetID_true $is_metabolite)
    table(all_ID_filter_HMDB_true$category_manual)
    
    
      
    table(all_ID_filter_HMDB$In_HMDB)
    tabyl(all_ID_filter_HMDB, In_HMDB, category_manual)
  }
  
  ## Artifacts 
  {
    all_ID_filter_artifact = all_ID_filter %>%
      filter(!grepl("Metabolite", category_manual)) 
    
    # PAVE performance
    table(all_ID_filter_artifact$feature...11)
    
    # NetID performance formula
    all_ID_filter_artifact_NetID_true = all_ID_filter_artifact %>%
      filter(ILP_result == 1) %>%
      filter(formula.x == `Ground truth`) 
    all_ID_filter_artifact_NetID_true_Yes = all_ID_filter_artifact_NetID_true %>%
      filter(is_metabolite == "Yes")
    table(all_ID_filter_artifact_NetID_true$is_metabolite)
    
    all_ID_filter_artifact_NetID_false = all_ID_filter_artifact %>%
      filter(ILP_result != 1 | formula.x != `Ground truth` | is.na(formula.x)) %>%
      rbind(all_ID_filter_artifact_NetID_true %>% filter(is_metabolite == "Yes"))
    table(all_ID_filter_artifact_NetID_false$is_metabolite)
    
    table(all_ID_filter$category_manual)
  }
  
  ## Artifact connections
  {
    artifact_connections = relation_list2 %>%
      filter(node1 %in% all_ID_filter$ID , node2 %in% all_ID_filter$ID) %>%
      filter(category != "biotransform") %>%
      filter(ILP_result != 0)
    table(artifact_connections$category)
  }
  
  ## unconnected peaks
  {
    ID_in_CPLEX = all_ID_filter %>% filter(ILP_result!=0) %>% pull(ID)
    all_unconnected = all_ID_filter %>%
      filter(!ID %in% ID_in_CPLEX) %>%
      distinct(ID, .keep_all = T)
  }
  
  ## all potential correct
  {
    all_labeled_peaks_filter_correct_potential = all_filter %>%
      filter(formula.x == `Ground truth`)
    
    all_labeled_peaks_filter_correct = all_labeled_peaks_filter_correct_potential %>%
      filter(ILP_result != 0)
    
    all_labeled_peaks_filter_wrong = all_ID_filter %>%
      filter(!ID %in% all_labeled_peaks_filter_correct$ID)
  }
  
    
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


# Evaluate all labeled peaks
{
  all_labeled_peaks = merge_result %>%
    filter(feature...11 != "Background") %>%
    arrange(-ILP_result) 
  
  all_labeled_peaks_distinct = all_labeled_peaks %>%
    distinct(ID, .keep_all=T)
  
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
# }
  
for(i in 1:length(NetID_files)){
  results = sapply(evaluate_summary[[i]], nrow)
  print(results)
}






selected_file = 1
test1 = evaluate_summary[[selected_file]]$Correct_potential %>%
  filter(!ID %in% evaluate_summary[[selected_file]]$Correct$ID) %>% 
  dplyr::select(ID) %>%
  inner_join(evaluate_summary[[selected_file]]$Total_potential_3e5)
result1 = evaluate_summary[[1]]$Total_potential_3e5 %>%
  dplyr::select(ID, ILP_result, formula.x)

result3 = evaluate_summary[[3]]$Total_potential_3e5 %>%
  dplyr::select(ID, ILP_result, formula.x)

result3_not1 = evaluate_summary[[3]]$Total_potential_3e5 %>%
  anti_join(result1)  %>%
  filter(formula.x == `Ground truth`, ILP_result ==1)

result1_not3 = evaluate_summary[[1]]$Total_potential_3e5 %>%
  anti_join(result3) %>%
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
  all_PAVE_unknowns_raw = merge_result %>%
    filter(feature...11 == "[]") %>%
    mutate(proposed_C = sapply(formula.x, elem_num_query, "C"),
           proposed_N = sapply(formula.x, elem_num_query, "N")) %>%
    arrange(-ILP_result) %>%
    distinct(ID, .keep_all =T)
  
  all_PAVE_unknowns_raw_unconnected = all_PAVE_unknowns_raw %>%
    filter(is.na(ILP_result) | ILP_result ==0)
  
  all_PAVE_unknowns = all_PAVE_unknowns_raw %>%
    filter(ILP_result !=0)
  
  
  
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




