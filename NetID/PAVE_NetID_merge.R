library(janitor)
# install.packages("readxl")
library(readxl)
library(readr)
library(dplyr)
library(enviPat)
library(fitdistrplus)


# for(work_dir in c("Lin_Yeast_Pos", "Wenyun_Yeast_neg", "Wenyun_yeast_pos")){
  

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# work_dir = "Lin_Yeast_Pos"
print(work_dir)
setwd(work_dir)

# Merge PAVE, NetID and previous assignment results ####
{
  # Read NetID results
  NetID_files = list.files(pattern = "[0-9]{14}.csv")
  NetID_edge_files = list.files(pattern = "edge.csv")
  formula_list_ls = edge_info_sum_ls = list()
  for(i in 1:length(NetID_files)){
    formula_list_ls[[length(formula_list_ls)+1]] = read.csv(NetID_files[i], stringsAsFactors = F) %>%
      mutate(pred_C = sapply(formula, elem_num_query, "C"),
             pred_N = sapply(formula, elem_num_query, "N"))
    edge_info_sum_ls[[length(edge_info_sum_ls)+1]] = read.csv(NetID_edge_files[i], stringsAsFactors = F)
  }
}
{
  formula_list2 = formula_list_ls[[1]]
  relation_list2 = edge_info_sum_ls[[1]]

  # Read PAVE results
  PAVE_file = list.files(pattern = "^pks.*.xlsx")
  data(isotopes)
  PAVE_result = read_excel(PAVE_file) %>%
    mutate(Formula = replace_na(Formula, "")) %>%
    mutate(Formula = check_chemform(isotopes, Formula)$new_formula)

  # Merge results 
  merge_result = formula_list2 %>%
    dplyr::select(ID,formula,is_metabolite, Artifact_assignment, Biotransform, ILP_result, pred_C, pred_N, Input_id) %>%
    merge(PAVE_result, by.x="Input_id", by.y = "id", all.y = T) 
}


# Background ####
{
  NetID_background = merge_result %>%
    filter(is.na(ID))
  PAVE_background = merge_result %>%
    filter(!is.na(Background))
}


## Yeast experiment ####
{
  HMDB_library = read.csv("../dependent/HMDB_CHNOPS_clean.csv")
  # Get all biotransformed metabolites and match CN number
  biotransform_formulas = merge_result %>%
    filter(is_metabolite %in% c("Yes", "Maybe")) %>%
    filter(pred_C == C, pred_N == N) %>%
    mutate(In_HMDB = formula %in% HMDB_library$MF)
  
  table(biotransform_formulas$In_HMDB)
  
  # excluding HMDB and metlin
  biotransform_formulas_new = biotransform_formulas %>%
    filter(!In_HMDB) 
  
  biotransform_formulas_new_inten = biotransform_formulas_new %>%
    filter(sig > log10(5e4))
  
  output_new_formulas = biotransform_formulas_new_inten %>%
    dplyr::select(-pred_C, - pred_N, -Biotransform)
  
  tabyl(output_new_formulas, is_metabolite, Feature)
  
  write.csv(biotransform_formulas, "output_new_formulas.csv", row.names = F)
}



## NetID assignment for groundtruth ####
{
  Groundtruth = read_excel("Evaluate_set_unique.xlsx")
  
  merge_result_groundtruth = merge_result %>%
    full_join(Groundtruth %>% dplyr::select(Index,`Ground truth`, category_manual, NOTE)) 
  
  merge_result_groundtruth_unique = merge_result_groundtruth %>%
    filter(!is.na(category_manual)) %>%
    arrange(-ILP_result) %>%
    distinct(Input_id, .keep_all = T) 
  
  metabolite_assignment = merge_result_groundtruth_unique %>%
    filter(grepl("Metabolite", category_manual))
  
  metabolite_assignment_correct = metabolite_assignment %>%
    filter(formula == `Ground truth`) 
  metabolite_assignment_wrong = metabolite_assignment %>%
    filter(formula != `Ground truth`) 
  
  artifact_assignment = merge_result_groundtruth_unique %>%
    filter(grepl("Artifact|Fragment", category_manual))
  
  artifact_assignment_correct = artifact_assignment %>%
    filter(formula == `Ground truth`) 
  artifact_assignment_wrong = artifact_assignment %>%
    filter(formula != `Ground truth` | is.na(formula)) 
}

## Measurement error and distribution 
{
  Groundtruth_assigned = Groundtruth %>%
    filter(!category_manual %in% c("No_label", "Wrong Mass", "Unidentified")) %>%
    filter(!grepl("artifact", `Ground truth`)) %>%
    mutate(cal.mz = formula_mz(`Ground truth`)) %>%
    mutate(mz_dif = mz - cal.mz + formula_mz("H1", 1)) %>%
    mutate(mz_dif_ppm = mz_dif/cal.mz * 1e6)
  
  adjust_mass = Groundtruth_assigned %>%
    inner_join(formula_list2 %>% dplyr::select(Input_id, msr_mass )) 
  
  fit_normal_data = Groundtruth_assigned$mz_dif_ppm[abs(Groundtruth_assigned$mz_dif_ppm + 0.5)<2]+0.5
  
  fitdistData = fitdistrplus::fitdist(fit_normal_data, "norm")
  summary(fitdistData)
  plot(fitdistData)
  shapiro.test(fit_normal_data)
  
}



# ## All sig > 3e5, has C/N label ####
# {
#   Evaluate_set = merge_result %>%
#     # full_join(Groundtruth %>% dplyr::select(Index, mz,rt, `Ground truth`, category_manual, NOTE)) %>%
#     mutate(Feature = replace_na(Feature, "Unknown")) %>%
#     filter(sig > log10(3e5)) %>%
#     filter(is.na(Background), is.na(NonBio))
#   
#   # number of peaks
#   Evaluate_set_unique = Evaluate_set %>%
#     arrange(-ILP_result) %>%
#     distinct(Input_id, .keep_all = T)
#   
#   table(Evaluate_set_unique$category_manual) %>% sum()
# 
# }


## Rest #####
{
# Evaluate Xi's data from annotation ####

evaluate_summary = list()
for(i in 1:length(NetID_files)){
 
  
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
    
    
    
    evaluate_summary[[i]] = list(Total_potential_3e5 = all_ID_filter,
                                 Correct_potential = all_labeled_peaks_filter_correct_potential,
                                 Correct = all_labeled_peaks_filter_correct,
                                 Wrong = all_labeled_peaks_filter_wrong,
                                 # Unconnected = all_unconnected_non_backgrounds,
                                 Unconnected_3e5 = all_unconnected
    )
  }
  
  
  ## full dataset non_background evaluation
  {
    signal_cutoff = 1e5
    
    all_ID_CPLEX = merge_result %>%
      arrange(ID, -ILP_result) %>%
      filter(ILP_result != 0) %>%
      filter(log10_inten > log10(signal_cutoff))
    table(all_ID_CPLEX$is_metabolite)
    
    all_ID_CPLEX_nonbg = all_ID_CPLEX %>%
      filter(feature...11 != "Background") %>%
      mutate(CN_match = pred_C == C & pred_N == N) %>%
      mutate(In_HMDB = formula.x %in% HMDB_library$MF)
    
    tabyl(all_ID_CPLEX_nonbg, CN_match,is_metabolite, In_HMDB)
    
    all_CNmatch_notHMDB = all_ID_CPLEX_nonbg %>%
      filter(CN_match, !In_HMDB) %>%
      filter(is_metabolite %in% c("Yes")) # , "Maybe"
    
    write.csv(all_CNmatch_notHMDB, "New_metabolite.csv")
    
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
      filter(log10_inten>log10(1e5))
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
  
  
  
}

for(i in 1:length(NetID_files)){
  results = sapply(evaluate_summary[[i]], nrow)
  print(results)
}


selected_file = 3
test1 = evaluate_summary[[selected_file]]$Correct_potential %>%
  filter(!ID %in% evaluate_summary[[selected_file]]$Correct$ID) %>% 
  dplyr::select(ID) %>%
  inner_join(evaluate_summary[[selected_file]]$Total_potential_3e5)
result3 = evaluate_summary[[3]]$Total_potential_3e5 %>%
  dplyr::select(ID, ILP_result, formula.x)

result2 = evaluate_summary[[2]]$Total_potential_3e5 %>%
  dplyr::select(ID, ILP_result, formula.x)

result2_not3 = evaluate_summary[[2]]$Total_potential_3e5 %>%
  anti_join(result3)  %>%
  filter(formula.x == `Ground truth`, ILP_result ==1)

result3_not2 = evaluate_summary[[3]]$Total_potential_3e5 %>%
  anti_join(result2) %>%
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





}
