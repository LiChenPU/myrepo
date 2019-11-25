library(janitor)
# install.packages("readxl")
library(readxl)
library(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(enviPat)
library(fitdistrplus)
library(lc8)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(scales)

my_palette = c(brewer.pal(4, "Set3"), rep("#666666", 50))
# Read files ####
{
  Pave_merge_ls = list()
  NetID_merge_ls = list()
  NetID_edge_merge_ls = list()
  merge_ls = list()
  input_files = c("Lin_Yeast_Neg"
                  # "Lin_Yeast_Pos"
                  # "Wenyun_Yeast_neg",
                  # "Wenyun_yeast_pos"
  )
  for(work_dir in input_files){
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
    # work_dir = "Lin_Yeast_Pos"
    print(work_dir)
    setwd(work_dir)
    
    PAVE_file = list.files(pattern = "^pks.*.xlsx")

    if(grepl("neg", work_dir, ignore.case = T)){
      ionization = -1
    } else if(grepl("pos", work_dir, ignore.case = T)){
      ionization = 1
    }
    
    PAVE_data = read_xlsx(PAVE_file) %>%
      mutate(origin = work_dir) %>%
      mutate(ionization = ionization) %>%
      mutate(neutral_mz = mz - ionization * formula_mz("H1", 1)) %>%
      mutate(Feature = replace_na(Feature, "Unknown"))
      
    Pave_merge_ls[[length(Pave_merge_ls)+1]] = PAVE_data
    
    NetID_files = list.files(pattern = "[0-9]{14}.csv")
    NetID_data = read.csv(NetID_files, stringsAsFactors = F)
    NetID_merge_ls[[length(NetID_merge_ls)+1]] = NetID_data
    
    NetID_edge_files = list.files(pattern = "edge.csv")
    NetID_edge_data = read.csv(NetID_edge_files, stringsAsFactors = F)
    NetID_edge_merge_ls[[length(NetID_edge_merge_ls)+1]] = NetID_edge_data
    
    # 
    # 
    # merge_data = merge(NetID_data %>% filter(ILP_result != 0), 
    #                    PAVE_data, by.x="Input_id", by.y = "id", all.y = T)
    # merge_ls[[length(merge_ls)+1]] = merge_data
  }
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

{
  NetID_data = NetID_merge_ls[[1]] %>%
    filter(ILP_result != 0)
  Groundtruth = read_xlsx("./Lin_Yeast_Neg/4-4-Table S10-Annotation of all peaks detected in S. cerevisiae and E. coli-xi_LC.xlsx",
                       sheet = "Yeast-neg-truth") %>%
    arrange(ID...1) %>%
    dplyr::select(-c(36:67))
    
    # rename_all(funs(paste0("old_",.))) %>%
    # mutate(origin = "Old_PAVE")
 
  
  all = merge(NetID_data, Groundtruth, by.x = "Input_id", by.y = "ID...1", all.y=T) %>%
    setNames(make.names(gsub("\\..*","", names(.)), unique = TRUE)) %>%
    mutate(ionization = ionization) %>%
    mutate(neutral_mz = mz.1 - ionization * formula_mz("H1", 1)) %>%
    # mutate(neutral_mz = mz) %>%
    mutate(Feature = replace_na(Feature, "Unknown"))
}

## Group peaks with simialr mz and rt ####
{
  # all = bind_rows(merge_ls) %>%
    # mutate(library_id = 1:nrow(.))
  
  ms_dif_ppm = 5e-6
  rt_dif_min = 2
  
  ##Group MS groups
  {
    
    all = all %>% 
      arrange(neutral_mz)
    mzs = all$neutral_mz
    
    count = 1
    MZ_group = rep(1,(length(mzs)))
    for(i in 2:length(mzs)){
      if(mzs[i]-mzs[i-1]>mzs[i-1]*ms_dif_ppm){
        count = count+1
      }
      MZ_group[i]=count
    }
    all = all %>%
      mutate(MZ_group = MZ_group)
  }
  
  ##Group RT similar groups based on MS groups
  {
    all = all %>%
      arrange(MZ_group, RT.1)
    rts = all$RT.1
    
    MZRT_group = rep(1,(length(rts)))
    
    count = 1
    for(i in 2:length(rts)){
      if(MZ_group[i]!=MZ_group[i-1] | rts[i]-rts[i-1]>rt_dif_min){
        count = count+1
      }
      MZRT_group[i]=count
    }
    all["MZRT_group"] = MZRT_group
  }
  
  all = all %>%
    arrange(ID)
  
  HMDB_library = read.csv("./dependent/HMDB_CHNOPS_clean.csv")
  all = all %>%
    mutate(pred_C = sapply(formula, elem_num_query, "C"),
           pred_N = sapply(formula, elem_num_query, "N")) %>%
    mutate(CN_match = (pred_C == C & pred_N == N)) %>%
    mutate(In_HMDB = formula %in% HMDB_library$MF) %>%
    mutate(Intensity_group = cut(sig, breaks = c(0, 5, Inf))) %>%
    mutate(transformation_category = case_when(
      is_metabolite == "Yes" ~ "Biotransformed only",
      is_metabolite == "No" ~ "Non-biotransformed only",
      is_metabolite == "Maybe" ~ "Likely bio or non-biotransformed",
      is.na(is_metabolite) ~ "Not assigned"
    ))
}

## Assign status based on NetID results ####
# {
#   write.csv(all, "all.csv", row.names = F)
# 
#   all_inten = all %>%
#     filter(sig > log10(1e5))
#   all_background_nonbio = all_inten %>%
#     # filter(is_na(ID) & Background %in% c("Background", "NonBio", "Artifact"), feature == "'Background'")
#     filter(feature %in% c("'Background'", "'Low_score'")) %>%
#     mutate(Status = feature)
# 
#   table(all_background_nonbio$Feature)
#   table(all_background_nonbio$feature)
# 
#   all_biological = all_inten %>%
#     filter(!Input_id %in% all_background_nonbio$Input_id)
#   table(all_biological$Feature)
#   table(all_biological$feature)
# 
#   {
#     all_adduct = all_biological %>%
#       filter(feature == "'Adduct'") %>%
#       mutate(Status = case_when(
#         grepl("Adduct", Artifact_assignment)~"Adduct",
#         Artifact_assignment!="" ~ "Other artifacts",
#         # Background == "Background" ~ "High background",
#         is.na(formula) ~ "No formula assigned",
#         Artifact_assignment=="" ~ "Biotransformed"
#       ))
# 
#     all_adduct_np = all_biological %>%
#       filter(feature == "'Adduct_np'") %>%
#       mutate(Status = case_when(
#         grepl("Adduct", Artifact_assignment)~"Adduct",
#         Artifact_assignment!="" ~ "Other artifacts",
#         # Background == "Background" ~ "High background",
#         is.na(formula) ~ "No formula assigned",
#         Artifact_assignment=="" ~ "Biotransformed"
#       ))
# 
#     all_dimer = all_biological %>%
#       filter(feature == "'Dimer'") %>%
#       mutate(Status = case_when(
#         grepl("Oligomer", Artifact_assignment)~"Oligomer",
#         Artifact_assignment!="" ~ "Other artifacts",
#         # Background == "Background" ~ "High background",
#         is.na(formula) ~ "No formula assigned",
#         Artifact_assignment=="" ~ "Biotransformed"
#       ))
# 
#     all_fragment = all_biological %>%
#       filter(feature == "'Fragment'") %>%
#       mutate(Status = case_when(
#         grepl("Fragment", Artifact_assignment)~"Fragment",
#         Artifact_assignment!="" ~ "Other artifacts",
#         # Background == "Background" ~ "High background",
#         is.na(formula) ~ "No formula assigned",
#         Artifact_assignment=="" ~ "Biotransformed"
#       ))
# 
# 
#     all_isotope = all_biological %>%
#       filter(feature == "'Isotope'") %>%
#       mutate(Status = case_when(
#         grepl("Isotope", Artifact_assignment)~"Isotope",
#         Artifact_assignment!="" ~ "Other artifacts",
#         # Background == "Background" ~ "High background",
#         is.na(formula) ~ "No formula assigned",
#         Artifact_assignment=="" ~ "Biotransformed"
#       ))
# 
#     all_multicharge = all_biological %>%
#       filter(feature == "'Multicharge'") %>%
#       mutate(Status = case_when(
#         grepl("Double_charge", Artifact_assignment)~"Double_charge",
#         Artifact_assignment!="" ~ "Other artifacts",
#         # Background == "Background" ~ "High background",
#         is.na(formula) ~ "No formula assigned",
#         Artifact_assignment=="" ~ "Biotransformed"
#       ))
# 
# 
#     all_metabolite = all_biological %>%
#       filter(feature == "'Metabolite'") %>%
#       mutate(Status = ifelse(Biotransform, "Biotransformed derived", "Biotransform_Not_matched"))
#   }
# 
# 
# 
#   all_unknown = all_biological %>%
#     filter(feature=="'Low_C'" | feature == "[]") %>%
#     mutate(Status = transformation_category) %>%
#     # mutate(Status = ifelse(Background == "Background", "High background", Status))
#     filter(T)
# 
#   all_bind = bind_rows(all_adduct,
#                        all_adduct_np,
#                        all_dimer,
#                        all_fragment,
#                        all_isotope,
#                        all_multicharge,
#                        all_metabolite,
#                        all_unknown) %>%
#     arrange(Input_id)
# 
#   write.csv(all_bind %>%
#               bind_rows(all) %>%
#               distinct(Input_id, .keep_all = T) %>%
#               arrange(Input_id)
#             , "all_bind.csv", row.names = F)
# 
#   table(all_bind$Status)
# 
# }



## Evaluate manual assignment results ####
{
  all_bio = all %>%
    filter(sig > log10(1e5)) %>%
    filter(!feature %in% c("'Background'", "'Low_score'")) %>%
    # filter(Formula_validated == "Y") %>%
    mutate(Formula_match = case_when(
      Formula_validated == "?" ~ "Ground truh unassigned",
      is.na(formula) ~ "No assignment",
      formula == Formula.1 ~ "Correct",
      formula != Formula.1 ~ "Incorrect"
    )) %>%
    mutate(Category3 = case_when(
      feature == "'Metabolite'" ~ "Metabolite",
      feature == "[]" ~ "Unknown",
      feature != "" ~ "Artifact"
    ))
  
  table(all_bio$Category3)
  tabyl(all_bio, Category3, Formula_match)
  
  all_bio_filter = all_bio %>%
    filter(Category3 == "Unknown") %>%
    filter(!is.na(Status_validated)) %>%
    filter(T)
  
  tabyl(all_bio_filter, Status_validated)
  
  all_bio_filter = all_bio %>%
    # filter(CN_match) %>%
    filter(feature %in% c("'Adduct'", "'Isotope'", "'Adduct_np'", "'Dimer'", "'Fragment'", "'Metabolite'")) %>%
    # filter(feature %in% c("'Dimer'", "'Fragment'", "'Metabolite'")) %>%
    # filter(is.na(CN_match) | !CN_match) %>%
    # filter(Formula_match %in% c("Incorrect","No assignment")) %>%
    # filter(Category3 == "Unknown") %>%
    # filter(Unknown_assignment == "Non-metabolite") %>%
    # filter(Status_validated == "Adduct") %>%
    # filter(Note != "NetID_mistake") %>%
    # filter(!is.na(Note)) %>%
    filter(T)
    
  tabyl(all_bio_filter, CN_match, feature)
  
  {
    NetID_lin_neg = NetID_merge_ls[[1]]
    NetID_edge_lin_neg = NetID_edge_merge_ls[[1]]
    
    match_parent_CN = function(selected_node){
      selected_formula = NetID_lin_neg %>%
        filter(ID == selected_node)
      selected_edge = NetID_edge_lin_neg %>%
        filter(node2 == selected_node) %>%
        filter(category != "biotransform", !grepl("\\[", category)) %>%
        filter(!category %in% c("Oligomer", "Heterodimer")) %>%
        filter(ILP_result != 0) %>%
        filter(T)
      
      selected_CN = all %>%
        filter(ID == selected_node) %>%
        dplyr::select(C,N)
      parent_CN = all %>%
        filter(ID %in% selected_edge$node1) %>%
        filter(CN_match) %>%
        # filter(ID == -1) %>%
        dplyr::select(C,N)
      
      matched_rows = match_df(parent_CN, selected_CN)
      return(nrow(matched_rows) != 0)
      
    }
    
    all_bio_adduct = all_bio %>%
      filter(feature %in% c("'Adduct'")) %>%
      mutate(CN_match_parent = sapply(ID, match_parent_CN))
    
    all_bio_adduct_filter = all_bio_adduct %>%
      filter(!CN_match_parent)
    tabyl(all_bio_adduct, CN_match_parent)
  }
  
  
  tabyl(all_bio_filter, CN_match, feature)
  
  tabyl(all_bio_filter$Status_validated)
  
  colnames(all_bio)
  tabyl(all_bio$Status)
  tabyl(all_bio$Formula_validated)
  tabyl(all_bio$Status_validated)
  
}

## Query NetID ####
{
  NetID_lin_neg = NetID_merge_ls[[1]]
  NetID_edge_lin_neg = NetID_edge_merge_ls[[1]]
  
  query_Input_id = 5225
  selected_formula = NetID_lin_neg %>%
    filter(Input_id == query_Input_id)
  selected_node = 1402
  selected_formula = NetID_lin_neg %>%
    filter(ID == selected_node)
  selected_edge = NetID_edge_lin_neg %>%
    filter(node1 == selected_node | node2 == selected_node)
  selected_ILP_id = 2960
  selected_ILP_edge = NetID_edge_lin_neg %>%
    filter(ILP_id1 == selected_ILP_id | ILP_id2 == selected_ILP_id)
}
