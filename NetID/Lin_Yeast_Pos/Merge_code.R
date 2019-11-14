



## Merge old PAVE ####
{
  new_PAVE = read.csv("raw_data.csv") %>%
    mutate(medRt = round(medRt, 2))
  old_PAVE = read_xlsx("../Lin_Yeast_Neg/4-4-Table S10-Annotation of all peaks detected in S. cerevisiae and E. coli-xi_LC.xlsx",
                       sheet = "yeast-pos") %>%
    rename_all(funs(paste0("old_",.))) %>%
    mutate(origin = "Old_PAVE") %>%
    mutate(old_RT = round(old_RT, 2))
  
  Merge_newold_PAVE = merge(new_PAVE, old_PAVE, by.x = c("medMz", "medRt"), by.y = c("old_mz", "old_RT"), all = T)%>%
    # distinct(mz.y, rt, .keep_all=T)
    filter(T)
  
  test = Merge_newold_PAVE %>%
    filter(is.na(metaGroupId) | is.na(old_feature))
  
  tabyl(Merge_newold_PAVE, old_feature, Feature)
  table(Merge_newold_PAVE$Feature)
  
  Merge_newold_PAVE_filter = Merge_newold_PAVE %>%
    filter(old_feature == "'Background'") %>%
    # filter(Feature == "Metabolite") %>%
    filter(Feature == "Unknown") %>%
    # filter(score.y > 0.5) %>%
    # filter(old_score > 0.5) %>%
    # filter(old_description == "'13C'") %>%
    filter(sig>log10(1e5)) %>%
    # mutate(score_dif = score.y - old_score) %>%
    arrange(-sig)%>%
    filter(T)
  
  hist(Merge_newold_PAVE_filter$score_dif)
  plot(Merge_newold_PAVE_filter$score.y, Merge_newold_PAVE_filter$old_score)
  
  temp = all_score %>%
    filter(formula == "C7H15N1O3")
  
}