library(janitor)
# install.packages("readxl")
library(readxl)
library(readr)
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
  input_files = c("Lin_Yeast_Neg", 
                  "Lin_Yeast_Pos",
                  "Wenyun_Yeast_neg",
                  "Wenyun_yeast_pos"
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
    
    
    
    merge_data = merge(NetID_data %>% filter(ILP_result != 0), 
                       PAVE_data, by.x="Input_id", by.y = "id", all.y = T)
    merge_ls[[length(merge_ls)+1]] = merge_data
  }
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}



## Group peaks with simialr mz and rt ####
{
  all = bind_rows(merge_ls) %>%
    mutate(library_id = 1:nrow(.))
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
      arrange(MZ_group, rt)
    rts = all$rt
    
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
    arrange(library_id)
  
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


## Fix score assignment ####
{
  score_cutoff = 0.75
  tabyl(all, Feature, origin)
  
  all_score = all %>%
    mutate(Feature = case_when(
      score.y < score_cutoff & Feature == "Unknown" ~ "NonBio",
      !is.na(Feature) ~ Feature
    ))
  tabyl(all, Feature, origin)
  tabyl(all_score, Feature, origin)
}

## Merge old PAVE ####
{
  new_PAVE = all_score %>%
    filter(origin == "Lin_Yeast_Neg")
  old_PAVE = read_xlsx("./Lin_Yeast_Neg/4-4-Table S10-Annotation of all peaks detected in S. cerevisiae and E. coli-xi_LC.xlsx",
                       sheet = "Yeast-neg") %>%
    rename_all(funs(paste0("old_",.))) %>%
    mutate(origin = "Old_PAVE")
  
  Merge_newold_PAVE = merge(new_PAVE, old_PAVE, by.x = c("mz.y", "rt"), by.y = c("old_mz", "old_RT"))%>%
    # distinct(mz.y, rt, .keep_all=T)
    filter(T)
  
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
## Ground truth for Lin data ####
{
  Merge_newold_PAVE %>% 
    arrange(old_ID) %>%
    write.csv("Merge_newold_PAVE.csv", row.names = F)
  
  Merge_newold_PAVE_filter = Merge_newold_PAVE %>%
    filter(!Feature %in% c("Background", "Nonbio") & !is.na(Feature)) %>%
    arrange(-sig)
  
  NetID_lin_neg = NetID_merge_ls[[1]]
  NetID_edge_lin_neg = NetID_edge_merge_ls[[1]]
  
  selected_node = 1157
  selected_formula = NetID_lin_neg %>%
    filter(ID == selected_node)
  selected_edge = NetID_edge_lin_neg %>%
    filter(node1 == selected_node | node2 == selected_node)
  selected_ILP_id = 1559
  selected_ILP_edge = NetID_edge_lin_neg %>%
    filter(ILP_id1 == selected_ILP_id | ILP_id2 == selected_ILP_id)
  
  
}

## Dataset Overlaps ####
{ 
  MZRT_group_filter = all_score %>%
    filter(Feature == "Metabolite") %>%
    # filter(! Feature %in% c("Artifact", "Background", "NonBio")) %>%
    # filter(Intensity_group == levels(Intensity_group)[2]) %>%
    pull(MZRT_group)
  
  all_filter = all_score %>%
    filter(Feature == "Metabolite") %>%
    filter(MZRT_group %in% MZRT_group_filter) 
    
  file_mzrtgroup = all_filter %>%
    mutate_if(is.character, as.factor) %>%
    group_by(origin, MZ_group, .drop = FALSE) %>%
    dplyr::select()
  
  coverage_matrix = matrix(nrow = nlevels(file_mzrtgroup[[1]]), 
                           ncol = nlevels(file_mzrtgroup[[1]]))
  
  
  temp_summary = t(table(file_mzrtgroup))
  colnames(coverage_matrix) = rownames(coverage_matrix) = colnames(temp_summary)
  
  for(i in 1:ncol(temp_summary)){
    temp = t(temp_summary[,i])
    temp[temp>1] = 1
    for(j in 1:ncol(temp_summary)){
      coverage_matrix[i,j] = temp %*% temp_summary[,j]
    }
  }
  print(coverage_matrix)
}

## PAVE summary ####
fig_ls = list()
for(origin_file in input_files){
{
  # origin_file = input_files[[1]]
  print(origin_file)
  final_step = all_score %>%
    filter(origin == origin_file) %>%
    filter(Feature %in% c("Unknown", "Metabolite"))
  
  level1 = final_step %>%
    mutate_if(is.character, as.factor) %>%
    group_by(Intensity_group, .drop = FALSE) %>%
    summarise(n=n())
  
  level2 = final_step %>%
    mutate_if(is.character, as.factor) %>%
    group_by(Intensity_group, Feature, .drop = FALSE) %>%
    summarise(n=n())
  
  level2_2 = final_step %>%
    mutate_if(is.character, as.factor) %>%
    group_by(Intensity_group, transformation_category, .drop = FALSE) %>%
    summarise(n=n())
  
  level3 = final_step %>%
    mutate_if(is.character, as.factor) %>%
    group_by(Intensity_group, Feature, transformation_category, .drop = FALSE) %>%
    summarise(n=n())
  # print(level3)
  # base::crossprod(table(test))
}

## ggplot2 
{
  fig_level2 = ggplot(level2, aes(x = Feature, y = n, fill = Intensity_group)) + 
    geom_bar(stat = "identity",
             position = position_dodge(),
             colour = "#333333"
    ) +
    labs(#titles = "Intensity group (log10)", 
         y = "Number") + 
    guides(fill = guide_legend(
      title = "Intensity group (log10)",
      reverse = F
    )) +
    scale_y_continuous(expand = c(0,0)) + 
    # facet_wrap(~cohorts) +
    theme_classic(base_size = 12 # edit font size for all non-data text
    ) +
    
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
          axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7),
          axis.ticks.x = element_blank()
    )
  print(fig_level2)
  
  
  fig_level3 = ggplot(level3, aes(x = Feature, y = n, fill = transformation_category)) + 
    geom_bar(stat = "identity",
             position = position_dodge(),
             colour = "#333333"
    ) +
    labs(title = "Intensity group (log10)", 
         y = "Number") + 
    guides(fill = guide_legend(
      title = "Formula assignment from NetID",
      reverse = F
    )) +
    scale_y_continuous(expand = c(0,0)) + 
    scale_fill_manual(values = my_palette[c(1:3,5)]) +
    facet_wrap(~Intensity_group, scales = "free_y") +
    theme_classic(base_size = 12 # edit font size for all non-data text
    ) +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
          axis.text.x = element_text(angle = 0, hjust = .5, vjust = .7),
          axis.ticks.x = element_blank()
    )
  print(fig_level3)
}
  
  
fig_ls[[length(fig_ls)+1]] = list(fig_level2, fig_level3)
}

names(fig_ls) = input_files

{
  ggarrange(
    ggarrange(fig_ls[[1]][[1]],fig_ls[[1]][[2]],ncol = 2, widths = c(1,2), labels = input_files[1]),
    ggarrange(fig_ls[[2]][[1]],fig_ls[[2]][[2]],ncol = 2, widths = c(1,2), labels = input_files[2]),
    ggarrange(fig_ls[[3]][[1]],fig_ls[[3]][[2]],ncol = 2, widths = c(1,2), labels = input_files[3]),
    ggarrange(fig_ls[[4]][[1]],fig_ls[[4]][[2]],ncol = 2, widths = c(1,2), labels = input_files[4]),
    nrow = 4) %>%
    ggexport(filename = "PAVE_NetID_merge.pdf", width = 15, height = 10)
}


