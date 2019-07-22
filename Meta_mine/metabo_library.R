
# Main ####

# Read files
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("./library")
  filenames = list.files(list.dirs(recursive = F), pattern = "mdata.csv", full.names = T)
  filenames = filenames[grepl("pos", filenames)]
  
  num_of_files = length(filenames)
  raw_ls = list()
  for(i in 1:num_of_files){
    filename=filenames[i]
    df_temp = read_csv(filename)
    df_temp = cbind(label = gsub("\\./*","",dirname(filename)),df_temp, stringsAsFactors = F)
    df_temp$medMz = round(df_temp$medMz, 4)
    df_temp$medRt = round(df_temp$medRt, 3)
    raw_ls[[i]]= df_temp
  }
  rm(df_temp)
}

# Search peak of interest based on medMz, medRt or formula
data_select_ls = lapply(raw_ls, filter_data, medMz = 0, formula = "C3H9N1O1")
# lapply(data_select_ls, View)
lapply(data_select_ls, row.names)
fig_ls = lapply(data_select_ls, plot_library_bar)
pdf("C3H9N1O1_short.pdf",onefile = TRUE, w = 20, h = 4)
print(fig_ls)
dev.off()




tissues = correlated_peaks(mdata = raw_ls[[1]],
                        target_id = 19977,
                        top_n = 3,
                        impute_options = c("threshold", "percentile"),
                        impute_options2 = c(T,F), # fill with random(T) or fixed(F)
                        normalize_options = c("row_mean"),
                        transform_options = c("log10",""),
                        scale_options = c("mean_center"),
                        print_pdf = F
                        )
pdf("C3H9N1O1_correlated.pdf",onefile = TRUE, w = 20, h = 10)
print(tissues$top$figure)
dev.off()


mdata = raw_ls[[1]]
# mdata = mdata[!is.na(mdata$library_match_name),]
mdata = mdata[!duplicated(mdata$ID),]
# mdata = mdata[mdata$log10_inten>4.5,]
# mdata = mdata[mdata$`_log10_FDR`>200,]
# mdata_row_name = paste(mdata$library_match_formula, mdata$medMz, mdata$medRt, sep="_")
mdata_row_name = paste(mdata$library_match_formula, mdata$library_match_name, mdata$ID, sep="_")
mdata_row_name = 1:nrow(mdata)
row.names(mdata) = mdata_row_name

mdata_pre_clean = mdata[,14:(ncol(mdata)-2)]
all_names = colnames(mdata_pre_clean)
if(length(grep("blank|blk", all_names, ignore.case = T))!=0){
  sample_names=all_names[-grep("blank|blk", all_names, ignore.case = T)]
} else {
  sample_names=all_names
}
cohort = stri_replace_last_regex(sample_names,'_\\d+|-\\d+|\\d+', '',
                                 stri_opts_regex(case_insensitive=T))
cohort = as.factor(cohort)
mdata_pre_clean = mdata_pre_clean[,sample_names]
mdata_clean = mdata_pre_clean %>% 
  data_impute(impute_method = "threshold", random = F) %>%
  data_normalize(nor_method = "") %>%
  data_transform(transform_method = "log10") %>%
  data_scale(scale_method = "mean_center")




data_topn_tissues = read.csv("../old_files/data_topn_tissues.csv", stringsAsFactors = F)
data_signif_counts = read.csv("../old_files/data_signif_counts.csv", stringsAsFactors = F)
data_topn_tissues = data_topn_tissues[!duplicated(data_topn_tissues),]


temp_ls = list()
for(i in 1:nrow(data_topn_tissues)){
  medMz = data_topn_tissues$medMz[i]
  medRt = data_topn_tissues$medRt[i]
  plot_data = filter_data(mdata, 
                          medMz = medMz,
                          medRt = medRt, 
                          formula = ""
                          )
  temp_ls[[length(temp_ls)+1]] = row.names(plot_data)[1]
}
topn_select = unlist(temp_ls)


my_plot_heatmap(mdata_clean = mdata_clean,
                           cohort = cohort,
                           imgName = "Test", 
                           format = "pdf", # pdf
                           dpi = 72,
                           width = NA, # define output graph width
                           palette = "RdBu",  # RdBu, gbr, heat, topo
                           viewOpt = "overview", # Detail
                           rowV = T, colV = T, # cluster by row/column
                           border = T, # border for each pixel
                           grp.ave = T, # group average
                           scale_ub = 4, scale_lb = -4 # heatmap scale
)

write.csv(mdata[topn_select,], "mdata_top.csv", row.names = F)
