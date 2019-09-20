


# Main ####


# Related functions
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('~/myrepo/Meta_mine/metabo_library_functions.R')


# Read files

{
  setwd("~/myrepo/Meta_mine/library")
  filenames = list.files(list.dirs(recursive = T), pattern = "mdata.csv", full.names = T)
  filenames = filenames[!grepl("merge_mdata.csv", filenames)]
  filenames = filenames[grepl("pos", filenames)]
  # filenames = filenames[!grepl("yeast", filenames)]
  filenames = filenames[grepl("yeast", filenames)]
  filenames = filenames[5]
  
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

## Cluster dataset to see if LC condition are similar ####
{
  dataset_dist_mtrx = matrix(0,
                             nrow = length(filenames),
                             ncol = length(filenames)
  )
  rownames(dataset_dist_mtrx) = colnames(dataset_dist_mtrx) = dirname(filenames)

  for(i in 1:length(filenames)){
    print(Sys.time())
    for(j in 1:length(filenames)){
      dataset_dist_mtrx[i,j] = DatasetDist(raw_ls[[i]], raw_ls[[j]],
                                           log10_inten_cutoff = 5,  # only compare peaks with intensity above the threshold
                                           merge_group_ppm_tol = 3 # mz within ppm_tol will merge into one mz group
      )
    }
  }
  dataset_dist_mtrx_sym = 0.5 *(dataset_dist_mtrx + t(dataset_dist_mtrx))
  pdf("median_RT_shift_between_dataset.pdf")
  pheatmap::pheatmap(dataset_dist_mtrx_sym)
  dev.off()
}



# Search peak of interest based on medMz, medRt or formula
{
  data_select_ls = lapply(raw_ls, filter_data, medMz = 0, formula = "C3H9N1O1")
  # lapply(data_select_ls, View)
  lapply(data_select_ls, row.names)
  fig_ls = lapply(data_select_ls, plot_library_bar)
  pdf(paste("C3H7N1O1","_",timestamp(), ".pdf",sep=""),onefile = TRUE, w = 20, h = 10)
  print(fig_ls)
  dev.off()
}



# Filter peaks based on existing mz RT intensity list
{
  setwd("../project/190828 Yeast comparison")
  ion_mode = 1
  inten_cutoff = 5e4
  # raw_peak_list = read.csv("yeast_Pave_neg.csv", stringsAsFactors = F)
  raw_peak_list = read.csv("yeast_Pave_pos_unknown.csv", stringsAsFactors = F)
  table(raw_peak_list$feature)
  
  peak_list = raw_peak_list %>%
    filter(feature =="[]") %>%
    filter(sig>log10(inten_cutoff)) %>%
    rename(medMz = mz, medRt = RT)
  
  target_mdata = raw_ls[[1]] %>%
    mutate(mean_inten2 = rowMeans(.[,grepl("12C14N.0ev",colnames(.))]))
  
  data_select_ls = list()
  find_i = c()
  i=2
  for(i in 1:nrow(peak_list)){
    medMz = peak_list$medMz[i] - 1.007276 * ion_mode
    medRt = peak_list$medRt[i]
    data_select = filter_data(target_mdata, 
                                 medMz = medMz, medRt = medRt, delta_rt = 1,
                                 formula = "")
    if(nrow(data_select)!=0){find_i = c(find_i, i)}
    data_select_ls[[i]] = data_select
  }
  mdata_filter = bind_rows(data_select_ls)
  
  mdata_filter2 = mdata_filter %>%
    filter(mean_inten2>inten_cutoff | mean_inten>inten_cutoff) 
  
  mdata_filter2unique = mdata_filter2 %>%
    distinct(ID)
  
  
  find_i2 = find_i[sapply(data_select_ls[find_i], function(x){max(x$mean_inten, x$mean_inten2) >inten_cutoff})]
  PAVE_filter_pos = peak_list[find_i2,]
  
  write.csv(mdata_filter2, "Yeast_WL_PAVEfiltered_pos_mdata_5e4.csv", row.names = F)
  write.csv(PAVE_filter_pos, "Yeast_PAVEfiltered_pos_mdata_5e4.csv", row.names = F)
}




# Search peak of interest based on peak list
{
  setwd("../project/190828 Yeast comparison")
  peak_list = read.csv("yeast_Pave_neg_unknown.csv", stringsAsFactors = F)
  # dir.create("library_plots", showWarnings=F)
  # setwd("./library_plots")
  ion_mode = -1

  for(i in 1: nrow(peak_list)){
    medMz = peak_list$medMz[i] - 1.007276 * ion_mode
    potential_formula = peak_list$formula[i]
    data_select_ls = lapply(raw_ls, filter_data, medMz = medMz, formula = "")
    if(max(sapply(data_select_ls,nrow))==0){next}
    test = lapply(data_select_ls, colnames)
    fig_ls = lapply(data_select_ls, plot_library_bar)
    
    pdf(paste(potential_formula,"_",timestamp(), ".pdf",sep=""),onefile = TRUE, w = 21, h = 14)
    print(fig_ls)
    dev.off()
  }
  setwd("..")
}

# find correlation peaks
{
  tissues = correlated_peaks(mdata = raw_ls[[1]],
                             target_id = 490,
                             top_n = 10,
                             impute_options = c("threshold", "percentile"),
                             impute_options2 = c(T,F), # fill with random(T) or fixed(F)
                             normalize_options = c("row_mean"),
                             transform_options = c("log10",""),
                             scale_options = c("mean_center"),
                             limit_to_library = F,
                             print_pdf = F
  )
  pdf("C3H9N1O1_correlated2.pdf",onefile = TRUE, w = 20, h = 10)
  print(tissues$top$figure)
  dev.off()
}

# Plot heatmap
{
  mdata = raw_ls[[3]]
  # mdata = mdata[!is.na(mdata$library_match_name),]
  mdata = mdata[!duplicated(mdata$ID),]
  # mdata = mdata[mdata$log10_inten>4.5,]
  # mdata = mdata[mdata$`_log10_FDR`>200,]
  # mdata_row_name = paste(mdata$library_match_formula, mdata$medMz, mdata$medRt, sep="_")
  # mdata_row_name = paste(mdata$library_match_formula, mdata$library_match_name, mdata$ID, sep="_")
  mdata_row_name = 1:nrow(mdata)
  row.names(mdata) = mdata_row_name
  
  mdata_pre_clean = mdata[,14:(ncol(mdata)-2)]
  # mdata_pre_clean = mdata[,grepl("Kidney", colnames(mdata))]
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
    data_transform(transform_method = "log2") %>%
    data_scale(scale_method = "")
  
  
  # method "ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"
  # distance "correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  # correlation is default as pearson correlation
  row_cluster = cluster_mat(mdata_clean, method = "complete", distance = "correlation")
  col_cluster = cluster_mat(mdata_clean, method = "complete", distance = "correlation")
  my_plot_heatmap(mdata_clean = mdata_clean,
                  cohort = cohort,
                  imgName = "Plot", 
                  format = "pdf", # pdf
                  dpi = 72,
                  width = NA, # define output graph width
                  palette = "RdBu",  # RdBu, gbr, heat, topo
                  viewOpt = "overview", # Detail
                  rowV = row_cluster, # cluster by row, "F" if not clutser
                  colV = col_cluster, # cluster by column, "F" if not clutser
                  border = T, # border for each pixel
                  grp.ave = F, # group average
                  scale_ub = 3, scale_lb = 3 # heatmap scale, auto if ub = lb
  )
}



# Sample correlation matrix plot
{
  cormat <- cor(mdata_clean)
  # cormat <- as.matrix(dist(t(mdata_clean), method = "minkowski", p = .5))
  # cormat <- as.matrix(dist(t(mdata_clean), method = "euclidean"))
  # cormat <- as.matrix(dist(t(mdata_clean), method = "manhattan"))
  rownames(cormat) <- colnames(cormat) <- colnames(mdata_clean)
  
  annotation <- data.frame(class = cohort)
  rownames(annotation) <- rownames(cormat)
  
  if(length(unique(cohort)) < 9){
    pal9 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
             "#FFFF33", "#A65628", "#F781BF", "#999999")
    dist.cols = pal9[1:length(unique(cohort))]
  } else {
    pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
              "#FFFF99", "#B15928")
    dist.cols <- colorRampPalette(pal12)(length(unique(cohort)))
  }
  # barplot(rep(1,length(dist.cols)), col=dist.cols)
  cols = dist.cols[as.numeric(cohort)]
  uniq.cols <- unique(cols)
  cls <- cohort
  # names(uniq.cols) <- unique(as.character((cls)))
  names(uniq.cols) <- unique(as.character(sort(cls)))
  ann_colors <- list(class = uniq.cols)
  
  ## Perform the cluster analysis
  pheatmap(cormat, annotation = annotation,
           annotation_color = ann_colors)
}


# tsne
{
  library(tsne)
  ecb = function(x,y){ plot(x,t='n'); text(x,labels=cohort, col=cols) }
  tsne_run = tsne(mdata_clean[1:150,], 
                  k=2,
                   epoch_callback = ecb,
                   perplexity=50) 
}
