
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
{
  data_select_ls = lapply(raw_ls, filter_data, medMz = 0, formula = "C7H15N1O1")
  # lapply(data_select_ls, View)
  lapply(data_select_ls, row.names)
  fig_ls = lapply(data_select_ls, plot_library_bar)
  pdf("C7H15N1O1.pdf",onefile = TRUE, w = 20, h = 10)
  print(fig_ls)
  dev.off()
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



