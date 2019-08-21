
# Import library ####
{
  # install.packages(c("readr","tidyr","dplyr","stringi","tictoc","profvis","ggplot2",
  #                    "ggrepel","RColorBrewer", "pheatmap"))
  # devtools::install_github("LiChenPU/Formula_manipulation")
  library(readr)
  library(tidyr)
  library(dplyr)
  library(stringi)
  library(tictoc)
  library(lc8)
  library(profvis)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  
}

# Function ####
# Timestamp ####
timestamp = function(){
  printtime = Sys.time()
  timestamp <- paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
  return(timestamp)
}

# plot_library_bar ####
plot_library_bar = function(data_select){
  if(nrow(data_select) ==0){return(NULL)}
  # data_select = data_select_ls[[4]]
  
  delColNames = c("ID","ILP_result","is_metabolite","_log10_FDR","library_match_formula",
                  "library_match_name","high_blank","goodPeakCount","maxQuality","parent",
                  "mean_inten","log10_inten")
  
  # if any colnames not found, use dplyr::select(-one_of(delColNames))
  data_plot = data_select %>%
    dplyr::select(-delColNames) %>%
    gather(key = names, value = number, -formula, -label, -medMz, -medRt) %>%
    mutate(cohort = stri_replace_last_regex(names,'\\d+|_\\d+|-\\d+|-a|-b|-c|_mean', '',stri_opts_regex(case_insensitive=T))) %>%
    filter(complete.cases(number))
  
  
  unique_cohort = unique(data_plot$cohort)
  data_ls = list(length = length(unique_cohort))
  select_col = c("label","formula","medMz","medRt","cohort","TIC","sd")
  for(i in 1:length(unique_cohort)){
    temp = data_plot[data_plot$cohort==unique_cohort[i],]
    temp_spread = spread(temp, key = names, value = number)
    if(ncol(temp_spread) == 6){
      temp_spread["TIC"] = temp_spread[,6]
      temp_spread["sd"] = 0
    } else{
      temp_spread["TIC"] = apply(temp_spread[,6:(ncol(temp_spread))], 1, mean)
      temp_spread["sd"] = apply(temp_spread[,6:(ncol(temp_spread)-1)], 1, sd)
    }
    
    data_ls[[i]] = temp_spread[, select_col]
  }
  
  data_bind = bind_rows(data_ls)
  data_bind$formula = factor(data_bind$formula, levels = unique(data_select$formula))
  data_bind = data_bind[order(data_bind$formula),]
  data_bind = unite(data_bind, "title", c(formula, medMz, medRt))
  # data_bind = data_bind[with(data_bind, order(-TIC)),]
  data_bind$title = factor(data_bind$title, levels = unique(data_bind$title))
  
  if(length(unique(unique_cohort)) < 9){
    pal9 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
             "#FFFF33", "#A65628", "#F781BF", "#999999")
    dist.cols = pal9[1:length(unique(unique_cohort))]
  } else {
    pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
              "#FFFF99", "#B15928")
    dist.cols <- colorRampPalette(pal12)(length(unique(unique_cohort)))
  }
  
  
  figure <- ggplot(data_bind, aes(x = cohort, y = TIC, fill = cohort))+
    geom_bar(position=position_dodge(),stat='identity', color = "black") +
    geom_errorbar(aes(ymin=TIC-sd, ymax=TIC+sd),width = 0.5, position=position_dodge(.9)) +
    ggtitle(data_select$label[1])+
    theme(legend.title=element_blank()) +
    facet_wrap(~title, scales = "free_y")+
    # scale_y_continuous(limits = c(-0.2,1.2),breaks = c(0,0.5,1)) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    # scale_fill_brewer(palette ="Set3", guide = guide_legend(reverse = TRUE))
    scale_fill_manual(values = dist.cols)
  
  return(figure)
}


# filter_data ####
filter_data = function(data,  medMz = 0,
                       medRt = 0,
                       delta_mz = 0.001,
                       delta_rt = 0.1,
                       formula = "")
{
  data_select = data
  if(medMz != 0){
    data_select = data_select[abs(data_select$medMz - medMz) < delta_mz,]
  }
  if(medRt != 0){
    data_select = data_select[abs(data_select$medRt - medRt) < delta_rt,]
  }
  if(formula != ""){
    data_select = data_select[data_select$formula == formula & !is.na(data_select$formula),]
  }
  return (data_select)
}

# data_impute - Substitute low signal and NA into small number ####
data_impute = function(pre_impute, 
                       impute_method = "threshold", # "min", "percentile"
                       random = F, # replacement is a fixed numberr or a random between 0 to threshold
                       inten_thrshld = 2500, 
                       data_percentile = 0.01  # Only used if method is "percentile
) 
{
  if(impute_method == "min"){
    inten_thrshld = min(pre_impute[pre_impute > 0], na.rm = T)
  } else if(impute_method == "percentile"){
    inten_thrshld = quantile(pre_impute[pre_impute > 0], data_percentile, na.rm = T)  
  } else if(impute_method == "threshold"){
    inten_thrshld = inten_thrshld
  }  
  subst_pstn = pre_impute<=inten_thrshld | is.na(pre_impute)
  post_impute = pre_impute
  if(random){
    post_impute[subst_pstn] = runif(sum(subst_pstn),min = 1, max = inten_thrshld)
  } else{
    post_impute[subst_pstn] = inten_thrshld
  }
  return(post_impute)
}


# Normalize ####
data_normalize = function(pre_norm, 
                          nor_method = "row_mean", 
                          sample_name = "", 
                          cohort_name = "")
{
  if(nor_method == "row_median"){
    row_median = apply(pre_norm, 1, median, na.rm=T)
    norm_factor = row_median
  } else if(nor_method == "row_mean"){
    row_mean = apply(pre_norm, 1, mean, na.rm=T)
    norm_factor = row_mean
  } else if(nor_method == "sample"){
    sample_pos = which(sample_name == colnames(pre_norm))
    norm_factor = pre_norm[,sample_pos]
  } else if(nor_method == "cohort"){
    col_name = colnames(pre_norm)
    cohort = stri_replace_last_regex(col_name,'_\\d+|-\\d+', '',stri_opts_regex(case_insensitive=T))
    sample_pos = which(cohort == cohort_name)
    norm_factor = rowMeans(pre_norm[,sample_pos])
  } else if(nor_method == "col_median"){
    col_median = apply(pre_norm, 2, median, na.rm=T)
    col_median = col_median/mean(col_median)
    norm_factor = col_median
    post_norm = sweep(pre_norm, 2, norm_factor, "/")
    return(post_norm)
  } else {
    print("Warnings: no normalize")
    return(pre_norm)
  }
  
  post_norm = sweep(pre_norm, 1, norm_factor, "/")
  return(post_norm)
}

# data_transform ####
data_transform = function(pre_trsf, transform_method = "log10")
{
  if(transform_method == "log10"){
    post_trsf = log10(pre_trsf)
  } else if(transform_method == "log2"){
    post_trsf = log2(pre_trsf)
  } else {
    print("Warnings: no transfrom")
    return(pre_trsf)
  }
  return(post_trsf)
}


# Scaling ####
data_scale = function(pre_scale, scale_method = "mean_center"){
  row_mean = rowMeans(pre_scale)
  
  if(scale_method == "mean_center"){
    post_scale = sweep(pre_scale, 1, row_mean, "-")
  } else if(scale_method == "auto"){
    post_scale = sweep(pre_scale, 1, row_mean, "-")
    row_sd = apply(pre_scale, 1, sd, na.rm = T)
    post_scale = sweep(pre_scale, 1, row_sd, "/")
  } else {
    print("Warnings: no scaling.")
    return(pre_scale)
  }
  return(post_scale)
}

# maxmin_n ####
maxmin_n <- function(m, n, option="max") {
  mlt = ifelse(option == "max", 1, -1)
  res <- order(-m * mlt)[seq_len(n)]
  list(values = m[res],
       position = res)
}
# correlated_peaks ####
correlated_peaks = function(mdata = raw_ls[[1]],
                            target_id = 75,
                            top_n = 3,
                            impute_options = c("threshold", "percentile"),
                            impute_options2 = c(T,F), # fill with random(T) or fixed(F)
                            normalize_options = c("row_mean"),
                            transform_options = c("log10",""),
                            scale_options = c("mean_center"),
                            limit_to_library = T,
                            print_pdf = F
)
{
  mdata_pre_clean = mdata[,14:(ncol(mdata)-2)]
  topn_list = list()
  bottomn_list = list()
  
  for(impute_method in impute_options){
    for(impute_method2 in impute_options2){
      for(normalize_method in normalize_options){
        for(transform_method in transform_options){
          for(scale_method in scale_options){
            mdata_clean = mdata_pre_clean %>% 
              data_impute(impute_method = impute_method, random = impute_method2) %>%
              data_normalize(nor_method = normalize_method) %>%
              data_transform(transform_method = transform_method) %>%
              data_scale(scale_method = scale_method)
            target = mdata_clean[target_id,]
            if(limit_to_library){
              mdata_clean = mdata_clean[!is.na(mdata$library_match_name),]
            }
            target_cor = cor(t(mdata_clean), t(target))
            topn_list[[length(topn_list)+1]]=maxmin_n(target_cor, top_n, "max")
            bottomn_list[[length(bottomn_list)+1]] = maxmin_n(target_cor, top_n, "min")
          }
        }
      }
    }
  }
  if(limit_to_library){
    mdata_in_library = mdata[!is.na(mdata$library_match_name),]
  } else {
    mdata_in_library = mdata
  }
  
  
  topn_summary = bind_rows(topn_list)
  topn_summary = topn_summary[with(topn_summary, order(-values)),]
  topn_summary[["ID"]] = row.names(mdata_in_library)[topn_summary$position]
  
  bottomn_summary = bind_rows(bottomn_list)
  bottomn_summary = bottomn_summary[with(bottomn_summary, order(values)),]
  bottomn_summary[["ID"]] = row.names(mdata_in_library)[bottomn_summary$position]
  
  data_select_topn = mdata_in_library[unique(topn_summary$position),]
  data_select_bottomn = mdata_in_library[unique(bottomn_summary$position),]
  
  figure_top = plot_library_bar(data_select_topn)
  figure_bottom = plot_library_bar(data_select_bottomn)
  
  if(print_pdf){
    pdf("correlated_peaks.pdf",onefile = T, width = 20, height = 10)
    print(figure_top)
    print(figure_bottom)
    dev.off()
  }
  return(list(top = list(values = topn_summary,
                         data_select = data_select_topn,
                         figure = figure_top),
              bottom = list(values = bottomn_summary,
                            data_select = data_select_bottomn,
                            figure = figure_bottom)
  ))
}

## my_plot_heatmap ####
my_plot_heatmap = function(raw_data,
                           cohort,
                           imgName = "Test", 
                           format = "pdf", # pdf, png, or jpeg
                           dpi = 72,
                           width = NA, # define output graph width
                           palette = "RdBu",  # RdBu, gbr, heat, topo
                           viewOpt = "overview", # detail, overview
                           rowV = T, colV = T, # cluster by row/column
                           border = T, # border for each pixel
                           grp.ave = F, # group average
                           scale_ub = 3, scale_lb = -3 # heatmap scale
)
{
  
  # imgName = "Test" 
  # format = "png" # pdf
  # dpi = 72
  # width = NA # define output graph width
  # palette = "RdBu"  # RdBu gbr heat topo
  # viewOpt = "overview" # Detail
  # rowV = T 
  # colV = T # cluster by row/column
  # border = T # border for each pixel
  # grp.ave = F # group average
  # scale_ub = 3 
  # scale_lb = -3 # heatmap scale
  # 
  
 
  imgName = paste0(imgName, "_dpi", dpi)
  if(file.exists(paste0(imgName, ".", format))){
    imgName = paste0(imgName, timestamp(), ".", format)
  }else{
    imgName = paste0(imgName, ".", format)
  }
     
     
     
     
  hc.dat <- as.matrix(raw_data)
  colnames(hc.dat) <- substr(colnames(hc.dat), 1, 32)
  if(!is.factor(cohort)){cohort = as.factor(cohort)}
  hc.cls <- cohort
  
  if (grp.ave) {
    lvs <- levels(hc.cls)
    my.mns <- matrix(ncol = ncol(hc.dat), nrow = length(lvs))
    for (i in 1:length(lvs)) {
      inx <- hc.cls == lvs[i]
      my.mns[i, ] <- apply(hc.dat[,inx], 1, mean)
    }
    rownames(my.mns) <- lvs
    colnames(my.mns) <- colnames(hc.dat)
    hc.dat <- my.mns
    hc.cls <- as.factor(lvs)
  }
  if (palette == "gbr") {
    colors <- colorRampPalette(c("green", "black", "red"),
                               space = "rgb")(256)
  } else if (palette == "heat") {
    colors <- heat.colors(256)
  } else if (palette == "topo") {
    colors <- topo.colors(256)
  } else if (palette == "gray") {
    colors <- colorRampPalette(c("grey90", "grey10"), space = "rgb")(256)
  } else {
    colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, 
                                                            "RdBu"))(256))
  }
  
  if (is.na(width)) {
    minW <- 630
    myW <- ncol(hc.dat) * 18 + 150
    if (myW < minW) {
      myW <- minW
    }
    w <- round(myW/72, 2)
  } else if (width == 0) {
    w <- 7.2
  } else {
    w <- 7.2
  }
  
  myH <- nrow(hc.dat) * 18 + 150
  h <- round(myH/72, 2)
  if (viewOpt == "overview") {
    if (is.na(width)) {
      if (w > 9) {
        w <- 9
      }
    }
    else if (width == 0) {
      if (w > 7.2) {
        w <- 7.2
      }
    }
    else {
      w <- 7.2
    }
    if (h > w) {
      h <- w
    }
  }
  if (grp.ave) {
    w <- ncol(hc.dat) * 25 + 300
    w <- round(w/72, 2)
  }
  if (border) {
    border.col <- "grey60"
  } else {
    border.col <- NA
  }
  if (format == "pdf") {
    pdf(file = imgName, width = w, height = h, bg = "white", 
        onefile = FALSE)
  } else {
    # will error if h is too large (due to large row number)
    Cairo::Cairo(file = imgName, 
                 unit = "in",
                 dpi = dpi,
                 width = w, height = h,
                 type = format, bg = "white"
    )
  }
  
  annotation <- data.frame(class = hc.cls)
  rownames(annotation) <- colnames(hc.dat)
  
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
  cols = dist.cols[as.numeric(hc.cls)]
  uniq.cols <- unique(cols)
  cls <- hc.cls
  # names(uniq.cols) <- unique(as.character((cls)))
  names(uniq.cols) <- unique(as.character(sort(cls)))
  ann_colors <- list(class = uniq.cols)
  
  if(scale_ub ==  scale_lb){breaksList = NA}
  breaksList = seq(scale_lb, scale_ub, by =(scale_ub-scale_lb)/length(colors))
  pheatmap::pheatmap(hc.dat, 
                     annotation = annotation, 
                     fontsize = 8, 
                     fontsize_row = 8, 
                     # clustering_distance_rows = smplDist, 
                     # clustering_distance_cols = smplDist, 
                     # clustering_method = clstDist, 
                     border_color = border.col, 
                     cluster_rows = rowV, 
                     cluster_cols = colV, 
                     # scale = scaleOpt, 
                     color = colors, 
                     annotation_colors = ann_colors,
                     breaks = breaksList)
  
  dev.off()
}
# cluster_mat from pheatmap
cluster_mat = function(mat, distance, method){
  if(!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
    stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
  }
  if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & class(distance) != "dist"){
    stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
  }
  if(distance[1] == "correlation"){
    d = as.dist(1 - cor(t(mat)))
  }
  else{
    if(class(distance) == "dist"){
      d = distance
    }
    else{
      d = dist(mat, method = distance)
    }
  }
  return(hclust(d, method = method))
}
## Cluster dataset to see if LC condition are similar ####
DatasetDist = function(dt1, dt2, 
                       log10_inten_cutoff = 4,  # only compare peaks with intensity above the threshold
                       merge_group_ppm_tol = 3 # mz within ppm_tol will merge into one mz group
)
{
  # dt1 = raw_ls[[1]] # liver in QE+
  # dt2 = raw_ls[[2]] # liver QE2
  # dt3 = raw_ls[[12]] # quad QE2
  
  if(identical(dt1, dt2)){return(0)}
  
  
  dt1_clean = dt1 %>%
    dplyr::select(ID, medMz, medRt, log10_inten) %>%
    filter(!duplicated(ID)) %>%
    filter(log10_inten > log10_inten_cutoff) %>%
    arrange(medMz)
  
  
  dt2_clean = dt2 %>%
    dplyr::select(ID, medMz, medRt, log10_inten) %>%
    filter(!duplicated(ID)) %>%
    filter(log10_inten > log10_inten_cutoff) %>%
    arrange(medMz)
  
  
  dt1_clean_medMz = dt1_clean$medMz
  nrow_dt1_clean = nrow(dt1_clean)
  
  
  i = 1
  count = 1
  dt1_clean_ls = list()
  while(i<=nrow_dt1_clean){
    temp_mz = dt1_clean_medMz[i]
    i_max = i
    while(i_max < nrow_dt1_clean){
      if(dt1_clean_medMz[i_max+1] - dt1_clean_medMz[i_max] <= dt1_clean_medMz[i_max] * merge_group_ppm_tol/1E6){
        i_max = i_max+1
      } else {
        break
      }
    }
    temp_dt1 = dt1_clean[i:i_max,]
    temp_df2_filter = dt2_clean$medMz>dt1_clean_medMz[i]*(1-merge_group_ppm_tol/1E6) & 
      dt2_clean$medMz<dt1_clean_medMz[i_max]*(1+merge_group_ppm_tol/1E6)
    temp_dt2 = dt2_clean[temp_df2_filter,]
    dt1_clean_ls[[count]] = list(dt1 = temp_dt1, dt2 = temp_dt2)
    i = i_max+1
    count = count+1
  }
  
  library(rdist)
  library(clue)
  test = unlist(lapply(dt1_clean_ls, function(dt_clean){
    # temp_dt1 = as.matrix(dt_clean$dt1[,4:5])
    # temp_dt2 = as.matrix(dt_clean$dt2[,4:5])
    temp_dt1 = as.matrix(dt_clean$dt1$medRt)
    temp_dt2 = as.matrix(dt_clean$dt2$medRt)
    if(nrow(temp_dt2)==0){return(NA)}
    
    x=cdist(temp_dt1, temp_dt2)
    if(ncol(x)<nrow(x)){x = t(x)}
    y=solve_LSAP(x)
    all_dist = x[cbind(seq_along(y), y)]
    # all_dist = 1
    return(median(all_dist))
  }))
  
  return(median(test, na.rm = T))
  
}
# legacy_code ####
{
  # Test row_mean and col_median effect ####
  # {  test = rbind(norm_row_mean_col_median=(norm_row_mean_col_median), norm_col_mean_row_mean=(norm_col_mean_row_mean))
  #   test = as.data.frame(t(test))
  #   test["cohort"] = stri_replace_last_regex(rownames(test),'_\\d+|-\\d+', '',stri_opts_regex(case_insensitive=T))
  #   unique_cohort = test$cohort
  #   if(length(unique(unique_cohort)) < 9){
  #     pal9 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  #              "#FFFF33", "#A65628", "#F781BF", "#999999")
  #     dist.cols = pal9[1:length(unique(unique_cohort))]
  #   } else {
  #     pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  #               "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  #               "#FFFF99", "#B15928")
  #     dist.cols <- colorRampPalette(pal12)(length(unique(unique_cohort)))
  #   }
  #   
  #   figure <- ggplot(test, aes(x=norm_row_mean_col_median, y=norm_col_mean_row_mean, color=cohort)) +
  #     geom_point() +
  #     theme(legend.title = element_blank()) +
  #     theme(legend.position = "right") +
  #     geom_abline(intercept = 0, slope = 1, alpha=0.5) +
  #     scale_color_manual(values = dist.cols) +
  #     scale_y_continuous(limits = c(0,1.2)) +
  #     scale_x_continuous(limits = c(0,1.2)) 
  #   # scale_x_continuous(limits = (0:2), breaks=seq(0, 2, .5))
  # }
  # my_cosine ####
  # my_cosine = function(x,y){
  #   if(length(x)!=length(y) | length(x) == 0){return (NULL)}
  #   x2s = sum(x * x)
  #   y2s = sum(y * y)
  #   xys = sum(x * y)
  #   cosine_similarity = xys / (sqrt(x2s)*sqrt(y2s))
  #   return(cosine_similarity)
  # }
  # cor(temp_x, temp_y, method = "pearson")
  # cor(temp_x, temp_y, method = "kendall")
  # cor(temp_x, temp_y, method = "spearman")
  # my_cosine(temp_x, temp_y)
  # Print_out intermediate number for data similarity search ####
  # topn_summary[["ID"]] = row.names(mdata_in_library)[topn_summary$position]
  # table_select = table(topn_summary$ID)
  # print(paste(impute_method,normalize_method,transform_method,scale_method, sep = "_"))
  # pstn = maxmin_n(target_cor, 12, 1)$position
  # print(pstn)
  # row.names(mdata_in_library)[pstn]
  # print(sum(table_select[row.names(mdata_in_library)[pstn]], na.rm = T))
  # data_select_topn["frequency"] = table_select[row.names(data_select_topn)]
  # figure_ls[[length(figure_ls)+1]] = figure
}