# 显示中文
# Sys.setlocale(category = "LC_ALL", locale = "Chinese")
# !diagnostics off

# Import library ####
{
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
# plot_library_bar ####
plot_library_bar = function(data_select){
  if(nrow(data_select) ==0){return(NULL)}
  # data_select = data_select_topn
  data_plot = data_select[,-c(2,4,5,6,7,8,9,12,13)]
  data_plot = data_plot[,-which(grepl("_inten",colnames(data_plot)))]
  data_gather = gather(data_plot, key = names, value = number, -formula, -label, -medMz, -medRt)
  data_gather["cohort"] = stri_replace_last_regex(data_gather$names,'\\d+|_\\d+|-\\d+|-a|-b|-c|_mean', '',stri_opts_regex(case_insensitive=T))
  data_gather = data_gather[complete.cases(data_gather),]
  
  unique_cohort = unique(data_gather$cohort)
  data_ls = list(length = length(unique_cohort))
  select_col = c("label","formula","medMz","medRt","cohort","TIC","sd")
  for(i in 1:length(unique_cohort)){
    temp = data_gather[data_gather$cohort==unique_cohort[i],]
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
                       formula = "C5H9N1O4")
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

# Substitute low signal and NA into small number ####
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
                          nor_method = "row_median", 
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
    return(pre_norm)
  } else {
    # print("incorrect input, nor_method options include row_median, row_mean, sample, cohort")
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
    return(pre_trsf)
  }
  return(post_trsf)
}


# Scaling ####
data_scale = function(pre_scale, scale_method = "mean_center"){
  row_mean = rowMeans(pre_scale)
  
  if(scale_method == "mean_center"){
    post_scale = sweep(pre_scale, 1, row_mean, "-")
  } else if(scale_method == "mean_center_then_sd"){
    post_scale = sweep(pre_scale, 1, row_mean, "-")
    row_sd = apply(pre_scale, 1, sd, na.rm = T)
    post_scale = sweep(pre_scale, 1, row_sd, "/")
  } else {
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
                            impute_options = c("threshold", "percentile"),
                            impute_options2 = c(T,F), # fill with random(T) or fixed(F)
                            normalize_options = c("row_mean"),
                            transform_options = c("log10",""),
                            scale_options = c("mean_center"),
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
            mdata_clean = mdata_clean[!is.na(mdata$library_match_name),]
            target_cor = cor(t(mdata_clean), t(target))
            topn_list[[length(topn_list)+1]]=maxmin_n(target_cor, 12, "max")
            bottomn_list[[length(bottomn_list)+1]] = maxmin_n(target_cor, 12, "min")
          }
        }
      }
    }
  }
  
  mdata_in_library = mdata[!is.na(mdata$library_match_name),]
  
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
# Main ####

# filenames = list.files(recursive = T)[grepl(".csv", list.files(recursive = T)) &
#                                       grepl("mdata", list.files(recursive = T))]

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
data_select_ls = lapply(raw_ls, filter_data, medMz = 0, formula = "C5H9N1O4")
lapply(data_select_ls, View)
lapply(data_select_ls, row.names)
fig_ls = lapply(data_select_ls, plot_library_bar)
print(fig_ls)



tissues = correlated_peaks(mdata = raw_ls[[1]],
                        target_id = 4971,
                        impute_options = c("threshold", "percentile"),
                        impute_options2 = c(T,F), # fill with random(T) or fixed(F)
                        normalize_options = c("row_mean"),
                        transform_options = c("log10",""),
                        scale_options = c("mean_center"),
                        print_pdf = F
                        )


