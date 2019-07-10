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
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}


filenames = list.files(recursive = T)[grepl(".csv", list.files(recursive = T)) &
                                      grepl("mdata", list.files(recursive = T))]


num_of_files = length(filenames)
raw_ls = list()
i=1
for(i in 1:num_of_files){
  filename=filenames[i]
  df_temp = read_csv(paste("./",filename, sep=""))
  df_temp = cbind(label = gsub("/.*","",filename),df_temp, stringsAsFactors = F)
  raw_ls[[i]]= df_temp
}
rm(df_temp)

{
  raw = bind_rows(raw_ls)
  ncol_raw = ncol(raw)
  nrow_raw = nrow(raw)
  # object.size(raw)
}

{
  medMz = 0
  medRt = 0
  delta_mz = 0.001
  delta_rt = 0.1
  formula = "C5H9N1O4"
  
  data_select = raw
  if(medMz != 0){
    data_select = data_select[abs(data_select$medMz - medMz) < delta_mz,]
  }
  if(medRt != 0){
    data_select = data_select[abs(data_select$medRt - medRt) < delta_rt,]
  }
  if(formula != ""){
    data_select = data_select[data_select$formula == formula & !is.na(data_select$formula),]
  }
}


plot_library_bar = function(data_select){
  data_select$medMz = round(data_select$medMz, 4)
  data_select$medRt = round(data_select$medRt, 3)
  data_plot = data_select[,-c(2,4,5,6,7,8,9,12,13)]
  # data_inten = data_plot[,14:ncol(data_plot)]
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
  # data_bind = unite(data_bind, "label", c(label, cohort))
  data_bind = unite(data_bind, "title", c(formula, medMz, medRt))
  data_bind = data_bind[with(data_bind, order(-TIC)),]
  
  figure <- ggplot(data_bind, aes(x = cohort, y = TIC, fill = label))+
    geom_bar(position=position_dodge(),stat='identity', color = "black") +
    geom_errorbar(aes(ymin=TIC-sd, ymax=TIC+sd),width = 0.5, position=position_dodge(.9)) +
    theme(legend.title=element_blank()) +
    facet_wrap(~title, scales = "free")+
    # scale_y_continuous(limits = c(-0.2,1.2),breaks = c(0,0.5,1)) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    scale_fill_brewer(palette = 'Set3', guide = guide_legend(reverse = TRUE))
  
  return(figure)
}



figure = plot_library_bar(data_select)
pdf("library.pdf",onefile = TRUE)
print(figure)
dev.off()





