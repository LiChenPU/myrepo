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
  formula = "C5H11N1O1"
  
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



lb=1000 # min abs value of flux to be included in the plot
dblog_trans <- function(){
  trans_new(name='dblog', transform = function(x) (log10(abs(x)+lb)-log10(lb))*sign(x),
            inverse = function(x) sign(x)*(10^(abs(x)+log10(lb))-lb))
}
plot_bar_scatter_graph = function(plot_data, inde = "Peak", bar_plot = F)
{
  oldw <- getOption("warn")
  options(warn = -1)
  plot_data = unite(plot_data, "mz_RT", "medMz", "medRt", sep = "@")
  
  plot_upper_limit = 10^7
  if(log10(max(plot_data$Ion_count, na.rm=T))>7){
    plot_upper_limit = 10^ceiling(log10(max(plot_data$Ion_count, na.rm=T)))
  }
  
  
  data_df = plot_data
  data_df_spread = spread(data_df, tissue, Ion_count, drop=T)
  data_df_spread = data_df_spread[, -3]
  data_df_old = data_df_spread[data_df_spread$age=="Old",]
  data_df_young = data_df_spread[data_df_spread$age=="Young",]
  
  scatter_data = as.data.frame(matrix(ncol = 8, nrow=length(unique(data_df$tissue))))
  colnames(scatter_data) = c("id", "tissue", "old", "old_lb", "old_ub", "young", "young_lb", "young_ub")
  scatter_data$id = data_df$mz_RT[1]
  scatter_data$tissue = unique(data_df$tissue)
  scatter_data$old = colMeans(data_df_old[,-c(1,2)], na.rm=T)
  old_var = sapply(data_df_old[,-c(1,2)], sd, T)
  scatter_data$old_lb = scatter_data$old - old_var
  scatter_data$old_ub = scatter_data$old + old_var
  scatter_data$young = colMeans(data_df_young[,-c(1,2)], na.rm=T)
  young_var = sapply(data_df_young[,-c(1,2)], sd, T)
  scatter_data$young_lb = scatter_data$young - young_var
  scatter_data$young_ub = scatter_data$young + young_var
  
  scatter_data$young_lb[scatter_data$young_lb<0]=0
  scatter_data$old_lb[scatter_data$old_lb<0]=0
  
  fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e\\+", "10^", l)
    # return this as an expression
    parse(text=l)
  }
  

    
  bar_data = scatter_data
  bar_data_gather = gather(bar_data, key = statistic_label, value = number, -id, -tissue)
  bar_data_separate = separate(bar_data_gather, col = "statistic_label", into = c("age", "statistic"), sep = "_" )
  bar_data_separate$statistic[is.na(bar_data_separate$statistic)] = "TIC"
  bar_data = spread(bar_data_separate, key = "statistic", value = "number")
  
  
  figure_bar <- ggplot(bar_data, aes(x = tissue, y = TIC, fill = age))+
    geom_bar(position=position_dodge(),stat='identity', color = "black") +
    geom_errorbar(aes(ymin=lb, ymax=ub),width = 0.5, position=position_dodge(.9)) +
    theme(legend.title=element_blank()) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    
    # facet_wrap(~mz_RT, scales = "free")+
    # geom_text_repel(aes(label=tissue),colour='red') +
    # scale_y_continuous(limits = c(-0.2,1.2),breaks = c(0,0.5,1)) +
    scale_y_continuous(trans = 'dblog',limit=c(0,plot_upper_limit), 
                       breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9), 
                       label = fancy_scientific)+
    scale_fill_brewer(palette = 'Set3', guide = guide_legend(reverse = TRUE))
  figure_bar
  options(warn = oldw)
  return(figure_bar)
}
















