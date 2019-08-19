library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(scales)
# install.packages("ggrepel")
library(ggrepel)
# install.packages("gridExtra")
library(gridExtra)
library(stringi)
library(MetaboAnalystR)
library(pheatmap)
library(dplyr)
library(profvis)
library(RColorBrewer)
library(janitor)

# function ####
## Graphic support ####
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
lb=1000 # min value to be included in the plot
# double log transform that enables plot both positive and negative value on a log scale 
dblog_trans <- function(){
  scales::trans_new(name='dblog', transform = function(x) (log10(abs(x)+lb)-log10(lb))*sign(x),
                    inverse = function(x) sign(x)*(10^(abs(x)+log10(lb))-lb))
}

## ranking ####
whichpartrev <- function(x, n=30) {
  which(x >= -sort(-x, partial=n)[n])
}
## scatter_plot ####
scatter_plot = function(plot_data, #expect 1 row of intensity and each column is a sample
                        cohortTable, #expect dataframe of cohort, each column is a cohort, at least 2 column
                        cohort1 = 1, #first column in cohortTable, each one is a dot in scatter
                        cohort2 = 2, #second column in cohortTable, used in x or y axis
                        xaxis = "WT", yaxis = "KO", #select which factor as x or y axis
                        fig_title = ""
)
{
  data_df = plot_data %>% gather(na.rm = T)
  data_df["cohort1"] = cohortTable[data_df$key, cohort1]
  data_df["cohort2"] = cohortTable[data_df$key, cohort2]
  data_df = data_df %>% filter(cohort1!="blank")
  cohort1Level = unique(data_df$cohort1)
  cohort2Level = unique(data_df$cohort2)
  mean_matrix = sd_matrix = matrix(0, 
                                   nrow=length(cohort1Level), 
                                   ncol=length(cohort2Level))
  for(i in 1:length(cohort1Level)){
    for(j in 1:length(cohort2Level)){
      mean_matrix[i,j] = mean(data_df$value[data_df$cohort1== cohort1Level[i] & data_df$cohort2== cohort2Level[j]], na.rm=T)
      sd_matrix[i,j] = sd(data_df$value[data_df$cohort1== cohort1Level[i] & data_df$cohort2== cohort2Level[j]], na.rm=T)
    }
  }
  
  scatter_data=data.frame(mean_matrix, (mean_matrix-sd_matrix),(mean_matrix+sd_matrix))
  colnames(scatter_data) = c(cohort2Level, paste(cohort2Level, "lb", sep="_"), paste(cohort2Level, "ub", sep="_"))
  rownames(scatter_data) = scatter_data["cohort1"] = cohort1Level
  

  plot_upper_limit = 10^7
  if(log10(max(plot_data, na.rm=T))>7){
    plot_upper_limit = 10^ceiling(log10(max(plot_data, na.rm=T)))
  }
  colnames_scatter_data = colnames(scatter_data)
  names(colnames_scatter_data) = colnames_scatter_data
  
  xStringName = colnames_scatter_data[xaxis]
  yStringName = colnames_scatter_data[yaxis]
  
  
  figure <- ggplot(scatter_data, aes_string(x=xStringName, y=yStringName)) +
    geom_point(aes(color = cohort1)) +
    ggtitle(fig_title) +
    # theme(legend.title = element_blank()) + 
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, alpha=0.5) +
    geom_text_repel(aes(label=cohort1),colour='red') +
    geom_errorbarh(aes_string(xmin = paste0(xStringName, "_lb"), xmax = paste0(xStringName, "_ub")), height = 0, alpha=0.5) +
    geom_errorbar(aes_string(ymin = paste0(yStringName, "_lb"), ymax = paste0(yStringName, "_ub")), width = 0, alpha=0.5) +
    scale_x_continuous(trans = 'dblog',limit=c(-10,plot_upper_limit), 
                       breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9),
                       labels = fancy_scientific) +
    scale_y_continuous(trans = 'dblog',limit=c(-100,plot_upper_limit), breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9), 
                       labels = fancy_scientific)
  return(figure)
}

## bar_plot ####
bar_plot = function(plot_data, #expect 1 row of intensity and each column is a sample
                    cohortTable, #expect dataframe of cohort, each column is a cohort, at least 2 column
                    cohort1 = 1, #first column in cohortTable, each factor is a group of bars
                    cohort2 = 2, #second column in cohortTable, each factor is a bar in group
                    xaxis = "Tissues", yaxis = "TIC", #Set xaxis and yaxis label
                    fig_title = "",
                    log_scale = T
)
{
  data_df = plot_data %>% gather(na.rm = T)
  
  if(length(cohort1)>1){
    data_df["cohort1"] = apply(cohortTable[data_df$key, cohort1], 1, paste, collapse = ".")
  } else {
    data_df["cohort1"] = cohortTable[data_df$key, cohort1]
  }
  
  if(length(cohort2)>1){
    data_df["cohort2"] = apply(cohortTable[data_df$key, cohort2], 1, paste, collapse = ".")
  } else {
    data_df["cohort2"] = cohortTable[data_df$key, cohort2]
  }
  
  # data_df = data_df %>% filter(cohort1!="blank")
  cohort1Level = unique(data_df$cohort1)
  cohort2Level = unique(data_df$cohort2)
  mean_matrix = sd_matrix = matrix(0, 
                                   nrow=length(cohort1Level), 
                                   ncol=length(cohort2Level))
  for(i in 1:length(cohort1Level)){
    for(j in 1:length(cohort2Level)){
      mean_matrix[i,j] = mean(data_df$value[data_df$cohort1== cohort1Level[i] & data_df$cohort2== cohort2Level[j]], na.rm=T)
      sd_matrix[i,j] = sd(data_df$value[data_df$cohort1== cohort1Level[i] & data_df$cohort2== cohort2Level[j]], na.rm=T)
    }
  }
  
  bar_data=data.frame(mean_matrix, (mean_matrix-sd_matrix),(mean_matrix+sd_matrix))
  colnames(bar_data) = c(cohort2Level, paste(cohort2Level, "lb", sep="_"), paste(cohort2Level, "ub", sep="_"))
  rownames(bar_data) = bar_data["cohort1"] = cohort1Level
  
  bar_data = bar_data %>%
    gather(key = statistic_label, value = number, -cohort1) %>%
    separate(col = "statistic_label", into = c("cohort2", "statistic"), sep = "_", fill="right") %>%
    dplyr::mutate(statistic = replace_na(statistic, "TIC")) %>%
    spread(key = "statistic", value = "number")
  

  plot_upper_limit = 10^7
  if(log10(max(plot_data, na.rm=T))>7){
    plot_upper_limit = 10^ceiling(log10(max(plot_data, na.rm=T)))
  }
  colnames_bar_data = colnames(bar_data)
  names(colnames_bar_data) = colnames_bar_data
  
  
  # oldw <- getOption("warn")
  # options(warn = -1)
  # options(warn = oldw)

  
  figure_bar <- ggplot(bar_data, aes_string(x = "cohort1", y = "TIC", fill = "cohort2"))+
    geom_bar(position=position_dodge(),stat='identity', color = "black") +
    geom_errorbar(aes(ymin=lb, ymax=ub),width = 0.5, position=position_dodge(.9)) +
    ggtitle(fig_title) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.title = element_blank()) +
    labs(x=xaxis, y=yaxis) +
    scale_fill_brewer(palette = 'Set3', guide = guide_legend())
  if(log_scale){
    figure_bar = figure_bar +
      scale_y_continuous(trans = 'dblog',limit=c(0,plot_upper_limit), 
                         breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9), 
                         label = fancy_scientific)
  }

  return(figure_bar)
}

## my_plot_heatmap ####
my_plot_heatmap = function(raw_data,
                           cohort,
                           imgName = "Test", 
                           format = "pdf", # pdf, png, or jpeg
                           dpi = 72,
                           width = NA, # define output graph width
                           palette = "RdBu",  # RdBu, gbr, heat, topo
                           viewOpt = "overview", # Detail
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
  
  printtime = Sys.time()
  matches <- paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
  imgName = paste(imgName, "_","dpi", dpi,"_", matches, ".", format, sep = "")
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

# main ####
## Load data ####
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  source("../metabo_library_functions.R")
  setwd("../library/pig_blood_pos")
  load("merge_mdata.RData")
  
  summaryTable_ls = split(summaryTable, summaryTable$MZRT_group)
  
  # Intensity dataframe
  data_all = intenTable
  
  # Statistics dataframe
  data_statistics = merge_mdata[,1:which(colnames(merge_mdata)=="Sum_logP")]
  
  # cohortTable
  sample_names = colnames(intenTable)[-1]
  blank_names = sample_names[grepl("blank|blk", sample_names, ignore.case = T)]
  
  test = str_split(sample_names, "_")
  n.obs <- sapply(test, length)
  seq.max <- seq_len(max(n.obs))
  mat <- t(sapply(test, "[", i = seq.max))
  
  cohortTable = as.data.frame(mat, stringsAsFactors =F)
  rownames(cohortTable) = sample_names
  cohortTable[grepl("blank|blk", rownames(cohortTable), ignore.case = T),] = "blank"
}

## Sample plot ####
{
  medMz = 545.9192
  medRt = 5.128
  
  
 
  targetMZRTGroupID = unique(filter_data(summaryTable, medMz = medMz, formula="")$MZRT_group)
  
  plot_data = intenTable[intenTable$MZRT_group==targetMZRTGroupID,-1]
  metaData_select = data_statistics[data_statistics$MZRT_group==targetMZRTGroupID,]
  
  ## Scatter plot
  temp_mz = round(metaData_select$medMz,4)
  temp_rt = metaData_select$medRt
  potential_formula = paste(unique(as.character(metaData_select[1,grepl("formula", colnames(metaData_select))])), collapse = "/")
  sample_scatter_plot = scatter_plot(plot_data = plot_data,
               cohortTable = cohortTable,
               cohort1 = 1, #first column in cohortTable, each one is a dot in scatter
               cohort2 = 4, #second column in cohortTable, used in x or y axis
               xaxis = "draw1", yaxis = "draw2", #select which factor as x or y axis
               fig_title = paste(temp_mz, temp_rt, potential_formula, sep="_")
  )
  
  ## Bar plot
  sample_bar_plot = bar_plot(plot_data = plot_data, #expect 1 row of intensity and each column is a sample
           cohortTable = cohortTable, #expect dataframe of cohort, each column is a cohort, at least 2 column
           cohort1 = 1, #first column in cohortTable, each factor is a group of bars
           cohort2 = c(3,4), #second column in cohortTable, each factor is a bar in group
           xaxis = "Tissues", yaxis = "TIC", #Set xaxis and yaxis label
           log_scale = F,
           fig_title = paste(temp_mz, temp_rt, potential_formula, sep="_")
  )
  
  ## heatmap
  cohort = paste(cohortTable$V1, cohortTable$V2, sep="_")
  sample_data = intenTable[1:200,-1]
  sample_data_normalize = sample_data %>%
    data_impute(impute_method = "threshold", random = F) %>%
    data_normalize(nor_method = "row_mean") %>%
    data_transform(transform_method = "log10") %>%
    data_scale(scale_method = "")
  
  my_plot_heatmap(sample_data_normalize, 
                  cohort,
                  format = "pdf")
  
}


## top10 sum logP ####
{
  # select
  topn_sum_logP = 10
  sum_logP = sapply(summaryTable_ls, function(x){
    x = x[!duplicated(x$ID),]
    sum(x$`_log10_FDR`, na.rm = T)
  })
  id_sum_logP = whichpartrev(sum_logP, topn_sum_logP)
  id_sum_logP = id_sum_logP[order(-sum_logP[id_sum_logP])]
  
  # plot
  data_select = intenTable[id_sum_logP, -1]
  metaData_select = data_statistics[id_sum_logP,]
  figure_list = list()
  for(i in 1:nrow(data_select)){
    plot_data = data_select[i,]
    temp_mz = round(metaData_select$medMz[i],4)
    temp_rt = metaData_select$medRt[i]
    potential_formula = paste(unique(as.character(metaData_select[i,grepl("formula", colnames(metaData_select))])), collapse = "/")
    # figure_list[[i]] = scatter_plot(plot_data = plot_data,
    #                                 cohortTable = cohortTable,
    #                                 cohort1 = 1, #first column in cohortTable, each one is a dot in scatter
    #                                 cohort2 = 4, #second column in cohortTable, used in x or y axis
    #                                 xaxis = "draw1", yaxis = "draw2", #select which factor as x or y axis
    #                                 fig_title = paste(temp_mz, temp_rt, potential_formula, sep="_"))
    # 
    ## Bar plot
    figure_list[[i]] = bar_plot(plot_data = plot_data, #expect 1 row of intensity and each column is a sample
                               cohortTable = cohortTable, #expect dataframe of cohort, each column is a cohort, at least 2 column
                               cohort1 = c(3,4), #first column in cohortTable, each factor is a group of bars
                               cohort2 = 1, #second column in cohortTable, each factor is a bar in group
                               xaxis = "pig_draw", yaxis = "TIC", #Set xaxis and yaxis label
                               log_scale = F,
                               fig_title = paste(temp_mz, temp_rt, potential_formula, sep="_")
    )
  }
  
  # print
  ml = marrangeGrob(figure_list, nrow=4, ncol=1, 
                    top = paste("Top", topn_sum_logP, "significant peaks."))
  pdf(paste("topn_sum_logP.pdf"), onefile = TRUE, width=10, height = 15)
  print(ml)
  dev.off()
  
}

## peaks with n signif ####
{
  # select
  signif_n = 3
  signif_level = 0.2
  signif_topn = 20
  stat_signif_num = sapply(summaryTable_ls, function(x){
    x = x[!duplicated(x$ID),]
    sum(x$`_log10_FDR` > -log10(signif_level),na.rm=T)
  })
  id_signif_n = stat_signif_num[stat_signif_num>signif_n]
  id_signif_n = id_signif_n[order(-id_signif_n)]
  id_signif_n = id_signif_n[1:min(length(id_signif_n), signif_topn)]
  id_signif_n = as.numeric(names(id_signif_n))
  
  # plot
  data_select = intenTable[id_signif_n, -1]
  metaData_select = data_statistics[id_signif_n,]
  figure_list = list()
  for(i in 1:nrow(data_select)){
    plot_data = data_select[i,]
    temp_mz = round(metaData_select$medMz[i],4)
    temp_rt = metaData_select$medRt[i]
    potential_formula = paste(unique(as.character(metaData_select[i,grepl("formula", colnames(metaData_select))])), collapse = "/")
    figure_list[[i]] = scatter_plot(plot_data = plot_data,
                                    cohortTable = cohortTable,
                                    cohort1 = 1, #first column in cohortTable, each one is a dot in scatter
                                    cohort2 = 2, #second column in cohortTable, used in x or y axis
                                    xaxis = "WT", yaxis = "KO", #select which factor as x or y axis
                                    fig_title = paste(temp_mz, temp_rt, potential_formula, sep="_")
                                    )
  }

  # print
  ml = marrangeGrob(figure_list, nrow=3, ncol=2, 
                    top = paste("Peaks significant in more than", signif_n, "different cohorts"))
  pdf(paste("signif_counts_plots.pdf"), onefile = TRUE, width=10, height = 15)
  print(ml)
  dev.off()
}

## top n signif in each tissue ####
{
  # select
  tissue_topn = 10
  stat_signif_num = list()
  i=1
  for(i in 1:length(unique(summaryTable$filename))){
    id_tissue_topn = summaryTable %>%
      distinct(ID, .keep_all=T) %>%
      filter(filename == unique(summaryTable$filename)[i]) %>%
      arrange(desc(`_log10_FDR`)) %>%
      top_n(tissue_topn,`_log10_FDR`)
    stat_signif_num[[i]] = id_tissue_topn$MZRT_group
  }
  
  
  # plot
  ml_list = list()
  for(j in 1:length(unique(summaryTable$filename))){
    data_select = intenTable[stat_signif_num[[j]], -1]
    metaData_select = data_statistics[stat_signif_num[[j]],]
    figure_list = list()
    for(i in 1:nrow(data_select)){
      plot_data = data_select[i,]
      temp_mz = round(metaData_select$medMz[i],4)
      temp_rt = metaData_select$medRt[i]
      potential_formula = paste(unique(as.character(metaData_select[i,grepl("formula", colnames(metaData_select))])), collapse = "/")
      figure_list[[i]] = scatter_plot(plot_data = plot_data,
                                      cohortTable = cohortTable,
                                      cohort1 = 1, #first column in cohortTable, each one is a dot in scatter
                                      cohort2 = 2, #second column in cohortTable, used in x or y axis
                                      xaxis = "WT", yaxis = "KO", #select which factor as x or y axis
                                      fig_title = paste(temp_mz, temp_rt, potential_formula, sep="_")
      )
    }
    ml = marrangeGrob(figure_list, nrow=3, ncol=2, 
                      top = paste("Top",tissue_topn, "peaks significant in", unique(summaryTable$filename)[j]))
    ml_list[[j]] = ml
  }
  
  # print
 
  pdf(paste("top_n_in_tissue"), onefile = TRUE, width=10, height = 15)
  print(ml_list)
  dev.off()
}

## heatmap ####
{
  
  cohort = paste(cohortTable$V1, cohortTable$V2, sep="_")
  raw_data = intenTable[1:200,-1]
  raw_data_normalize = raw_data %>%
    data_impute(impute_method = "threshold", random = F) %>%
    data_normalize(nor_method = "row_mean") %>%
    data_transform(transform_method = "log10") %>%
    data_scale(scale_method = "")
  
  my_plot_heatmap(raw_data_normalize, 
                  cohort,
                  format = "pdf")
  
  # sum_logP_plot_data = intenTable[id_sum_logP,]
  # plot_data_note = 
  # 
  # signif_n_plot_data = intenTable[id_signif_n, -1]
  # signif_n_plot_data["Note"] = paste0("signif_n=",stat_signif_num[id_signif_n])
}




  
  
## Volcano plot ####
  

mSet<-Volcano.Anal(mSet, paired = F, fcthresh = 2.0, cmpType = 1, percent.thresh = 0.75,
                   nonpar = F, threshp = 0.1, equal.var	= TRUE, pval.type = "raw")
mSet<-PlotVolcano(mSet, "volcano_0_",1, "png", 72, width=NA)
test = mSet$analSet$volcano$p.log

function (mSetObj = NA, paired = FALSE, fcthresh, cmpType, percent.thresh, 
          nonpar = F, threshp, equal.var = TRUE, pval.type = "raw") 
{
  paired = F 
  fcthresh = 2.0 
  cmpType = 0 
  percent.thresh = 0.75
  nonpar = F 
  threshp = 0.1 
  equal.var	= TRUE 
  pval.type = "raw"
  
  mSetObj <- mSet 
  p.value = 10^-(mSet$analSet$volcano$p.log)
    
  inx.p <- p.value <= threshp
  p.log <- -log10(p.value)
  fcthresh = ifelse(fcthresh > 1, fcthresh, 1/fcthresh)
  max.xthresh <- log(fcthresh, 2)
  min.xthresh <- log(1/fcthresh, 2)
  res <- GetFC(mSetObj, F, cmpType)
  fc.log <- res$fc.log
  fc.all <- res$fc.all
  inx.up <- fc.log > max.xthresh
  inx.down <- fc.log < min.xthresh
  # if (paired) {
  #   res <- GetFC(mSetObj, T, cmpType)
  #   count.thresh <- round(nrow(mSetObj$dataSet$norm)/2 * 
  #                           percent.thresh)
  #   mat.up <- res >= max.xthresh
  #   mat.down <- res <= min.xthresh
  #   count.up <- apply(mat.up, 2, sum)
  #   count.down <- apply(mat.down, 2, sum)
  #   fc.all <- rbind(count.up, count.down)
  #   inx.up <- count.up >= count.thresh
  #   inx.down <- count.down >= count.thresh
  #   colnames(fc.all) <- colnames(mSetObj$dataSet$norm)
  #   rownames(fc.all) <- c("Count (up)", "Count (down)")
  #   max.xthresh <- count.thresh
  #   min.xthresh <- -count.thresh
  # }
  inx.imp <- (inx.up | inx.down) & inx.p
  if (paired) {
    sig.var <- cbind(fc.all[1, ][inx.imp, drop = F], fc.all[2, 
                                                            ][inx.imp, drop = F], p.value[inx.imp, drop = F], 
                     p.log[inx.imp, drop = F])
    if (pval.type == "fdr") {
      colnames(sig.var) <- c("Counts (up)", "Counts (down)", 
                             "p.adjusted", "-log10(p)")
    }
    else {
      colnames(sig.var) <- c("Counts (up)", "Counts (down)", 
                             "raw.pval", "-log10(p)")
    }
    dif.count <- abs(sig.var[, 1] - sig.var[, 2])
    ord.inx <- order(dif.count, sig.var[, 4], decreasing = T)
    sig.var <- sig.var[ord.inx, , drop = F]
    sig.var[, c(3, 4)] <- signif(sig.var[, c(3, 4)], 5)
  }
  else {
    sig.var <- cbind(fc.all[inx.imp, drop = F], fc.log[inx.imp, 
                                                       drop = F], p.value[inx.imp, drop = F], p.log[inx.imp, 
                                                                                                    drop = F])
    if (pval.type == "fdr") {
      colnames(sig.var) <- c("FC", "log2(FC)", "p.ajusted", 
                             "-log10(p)")
    }
    else {
      colnames(sig.var) <- c("FC", "log2(FC)", "raw.pval", 
                             "-log10(p)")
    }
    ord.inx <- order(sig.var[, 4], abs(sig.var[, 2]), decreasing = T)
    sig.var <- sig.var[ord.inx, , drop = F]
    sig.var <- signif(sig.var, 5)
  }
  fileName <- "volcano.csv"
  write.csv(signif(sig.var, 5), file = fileName)
  volcano <- list(raw.threshx = fcthresh, raw.threshy = threshp, 
                  paired = paired, max.xthresh = max.xthresh, min.xthresh = min.xthresh, 
                  thresh.y = -log10(threshp), fc.all = fc.all, fc.log = fc.log, 
                  fc.log.uniq = jitter(fc.log), inx.up = inx.up, inx.down = inx.down, 
                  p.log = p.log, inx.p = inx.p, sig.mat = sig.var)
  mSetObj$analSet$volcano <- volcano
  return(.set.mSet(mSetObj))
}




{
  mSetObj = mSet
  imgName = "volcano_0_"
  plotLbl = 1
  format = "png"
  dpi = 72
  width = NA
  imgName = paste(imgName, "dpi", dpi, ".", format, sep = "")
  
  
  if (is.na(width)) {
    w <- 10
  }
  else if (width == 0) {
    w <- 8
  }
  else {
    w <- width
  }
  h <- w * 6/10
  mSetObj$imgSet$volcano <- imgName
  Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, 
               height = h, type = format, bg = "white")
  par(mar = c(5, 5, 3, 4))
  vcn <- mSetObj$analSet$volcano
  MyGray <- rgb(t(col2rgb("black")), alpha = 40, maxColorValue = 255)
  MyHighlight <- rgb(t(col2rgb("magenta")), alpha = 80, maxColorValue = 255)
  
  
  # if (vcn$paired) {
  #   xlim <- c(-nrow(mSetObj$dataSet$norm)/2, nrow(mSetObj$dataSet$norm)/2) * 
  #     1.2
  #   fc.all <- apply(vcn$fc.all, 2, function(x) {
  #     if (x[1] > x[2]) {
  #       return(x[1])
  #     }
  #     else {
  #       return(-x[2])
  #     }
  #   })
  #   hit.inx <- vcn$inx.p & (vcn$inx.up | vcn$inx.down)
  #   plot(fc.all, vcn$p.log, xlim = xlim, pch = 20, cex = ifelse(hit.inx, 
  #                                                               1.2, 0.8), col = ifelse(hit.inx, MyHighlight, MyGray), 
  #        xlab = "Count of Significant Pairs", ylab = "-log10(p)")
  #   sig.upInx <- vcn$inx.p & vcn$inx.up
  #   p.topInx <- GetTopInx(vcn$p.log, 5, T) & vcn$inx.up
  #   fc.rtInx <- GetTopInx(vcn$fc.all[1, ], 5, T)
  #   lblInx <- p.topInx & sig.upInx & fc.rtInx
  #   if (plotLbl & sum(lblInx, na.rm = T) > 0) {
  #     text.lbls <- substr(colnames(mSetObj$dataSet$norm)[lblInx], 
  #                         1, 14)
  #     text(vcn$fc.all[1, lblInx], vcn$p.log[lblInx], labels = text.lbls, 
  #          pos = 4, col = "blue", srt = 30, xpd = T, cex = 0.8)
  #   }
  #   sig.dnInx <- vcn$inx.p & vcn$inx.down
  #   p.topInx <- GetTopInx(vcn$p.log, 5, T) & vcn$inx.down
  #   fc.leftInx <- GetTopInx(vcn$fc.all[2, ], 5, T) & vcn$inx.down
  #   lblInx <- p.topInx & sig.dnInx & fc.leftInx
  #   if (plotLbl & sum(lblInx, na.rm = T) > 0) {
  #     text.lbls <- substr(colnames(mSetObj$dataSet$norm)[lblInx], 
  #                         1, 14)
  #     text(-vcn$fc.all[2, lblInx], vcn$p.log[lblInx], 
  #          labels = text.lbls, pos = 2, col = "blue", srt = -30, 
  #          xpd = T, cex = 0.8)
  #   }
  # }
  # else {
    imp.inx <- (vcn$inx.up | vcn$inx.down) & vcn$inx.p
    plot(vcn$fc.log, vcn$p.log, pch = 20, cex = ifelse(imp.inx, 
                                                       1.2, 0.7), col = ifelse(imp.inx, MyHighlight, MyGray), 
         xlab = "log2 (FC)", ylab = "-log10(p)")
    sig.inx <- imp.inx
    p.topInx <- GetTopInx(vcn$p.log, 5, T) & (vcn$inx.down)
    fc.leftInx <- GetTopInx(vcn$fc.log, 5, F)
    lblInx <- sig.inx & (p.topInx | fc.leftInx)
    lblInx = c(1,2,3,4,5)
    if (plotLbl & sum(lblInx, na.rm = T) > 0) {
      text.lbls <- substr(colnames(mSetObj$dataSet$norm)[lblInx], 
                          1, 14)
      text(vcn$fc.log[lblInx], vcn$p.log[lblInx], labels = text.lbls, 
           pos = 2, col = "blue", srt = -30, xpd = T, cex = 0.8)
    }
    p.topInx <- GetTopInx(vcn$p.log, 5, T) & (vcn$inx.up)
    fc.rtInx <- GetTopInx(vcn$fc.log, 5, T)
    lblInx <- sig.inx & (p.topInx | fc.rtInx)
    if (plotLbl & sum(lblInx, na.rm = T) > 0) {
      text.lbls <- substr(colnames(mSetObj$dataSet$norm)[lblInx], 
                          1, 14)
      text(vcn$fc.log[lblInx], vcn$p.log[lblInx], labels = text.lbls, 
           pos = 4, col = "blue", srt = 30, xpd = T, cex = 0.8)
    }
  # }
  abline(v = vcn$max.xthresh, lty = 3)
  abline(v = vcn$min.xthresh, lty = 3)
  abline(h = vcn$thresh.y, lty = 3)
  axis(4)
  dev.off()
  return(.set.mSet(mSetObj))
}




