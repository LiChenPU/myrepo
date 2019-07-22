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

## Data_clean ####
data_clean = function(giant_df){
  giant_colname = colnames(giant_df)
  giant_colname = gsub("\\.", "_", giant_colname)
  giant_colname = gsub("old", "Old", giant_colname)
  giant_colname = gsub("young", "Young", giant_colname)
  colnames(giant_df) = giant_colname
  giant_df$medMz = round(giant_df$medMz, digit = 5)
  giant_df$medRt = round(giant_df$medRt, digit = 3)
  
  return(giant_df)
}
## Intensity dataframe ####
intensity_dataframe = function(giant_df)
{
  giant_colname = colnames(giant_df)
  data_raw = giant_df[,c(1,2,grep("[[:digit:]]",giant_colname))]
  data_gather <- gather(data_raw, key = tissue_age_id, value = Ion_count, -"medMz", -"medRt")
  data_separate = separate(data_gather, col = "tissue_age_id", into = c("tissue","age","mice_id"), remove = T)
  data_separate$age = as.factor(data_separate$age)
  data_separate$tissue = as.factor(data_separate$tissue)
  data_separate$mice_id = as.factor(data_separate$mice_id)
  data_all = data_separate
  return(data_all)
}

## normalization_each_tissue ####
normalization_each_tissue = function(data_all)
{
  data_metaboanalyst = data_all
  data_metaboanalyst$Ion_count[is.na(data_metaboanalyst$Ion_count)] = 0
  data_metaboanalyst$Ion_count[data_metaboanalyst$Ion_count==0] = runif(length(data_metaboanalyst$Ion_count[data_metaboanalyst$Ion_count==0]), 0, 2500)
  test = spread(data_metaboanalyst, key=mice_id, value=Ion_count)
  test_rowmean = rowMeans(test[,5:ncol(test)], na.rm = T)
  test_rowmean_matrix = matrix(test_rowmean, nrow=2, byrow=F)
  test_rowmean_matrix[1,] = test_rowmean_matrix[2,]
  test_rowmean2 = as.numeric(test_rowmean_matrix)
  test[,5:ncol(test)] =sweep(test[,5:ncol(test)],1,test_rowmean2,FUN = "/")
  
  data_metab_unite = unite(test, "tissue_age", tissue, age)
  data_metab_normalize = gather(data_metab_unite, key = mice_id, value = Norm_ion_count, -medMz, -medRt, -tissue_age, na.rm = T)
  data_metab_normalize = unite(data_metab_normalize, "tissue_age_id", tissue_age, mice_id)
  data_metab_spread = spread(data_metab_normalize, key = tissue_age_id, value = Norm_ion_count)
  return(data_metab_spread)
}
## plot_bar_scatter_graph ####
# double log transform that enables plot both positive and negative value on a log scale 
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
  
  # fancy_scientific <- function(l) {
  #   # turn in to character string in scientific notation
  #   l <- format(l, scientific = TRUE)
  #   # quote the part before the exponent to keep all the digits
  #   l <- gsub("^(.*)e", "'\\1'e", l)
  #   # turn the 'e+' into plotmath format
  #   l <- gsub("e", "%*%10^", l)
  #   # return this as an expression
  #   parse(text=l)
  # }
  
  
  ## Scatter plot
  {
    figure <- ggplot(scatter_data, aes(x=young, y=old, color = tissue)) +
      geom_point() +
      ggtitle(paste(inde,"_", plot_data$mz_RT[1], sep="")) +
      # theme(legend.title = element_blank()) + 
      theme(legend.position = "none") +
      geom_abline(intercept = 0, slope = 1, alpha=0.5) +
      geom_text_repel(aes(label=tissue),colour='red') +
      geom_errorbarh(aes(xmin = young_lb, xmax = young_ub), height = 0, alpha=0.5) +
      geom_errorbar(aes(ymin = old_lb, ymax = old_ub), width = 0, alpha=0.5) +
      scale_x_continuous(trans = 'dblog',limit=c(-10,plot_upper_limit), 
                         breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9),
                         labels = fancy_scientific) +
      scale_y_continuous(trans = 'dblog',limit=c(-100,plot_upper_limit), breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9), 
                         labels = fancy_scientific)
    figure
  }
  
  # Bar plot
  if(bar_plot)
  {
    
    
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
    return(list(scatter = figure, bar = figure_bar))
  }
  
  options(warn = oldw)
  return(figure)
  
  # pdf(paste(inde,"_", plot_data$mz_RT[1],".pdf", sep=""),width = 6.5, height = 5)
  # print(figure_bar)
  # print(figure)
  # print(figure_bar2)
  # dev.off()
}


## plot_scatter_by_tissue ####
plot_scatter_by_tissue = function(data_all, data_statistics, top_n = 9, nrow=3, ncol=3)
{
  tissue = levels(data_all$tissue)
  setwd(project_dir)
  
  ml_list = list()
  for(j in 1:length(tissue)){
    data_select = data_statistics[with(data_statistics, order(-eval(parse(text = paste(tissue[j],"_y",sep=""))))),]
    data_select = data_select[1:top_n,]
    data_select = data_select[!is.na(data_select[paste(tissue[j],"_y",sep="")]),]
    
    figure_list = list()
    for(i in 1:nrow(data_select)){
      medMz = data_select$medMz[i]
      medRt = data_select$medRt[i]
      plot_data = data_all[data_all$medMz == medMz & data_all$medRt == medRt,]
      figure_list[[length(figure_list)+1]] = plot_bar_scatter_graph(plot_data = plot_data, i)
      
    }
    
    ml <- marrangeGrob(figure_list, nrow=nrow, ncol=ncol, top = tissue[j])
    ml_list[[length(ml_list)+1]] = ml
    
  }
  
  pdf("tissue_plot.pdf",onefile = TRUE)
  print(ml_list)
  dev.off()
}

## plot_scatter_by_tissue_num ####
plot_scatter_by_tissue_num =function(data_all, data_statistics, signif_counts = 2, signif_level = 0.05)
{
  setwd(project_dir)
  
  tissue = levels(data_all$tissue)
  signif_level = -log10(signif_level)
  
  stat_col = grepl("_y", colnames(data_statistics))
  data_statistics["signif_counts"] = apply(data_statistics, 1, function(x) {
    x = as.numeric(x[stat_col])
    non_na = !is.na(x)
    signif = sum(x[non_na]>signif_level, na.rm = T)
    return(signif)
  })
  
  data_signif_counts = data_statistics[with(data_statistics, order(-signif_counts, -sum_log_p)),]
  data_signif_counts = data_signif_counts[data_signif_counts$signif_counts>=signif_counts,]
  
  ml_list = list()
  for(j in unique(data_signif_counts$signif_counts)){
    data_select = data_signif_counts[data_signif_counts$signif_counts == j,]
    figure_list = list()
    for(i in 1:nrow(data_select)){
      medMz = data_select$medMz[i]
      medRt = data_select$medRt[i]
      plot_data = data_all[data_all$medMz == medMz & data_all$medRt == medRt,]
      figure_list[[length(figure_list)+1]] = plot_bar_scatter_graph(plot_data = plot_data, i)
    }
    ml_list[[length(ml_list)+1]] <- marrangeGrob(figure_list, nrow=3, ncol=3, top = paste("Peaks significant in", j, "different tissues."))
  }
  pdf(paste("signif_counts_plots.pdf"),onefile = TRUE)
  print(ml_list)
  dev.off()
}


## my_plot_heatmap ####
my_plot_heatmap = function(raw_data,
                           cohort,
                           imgName = "Test", 
                           format = "png", # pdf
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
  printtime = Sys.time()
  matches <- paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
  imgName = paste(imgName, "_","dpi", dpi,"_", matches, ".", format, sep = "")
  hc.dat <- as.matrix(raw_data)
  colnames(hc.dat) <- substr(colnames(hc.dat), 1, 32)
  hc.cls <- cohort
  
  if (grp.ave) {
    lvs <- levels(hc.cls)
    my.mns <- matrix(ncol = ncol(hc.dat), nrow = length(lvs))
    for (i in 1:length(lvs)) {
      inx <- hc.cls == lvs[i]
      my.mns[i, ] <- apply(hc.dat[inx, ], 2, mean)
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
    myW <- nrow(hc.dat) * 18 + 150
    if (myW < minW) {
      myW <- minW
    }
    w <- round(myW/72, 2)
  } else if (width == 0) {
    w <- 7.2
  } else {
    w <- 7.2
  }
  
  myH <- ncol(hc.dat) * 18 + 150
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
    w <- nrow(hc.dat) * 25 + 300
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
  rownames(annotation) <- rownames(hc.dat)
  
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
  pheatmap::pheatmap(t(hc.dat), annotation = annotation, 
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
### my_PlotSubHeatMap - get PlotSubHeatMap function work ####
my_PlotSubHeatMap = function (mSetObj = NA, imgName, format = "png", dpi = 72, width = NA, 
                              dataOpt, scaleOpt, smplDist, clstDist, palette, method.nm, 
                              top.num, viewOpt, rowV = T, colV = T, border = T, grp.ave = F) 
{
  var.nms = colnames(mSetObj$dataSet$norm)
  if (top.num < length(var.nms)) {
    if (method.nm == "tanova") {
      if (mSetObj$dataSet$cls.num == 2) {
        if (is.null(mSetObj$analSet$tt)) {
          Ttests.Anal(mSetObj)
          #mSetObj <- .get.mSet(mSetObj)
        }
        #var.nms <- names(sort(mSetObj$analSet$tt$p.value))[1:top.num]
        var.nms <- names(sort(mSetObj$analSet$aov$p.value))[1:top.num]
      }
      else {
        if (is.null(mSetObj$analSet$aov)) {
          ANOVA.Anal(mSetObj)
          #mSetObj <- .get.mSet(mSetObj)
        }
        var.nms <- names(sort(mSetObj$analSet$aov$p.value))[1:top.num]
      }
    }
    else if (method.nm == "cor") {
      if (is.null(mSetObj$analSet$cor.res)) {
        Match.Pattern(mSetObj)
        #mSetObj <- .get.mSet(mSetObj)
      }
      cor.res <- mSetObj$analSet$cor.res
      ord.inx <- order(cor.res[, 3])
      cor.res <- cor.res[ord.inx, ]
      ord.inx <- order(cor.res[, 1])
      cor.res <- cor.res[ord.inx, ]
      var.nms <- rownames(cor.res)[1:top.num]
    }
    else if (method.nm == "vip") {
      if (is.null(mSetObj$analSet$plsda)) {
        PLSR.Anal(mSetObj)
        PLSDA.CV(mSetObj)
        #mSetObj <- .get.mSet(mSetObj)
      }
      vip.vars <- mSetObj$analSet$plsda$vip.mat[, 1]
      var.nms <- names(rev(sort(vip.vars)))[1:top.num]
    }
    else if (method.nm == "rf") {
      if (is.null(analSet$rf)) {
        RF.Anal(mSetObj)
        #mSetObj <- .get.mSet(mSetObj)
      }
      var.nms <- GetRFSigRowNames()[1:top.num]
    }
  }
  var.inx <- match(var.nms, colnames(mSetObj$dataSet$norm))
  PlotHeatMap(mSetObj, imgName, format, dpi, width, dataOpt, 
              scaleOpt, smplDist, clstDist, palette, viewOpt, rowV, 
              colV, var.inx, border, grp.ave)
}






# main ####

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# setwd("./2liver_quad_QE+pos")
project_dir = getwd()

giant_df = read.csv("result.csv")

# data cleaning
giant_df = data_clean(giant_df)
giant_colname = colnames(giant_df)
data_raw = giant_df[,c(1,2,grep("[[:digit:]]",giant_colname))]

# Intensity dataframe 
data_all = intensity_dataframe(giant_df)

# Statistics dataframe
data_statistics = giant_df[,1:which(colnames(giant_df)=="mean_log_p")]

# Plot scatter graph for top hits in each tissue
# plot_scatter_by_tissue(data_all, data_statistics, top_n = 9, nrow=3, ncol=3)

# Plot scatter graph for top hits found in > 5 tissue 
# plot_scatter_by_tissue_num(data_all, data_statistics, signif_counts = 2, signif_level = 0.05)

# Select peaks
{
  
  signif_level = 0.1
  signif_counts = 5
  signif_level = -log10(signif_level)
  data_statistics["signif_counts"] = apply(data_statistics, 1, function(x) {
    x = as.numeric(x[13:22])
    non_na = !is.na(x)
    signif = sum(x[non_na]>signif_level, na.rm = T)
    return(signif)
  })
  data_signif_counts = data_statistics[with(data_statistics, order(-signif_counts, -sum_log_p)),]
  data_signif_counts = data_signif_counts[data_signif_counts$signif_counts>=signif_counts,]
  
  top_n = 5
  tissue = levels(data_all$tissue)
  data_list = list()
  for(j in 1:length(tissue)){
    data_select = data_statistics[with(data_statistics, order(-eval(parse(text = paste(tissue[j],"_y",sep=""))))),]
    data_select = data_select[1:top_n,]
    data_list[[length(data_list)+1]] = data_select[!is.na(data_select[paste(tissue[j],"_y",sep="")]),]
  }
  data_select = bind_rows(data_list)
  data_select = rbind(data_select, data_signif_counts)
  data_select = data_select[!duplicated(data_select),]
  
  figure_list = list()
  for(i in 1:nrow(data_select)){
    medMz = data_select$medMz[i]
    medRt = data_select$medRt[i]
    plot_data = data_all[data_all$medMz == medMz & data_all$medRt == medRt,]
    figure_list[[length(figure_list)+1]] = plot_bar_scatter_graph(plot_data = plot_data, bar_plot = T)
  }
  
  ml_ls = list()
  for(i in 1:length(figure_list)){
    ml_ls[[i]] = marrangeGrob(figure_list[[i]], nrow = 2, ncol = 2)
  }
  
  pdf("mdata_select.pdf",width = 6.5, height = 5)
  print(ml_ls)
  dev.off()
}

g_vertex_ILP = read.csv("g_vertex_ILP.txt", stringsAsFactors = F)

# Mdata = normalization_each_tissue(data_all)
Mdata = data_raw
Mdata_mz = Mdata$medMz
Mdata_RT = Mdata$medRt
row_id = numeric(length = nrow(data_select))
formula_id = rep("NA", nrow(data_select))
for(i in 1:nrow(data_select)){
  medRT = data_select$medRt[i]
  medMz = data_select$medMz[i]
  temp_mz = which(Mdata_mz == medMz)
  temp_RT = which(Mdata_RT == medRT)
  row_id[i] = intersect(temp_mz, temp_RT)
  temp_formula = g_vertex_ILP$formula[abs(g_vertex_ILP$mz - medMz) <0.001 &
                                        abs(g_vertex_ILP$RT - medRT) < 0.1]
  if(length(temp_formula)==1){formula_id[i] = temp_formula}
}

{
  data_signif = Mdata[row_id,]
  # data_signif = Mdata
  sample_names = colnames(data_signif)[3:length(colnames(data_signif))]
  sample_cohort = stri_replace_last_regex(sample_names,'_\\d+|-\\d+', '',stri_opts_regex(case_insensitive=T))
  data_signif["formula"] = formula_id
  data_signif["Merge_ID"] = with(data_signif, paste(formula, medMz, medRt, sep="_" ))
  data_signif["Merge_ID"] = with(data_signif, paste(medMz, medRt, sep="_" ))
  MA_output = data_signif[,c("Merge_ID",sample_names)]
  MA_output = rbind(c("cohort",sample_cohort),MA_output)
  
  write.csv(MA_output, file="MetaboAnalyst_file.csv", row.names=F)
}

# Legacy code ####
# # Plot scatter graph with given mz and RT
# {
#   medMz = 131.09462
#   medRt = 5.606
#   plot_data = data_all[abs(data_all$medMz - medMz) <0.001 &
#                          abs(data_all$medRt - medRt) < 0.1,]
#   figure_select = plot_bar_scatter_graph(plot_data = plot_data, bar_plot = T)
#   figure_select
# }
# profvis::profvis(for(i in 1:100){}) ####

## My_plot_heatmap ####
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "MetaboAnalyst_file.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
# mSet<-Normalization(mSet, "None", "LogNorm", "MeanCenter", ratio=FALSE, ratioNum=20)
raw_data = mSet$dataSet$norm
cohort = mSet$dataSet$cls


my_plot_heatmap(raw_data, cohort)
  
  
## Volcano plot
  

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




