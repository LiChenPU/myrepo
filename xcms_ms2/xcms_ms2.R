library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(lc8)
library(mzR)

library(microbenchmark)
# library(readr)


# functions ####
# double log transform that enables plot both positive and negative value on a log scale 
# min abs value of flux to be included in the plot
lb=1000
dblog_trans <- function(){
  scales::trans_new(name='dblog', transform = function(x) (log10(abs(x)+lb)-log10(lb))*sign(x),
                    inverse = function(x) sign(x)*(10^(abs(x)+log10(lb))-lb))
}
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

## predict formula ####
my_pred_formula=function(mz = df$mz, inten = df$inten, ion_mode = ion_mode,
                         parent_formula = "C99H100N15O15S3P3", N_rule = F, 
                         ppm=15, db_max=8){
  
  C_range = 0:(elem_num_query(parent_formula, "C"))
  H_range = 0:(elem_num_query(parent_formula, "H"))
  N_range = 0:(elem_num_query(parent_formula, "N"))
  O_range = 0:(elem_num_query(parent_formula, "O"))
  P_range = 0:(elem_num_query(parent_formula, "P"))
  S_range = 0:(elem_num_query(parent_formula, "S"))
  
  predict_formula = character(length(mz))
  #i=16
  for(i in 1:length(mz)){
    temp = mz_formula(mz[i], charge = ion_mode,  N_rule = N_rule, 
                      C_range=C_range, H_range = H_range, N_range=N_range, O_range=O_range, P_range = P_range, S_range = S_range, ppm=ppm, db_max = db_max)
    if(!is.data.frame(temp)){
      predict_formula[i]=round(mz[i],digits = 3)
      next
    }
    temp["abs"]=abs(temp$differ)
    temp = temp[with(temp, order(abs)),]
    predict_formula[i]=temp$formula[1]
  }
  return(predict_formula)
}

## filter_MS2_Spec ####
filter_MS2_Spec = function(MS2ScanData = MS2ScanData,
                           peak_list = peak_list,
                           spec_all = spec_all,
                           targetMzError = 10E-6,
                           targetRtError = 0.3
                           )
{
  targetMzError = 10E-6
  targetRtError = 0.3
  targetMS2Spectra_ls = list()
  for(i in 1:nrow(peak_list)){
    
    targetMz = peak_list$medMz[i]
    targetRt = peak_list$medRt[i]

    targetMS2Scans = MS2ScanData %>%
      filter(abs(precursorMZ - targetMz) < targetMzError*targetMz) %>%
      filter(abs(retentionTime - targetRt*60) < targetRtError*60)
    
    if(nrow(targetMS2Scans)==0){next}
    
    targetMS2Scans = targetMS2Scans %>%
      arrange(-precursorIntensity, -totIonCurrent) %>%
      distinct(fileIdx, .keep_all=T)
    targetMS2Spectra = spec_all[targetMS2Scans$spectrum]
    targetMS2Spectra_ls[[length(targetMS2Spectra_ls) + 1]] = targetMS2Spectra
    names(targetMS2Spectra_ls)[length(targetMS2Spectra_ls)] = paste(i,targetMz, targetRt, sep="_")  
  }
  
  return(targetMS2Spectra_ls)
}

## plot_MS2_spec ####
plot_MS2_spec = function(targetMS2Spectra,
                         show_mz_formula = "formula",
                         fileNames = fileNames, 
                         top_n_peaks = 10
                         )
{
  fig_ls = list()
  for(i in 1:length(targetMS2Spectra)){
    
    temp_spec = targetMS2Spectra[[i]]
    temp_mzs = round(mz(temp_spec),5)
    temp_intens = intensity(temp_spec)
    
    df = as.data.frame(cbind(mz=temp_mzs, inten=temp_intens))
    df = df %>%
      filter(inten > max(inten)*.05 & inten >5E3) %>%
      arrange(-inten) %>%
      top_n(top_n_peaks, inten)
    
    if(show_mz_formula == "formula"){
      ion_mode = polarity(temp_spec)
      df["pred_formula"] = my_pred_formula(df$mz, df$inten, ion_mode = ion_mode)
      ms2Plot = ggplot(df, aes(x=mz, y=inten, ymax = inten, ymin = 0)) +
        geom_linerange() + 
        ggtitle(paste(fileNames[temp_spec@fromFile], "RT =", round(temp_spec@rt/60,3))) + 
        # geom_text_repel(aes(label=mz),colour='red') +
        geom_text_repel(aes(label=pred_formula),colour='blue') +
        xlim(50, ceiling(max(df$mz)/50)*50)
    } else {
      ms2Plot = ggplot(df, aes(x=mz, y=inten, ymax = inten, ymin = 0)) +
        geom_linerange() + 
        ggtitle(paste(fileNames[temp_spec@fromFile], "RT =", round(temp_spec@rt/60,3))) + 
        geom_text_repel(aes(label=mz),colour='red') +
        # geom_text_repel(aes(label=pred_formula),colour='blue') +
        xlim(50, ceiling(max(df$mz)/50)*50)
    }
    
    fig_ls[[length(fig_ls)+1]] = ms2Plot
  }
  return(fig_ls)
}

## print_MS2_spec ####
print_MS2_spec = function(plotsMS2Spectra_ls, 
                          nrow = 4,
                          ncol = 2, 
                          outputFileName = "all_ms2")
{
  ml_list = list()
  for(i in 1:length(plotsMS2Spectra_ls)){
    ml <- marrangeGrob(plotsMS2Spectra_ls[[i]], nrow=nrow, ncol=ncol, 
                       top = names(plotsMS2Spectra_ls[i]))
    ml_list[[length(ml_list)+1]] = ml
  }
  
  if(file_test("-f", paste(outputFileName, ".pdf", sep=""))){
    printtime = Sys.time()
    timestamp <- paste(unlist(regmatches(printtime, gregexpr("[[:digit:]]+", printtime))),collapse = '')
    outputFileName = paste(outputFileName,"_", timestamp, ".pdf", sep="")
  } else {
    outputFileName = paste(outputFileName, ".pdf", sep="")
  }
  

  pdf(outputFileName,
      onefile = TRUE,
      width = 6.5,
      height = 10.5)
  
  print(ml_list)
  dev.off()
}

## calculate MS2 spec similarity ####
spec2mzIntensity = function(spec, top_n_peaks=10){
  mz1 = mz(spec)
  inten1 = intensity(spec)
  inten1 = inten1/max(inten1)
  spec_df = as.data.frame(cbind(mz=mz1, inten=inten1))
  spec_df = spec_df %>%
    arrange(inten) %>%
    top_n(top_n_peaks, inten) %>%
    arrange(mz)
  return(spec_df)
}

mergeMzIntensity = function(spec1_df, spec2_df, ppmTol = 10E-6){
  mz1 = spec1_df$mz
  mz2 = spec2_df$mz
  mzs = c(mz1, mz2)
  mzs = sort(mzs)
  inten1 = spec1_df$inten
  inten2 = spec2_df$inten
  keeps = c(1)
  count = 1
  MS_group = rep(1,(length(mzs)))
  for(i in 2:length(mzs)){
    if(mzs[i]-mzs[i-1]>mzs[i-1]*ppmTol){
      count = count+1
      keeps=c(keeps,i)
    }
    MS_group[i]=count
  }
  
  mz_result = mzs[keeps]
  count_mzs = 2
  count_keep = 1
  i = j = 1
  intens_mat = matrix(0, ncol=3, nrow = length(mz_result))
  while(count_mzs <= length(mzs)){
    if(MS_group[count_mzs] == MS_group[count_mzs-1]){
      intens_mat[count_keep,2] = inten1[i]
      intens_mat[count_keep,3] = inten2[j]
      i = i+1
      j = j+1
      count_mzs = count_mzs+2
    } else if(is.na(mz1[i])){
      intens_mat[count_keep,3] = inten2[j]
      count_mzs = count_mzs+1
      j=j+1
    } else if(mz1[i] == mzs[count_mzs-1]){
      intens_mat[count_keep,2] = inten1[i]
      i=i+1
      count_mzs = count_mzs+1
    } else {
      intens_mat[count_keep,3] = inten2[j]
      j=j+1
      count_mzs = count_mzs+1
    }
    count_keep = count_keep+1
  }
  # Handle boundary
  if(!is.na(mz1[i])){intens_mat[count_keep,2] = inten1[i]}
  if(!is.na(mz2[j])){intens_mat[count_keep,3] = inten2[j]}
  
  intens_mat[,1] = mz_result
  return(intens_mat)
}

# mergeMzIntensity_deprec = function(spec1_df, spec2_df, ppmTol = 10E-6){
#   # if(identical(spec1_df, spec2_df))
#   colnames(spec1_df)[2] = "inten1"
#   colnames(spec2_df)[2] = "inten2"
#   spec_df = merge(spec1_df, spec2_df, all=T, by = "mz")
#   
#   MS_group = rep(1,(nrow(spec_df)))
#   mzs = spec_df$mz
#   keeps = c(1)
#   count = 1
#   for(i in 2:length(mzs)){
#     if(mzs[i]-mzs[i-1]>mzs[i-1]*ppmTol){
#       count = count+1
#       keeps=c(keeps,i)
#     }
#     MS_group[i]=count
#   }
#   spec_df[is.na(spec_df)]=0
#   
#   k_max=k_min=1
#   while (k_max <= length(MS_group)){
#     k_min = k_max
#     while (MS_group[k_min] == MS_group[k_max]){
#       k_max = k_max+1
#       if(k_max > length(MS_group)){break}
#     }
#     # spec_df$mz[k_min]=max(spec_df$mz[k_min:(k_max-1)], na.rm = T)
#     spec_df$inten1[k_min]=max(spec_df$inten1[k_min:(k_max-1)], na.rm = T)
#     spec_df$inten2[k_min]=max(spec_df$inten2[k_min:(k_max-1)], na.rm = T)
#   }
#   spec_df = spec_df[keeps,]
#   return(spec_df)
# }

DotProduct = function(a, b){
  DP = a%*%b / sqrt(a%*%a * b%*%b)
  return(as.numeric(DP))
}

## cluster_mat ####
cluster_mat = function(mat, distance = "euclidean", method = "complete"){
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
# Main ####
## Read files ####
{
  setwd("C:/study/data/exactive/190731 Melanie young old mice MS2/")
  
  list.files()
  peak_list = read.csv("select_peak_list.csv", stringsAsFactors = F)
  
  
  mzXML <- dir(full.names = TRUE, recursive = F, pattern = ".mzXML")
  pd <- data.frame(sample_name = sub(basename(mzXML), pattern = ".mzXML",
                                     replacement = "", fixed = TRUE),
                   stringsAsFactors = FALSE)
  raw_data <- readMSData(files = mzXML, pdata = new("NAnnotatedDataFrame", pd),
                         mode = "onDisk")
}

## Extract spectra and MS2 data #### 
{
  spec_all = xcms::spectra(raw_data)
  scanData = raw_data@featureData@data
  MS2ScanData = scanData %>%
    filter(msLevel==2) %>%
    filter(totIonCurrent>1E5) 
  unique(MS2ScanData$precursorMZ)
 
  # Update precursorIntensity ###W
  # Previous precusorIntensity is calculated by 
  # 1) go to previous full scan
  # 2) look for mz fall within precusorMz and 2 mass unit error
  # 3) take the maximum intensity number
  # Here the precusorIntensity is updated by 
  # 1) go into the exact MS2 scan spectrum
  # 2) look for mz fall within precusorMz and ~ 10ppm error
  # 3) take the intensity number (generally there is only one peak within 10ppm)
  
  targetMS2Spectra = spec_all[MS2ScanData$spectrum]
  targetMzError = 10E-6
  
  mzs_targetMz_position = lapply(targetMS2Spectra, function(x){
    targetMz = x@precursorMz
    mzs = mz(x)
    which(abs(mzs - targetMz) < targetMzError*targetMz & min(abs(mzs - targetMz)))
  })  
  intens = lapply(targetMS2Spectra, intensity)
  intens_targetMz = mapply(function(X,Y){
    if(length(Y)==0){return(0)}
    X[Y]
  }, X=intens, Y=mzs_targetMz_position)
  MS2ScanData$precursorIntensity = intens_targetMz
}


fileNames = raw_data@phenoData@data$sample_name


## Examine all MS2 data in one file ####
{
  plot_upper_limit = 1E8
  fig_TIC_basepeak = ggplot(MS2ScanData, aes(x = totIonCurrent, 
                          y = basePeakIntensity,
                          color = as.factor(precursorMZ)))+
    geom_point() +
    geom_abline(intercept = 0, slope = 1, alpha=0.5) +
    scale_x_continuous(trans = 'dblog',limit=c(-10,plot_upper_limit), 
                       breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9),
                       labels = fancy_scientific) +
    scale_y_continuous(trans = 'dblog',limit=c(-100,plot_upper_limit), breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9), 
                       labels = fancy_scientific)
  
  fig_TIC_precursor = ggplot(MS2ScanData, aes(x = totIonCurrent, 
                          y = precursorIntensity,
                          color = as.factor(precursorMZ)))+
    geom_point() +
    geom_abline(intercept = 0, slope = 1, alpha=0.5) +
    scale_x_continuous(trans = 'dblog',limit=c(-10,plot_upper_limit), 
                       breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9),
                       labels = fancy_scientific) +
    scale_y_continuous(trans = 'dblog',limit=c(-100,plot_upper_limit), breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9), 
                       labels = fancy_scientific)
  
  fig_RT_precursor = ggplot(MS2ScanData, aes(x = retentionTime, 
                          y = precursorMZ, 
                          size = log10(totIonCurrent)))+
    geom_point(alpha=0.2)
  
  fig_overview = list(fig_TIC_basepeak, fig_TIC_precursor, fig_RT_precursor)
}



## Run functions ####
{
  # filter MS2 based on peak_list
  targetMS2Spectra_ls = filter_MS2_Spec(MS2ScanData, 
                                        peak_list,
                                        spec_all,
                                        targetMzError = 10E-6,
                                        targetRtError = 0.3
  )
  # Plot MS2 spectra, select whether to lable mz or formula on spectra
  plotsMS2Spectra_ls = lapply(targetMS2Spectra_ls, plot_MS2_spec, 
                              show_mz_formula = "mz",
                              fileNames = fileNames, 
                              top_n_peaks = 10)
  
  plotsMS2Spectra_ls = lapply(targetMS2Spectra_ls, plot_MS2_spec,
                              show_mz_formula = "formula",
                              fileNames = fileNames,
                              top_n_peaks = 10)
  
  # Print out plots
  print_MS2_spec(plotsMS2Spectra_ls, 
                 nrow = 4,
                 ncol = 2, 
                 outputFileName = "all_ms2_formula")
  
  
  # calculate dot product distance between spectra
  testMs2Spectra = unlist(targetMS2Spectra_ls)
  spec_df_ls = lapply(testMs2Spectra, spec2mzIntensity, top_n_peaks = 10)
  spec_DP_matrix = matrix(1, nrow = length(spec_df_ls), ncol = length(spec_df_ls))
  for(i in 1:length(spec_df_ls)){
    for(j in 1:length(spec_df_ls)){
      spec_merge_df = mergeMzIntensity(spec_df_ls[[i]], spec_df_ls[[j]], ppmTol = 10E-6)
      spec_DP_matrix[i,j] = DotProduct(spec_merge_df[,2], spec_merge_df[,3])
    }
  }
  
  spec_DP_matrix_sym = 0.5 *(spec_DP_matrix + t(spec_DP_matrix))
  pheatmap::pheatmap(spec_DP_matrix_sym, clustering_method = "average")
  
  
  hc.cls = factor(rep(names(targetMS2Spectra_ls), lapply(targetMS2Spectra_ls, length)))
  annotation <- data.frame(class = hc.cls)
  rownames(annotation) <- unlist(lapply(targetMS2Spectra_ls, names))
  colnames(spec_DP_matrix_sym) = rownames(spec_DP_matrix_sym) = unlist(lapply(targetMS2Spectra_ls, names))

  if(length(unique(hc.cls)) < 9){
    pal9 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
             "#FFFF33", "#A65628", "#F781BF", "#999999")
    dist.cols = pal9[1:length(unique(hc.cls))]
  } else {
    pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
              "#FFFF99", "#B15928")
    dist.cols <- colorRampPalette(pal12)(length(unique(hc.cls)))
  }
  # barplot(rep(1,length(dist.cols)), col=dist.cols)
  cols = dist.cols[as.numeric(hc.cls)]
  uniq.cols <- unique(cols)
  cls <- hc.cls
  # names(uniq.cols) <- unique(as.character((cls)))
  names(uniq.cols) <- unique(as.character(sort(cls)))
  ann_colors <- list(class = uniq.cols)
  
  test = cluster_mat(spec_DP_matrix_sym, "euclidean", "complete")
  
  spec_DP_matrix_cls = spec_DP_matrix_sym[test$order, test$order]
  
  pheatmap::pheatmap(spec_DP_matrix_cls, 
                     annotation = annotation, 
                     # fontsize = 8, 
                     # fontsize_row = 8, 
                     # clustering_distance_rows = smplDist,
                     # clustering_distance_cols = smplDist, 
                     # clustering_method = clstDist, 
                     # border_color = border.col, 
                     # cluster_rows = T,
                     # cluster_cols = T,
                     # scale = scaleOpt, 
                     # color = colors,
                     annotation_colors = ann_colors)
  
}








targetMS2Spectra[[row.names(spec_DP_matrix_cls)[1]]]
  
check_identical = matrix(nrow=105, ncol=20)
i=7
system.time(
for(i in 1:length(spec_df_ls)){
  for(j in 1:20){
    spec1_df = spec_df_ls[[j]]
    spec2_df = spec_df_ls[[i]]
    spec_df1 = mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6)
    # spec_df2 = mergeMzIntensity2(spec1_df, spec2_df, ppmTol = 10E-6)
    
    
    check_identical[i,j] = (all(identical(spec_df1[,1], spec_df2[,1]),
                                identical(spec_df1[,2], spec_df2[,2]),
                                identical(spec_df1[,3], spec_df2[,3])))
  }
}
)


profvis::profvis({for(repeating in 1:100){
  spec_df1 = mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6)
  spec_df2 = mergeMzIntensity2(spec1_df, spec2_df, ppmTol = 10E-6)
}})
mbm = microbenchmark(
  "m1" = mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6),
  "m2" = mergeMzIntensity2(spec1_df, spec2_df, ppmTol = 10E-6)
)

autoplot(mbm)







# Legacy code ####
# Peak picking - ms is not accurate
# mfp <- MatchedFilterParam(binSize = .6)
# xod_all <- findChromPeaks(raw_data, param = mfp)
# xod <- filterFile(xod_all, file = 2)
# plotChromPeaks(xod)
# peakData = xod@msFeatureData$chromPeaks

# mzs_all <- mz(raw_data)
# mzs_by_file <- split(mzs_all, f = fromFile(raw_data))
# intens_all = intensity(raw_data)
# intens_by_file <- split(intens_all, f = fromFile(raw_data))

# spec = split(spec_all, f = fromFile(raw_data))
# raw = filterFile(raw_data, file = 2)
# mzs = mzs_by_file[[2]]
# intens = intens_by_file[[2]]
