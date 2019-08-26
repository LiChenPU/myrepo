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
library(cowplot)
library(gridGraphics)
library(lc8)
library(mzR)
library(ChemmineR)
library(ChemmineOB)
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
my_pred_formula=function(mz = df$mz, inten = df$inten, ion_mode,
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
                      C_range=C_range, H_range = H_range, N_range=N_range, O_range=O_range, P_range = P_range, S_range = S_range, 
                      ppm=ppm, db_max = db_max)
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
plot_MS2_spec = function(MS2Spectra,
                         show_mz_formula = "formula",
                         top_n_peaks = 10,
                         exp_inten_cutoff = 5000,
                         ion_mode)
{
  # MS2Spectra = library_files[[2]][[4403]]
  # MS2Spectra =expMS2Spectra_ls[[1]][[1]]
  if(class(MS2Spectra) == "Spectrum2"){
    temp_spec = MS2Spectra
    temp_mzs = round(mz(temp_spec),4)
    temp_intens = intensity(temp_spec)
    
    temp_caption = try(paste(fns[temp_spec@fromFile], "RT =", round(temp_spec@rt/60,3)),T)
    if(inherits(temp_caption, "try-error")){temp_caption = paste(temp_spec@fromFile, "RT =", round(temp_spec@rt/60,3))}
    # ion_mode = polarity(temp_spec)
    df = as.data.frame(cbind(mz=temp_mzs, inten=temp_intens))
    df = df %>%
      arrange(-inten) %>%
      top_n(top_n_peaks, inten) %>% 
      filter(inten > exp_inten_cutoff)
  }
  if(class(MS2Spectra) == "list"){
    temp_mzs = round(MS2Spectra$spectrum[,1],4)
    temp_intens = MS2Spectra$spectrum[,2]
    if(is.data.frame(MS2Spectra$external_id)){
      temp_id = paste0(MS2Spectra$external_id[1,], collapse = ":")
    } else{
      temp_id = MS2Spectra$external_id
    }
    temp_caption = paste(temp_id, MS2Spectra$formula)
    df = as.data.frame(cbind(mz=temp_mzs, inten=temp_intens))
    df = df %>%
      arrange(-inten) %>%
      top_n(top_n_peaks, inten)
  }

  
  if(show_mz_formula == "formula"){
    df["pred_formula"] = my_pred_formula(df$mz, df$inten, ion_mode = ion_mode)
    ms2Plot = ggplot(df, aes(x=mz, y=inten, ymax = inten, ymin = 0)) +
      geom_linerange() + 
      ggtitle(temp_caption) + 
      # geom_text_repel(aes(label=mz),colour='red') +
      geom_text_repel(aes(label= pred_formula),
                      colour='blue', 
                      box.padding = 0.3, size = 3) +
      xlim(50, ceiling(max(df$mz)/50)*50)
  } else {
    ms2Plot = ggplot(df, aes(x=mz, y=inten, ymax = inten, ymin = 0)) +
      geom_linerange() + 
      ggtitle(temp_caption) + 
      geom_text_repel(aes(label=mz),
                      colour='red',
                      box.padding = 0.3, size = 3) +
      # geom_text_repel(aes(label=pred_formula),colour='blue') +
      xlim(50, ceiling(max(df$mz)/50)*50)
  }
  
  return(ms2Plot)
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

## print_txt_expMS2Spec
print_txt_expMS2Spec = function(targetSpec, top_n_inten = 15){
  spec = as.data.frame(cbind(mz(targetSpec), intensity(targetSpec))) 
  colnames(spec) = c("mz", "intensity")
  spec = spec %>%
    arrange(intensity) %>%
    top_n(top_n_inten, intensity) %>%
    arrange(mz)
  write.table(spec,"spec.txt",row.names = F, col.names = F, sep=" ")
}
## Update precursorIntensity ####
updatePrecursorIntensity = function(MS2ScanData, spec_all, targetMzError = 10E-6){
  # Previous precusorIntensity is calculated by 
  # 1) go to previous full scan
  # 2) look for mz fall within precusorMz and 2 mass unit error
  # 3) take the maximum intensity number
  # Here the precusorIntensity is updated by 
  # 1) go into the exact MS2 scan spectrum
  # 2) look for mz fall within precusorMz and ~ 10ppm error
  # 3) take the intensity number (generally there is only one peak within 10ppm)
  targetMS2Spectra = spec_all[MS2ScanData$spectrum]
  mzs_targetMz_position = lapply(targetMS2Spectra, function(x){
    targetMz = x@precursorMz
    mzs = mz(x)
    which(abs(mzs - targetMz) < targetMzError*targetMz & min(abs(mzs - targetMz)))
  })  
  intens = lapply(targetMS2Spectra, intensity)
  intens_targetMz = mapply(function(X,Y){
    if(length(Y)==0){return(0)}
    if(length(Y)>1){return(max(X[Y]))}
    X[Y]
  }, X=intens, Y=mzs_targetMz_position)
  MS2ScanData$precursorIntensity = intens_targetMz
  return(MS2ScanData)
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
  
  # spec1_df = MoNA_MS2_pos_spectra[[i]]
  # spec2_df = spec_df_test
  # true_merge = merge(as.data.frame(spec1_df),spec2_df, by="mz", all=T)
  mz1 = spec1_df[,1]
  mz2 = spec2_df[,1]
  mzs = c(mz1, mz2)
  mzs = sort(mzs)
  inten1 = spec1_df[,2]
  inten2 = spec2_df[,2]
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
mergeMzIntensity_backup = function(spec1_df, spec2_df, ppmTol = 10E-6){
  # if(identical(spec1_df, spec2_df))
  colnames(spec1_df)[2] = "inten1"
  colnames(spec2_df)[2] = "inten2"
  spec_df = merge(spec1_df, spec2_df, all=T, by = "mz")
  
  MS_group = rep(1,(nrow(spec_df)))
  mzs = spec_df$mz
  keeps = c(1)
  count = 1
  for(i in 2:length(mzs)){
    if(mzs[i]-mzs[i-1]>mzs[i-1]*ppmTol){
      count = count+1
      keeps=c(keeps,i)
    }
    MS_group[i]=count
  }
  spec_df[is.na(spec_df)]=0
  
  k_max=k_min=1
  while (k_max <= length(MS_group)){
    k_min = k_max
    while (MS_group[k_min] == MS_group[k_max]){
      k_max = k_max+1
      if(k_max > length(MS_group)){break}
    }
    # spec_df$mz[k_min]=max(spec_df$mz[k_min:(k_max-1)], na.rm = T)
    spec_df$inten1[k_min]=max(spec_df$inten1[k_min:(k_max-1)], na.rm = T)
    spec_df$inten2[k_min]=max(spec_df$inten2[k_min:(k_max-1)], na.rm = T)
  }
  spec_df = spec_df[keeps,]
  return(spec_df)
}
DotProduct = function(a, b){
  DP = a%*%b / sqrt(a%*%a * b%*%b)
  return(as.numeric(DP))
}

## search_MS2_library ####
search_MS2_library = function(MS2_library_spectra,
                              target_spectrum_df,
                              top_n_spectra = 100)
{
  test_score = rep(0,length(MS2_library_spectra))
  error_message = c()
  for(i in 1:length(MS2_library_spectra)){
    # if(i%%10000 ==0){print(i); print(Sys.time())}
    spec_merge_df = try(mergeMzIntensity(MS2_library_spectra[[i]], target_spectrum_df, ppmTol = 10E-6), silent = T)
    if(inherits(spec_merge_df, "try-error")){
      error_message=c(error_message,i)
      spec_merge_df = mergeMzIntensity_backup(MS2_library_spectra[[i]], target_spectrum_df, ppmTol = 10E-6)
    }
    test_score[i] = DotProduct(spec_merge_df[,2], spec_merge_df[,3])
  }
  
  whichpartrev <- function(x, n=30) {
    which(x >= -sort(-x, partial=n)[n])
  }
  
  id = whichpartrev(test_score, top_n_spectra)
  score = test_score[id]
  
  return(data.frame(id = id, score =score))
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

## my_SMILES2structure ####
my_SMILES2structure =function(SMILES){
  SDF = try(smiles2sdf(SMILES), silent=T)
  if(inherits(SDF, "try-error")){
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x=0.5, y=0.5, paste("Invalid SMILES"))
    return(0)
  }
  tryError = try(ChemmineR::plotStruc(SDF[[1]]), silent=T)
  if(inherits(tryError, "try-error")){
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x=0.5, y=0.5, paste("Invalid SDF"))
    return(0)
  }
}
