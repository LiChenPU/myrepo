# Import ####
{
  library(xcms)
  library(readr)
  library(pheatmap)
  library(dplyr)
  # library(RColorBrewer)
  # library(pander)
  # library(magrittr)
  # library(ggplot2)
  # library(ggrepel)
  # library(gridExtra)
  # library(cowplot)
  # library(gridGraphics)
  # library(lc8)
  # library(mzR)
  # library(proxy)
  # library(rdist)
  # library(clue)
  
  
  library(dtw)
  library(dtwclust)
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  
}

# Function ####
DotProduct = function(a, b){
  min_len = min(length(a),length(b))
  a = a[1:min_len]
  b = b[1:min_len]
  DP = a%*%b / sqrt(a%*%a * b%*%b)
  return(as.numeric(DP))
}

## Chrom ####
Chrom = function(scanData,
                 mz_i=NA, rt_i=NA, 
                 ion_mode = -1,
                 rt_tol = 1, # min
                 mz_tol = 20, # ppm
                 compound_name="",
                 normalized = T,
                 spline = T){
  
  # profvis::profvis({for(repeating_num in 1:100){
  # mz_i=known_rt_IV$mass[1]
  # rt_i=known_rt_IV$rt[1]
  # ion_mode = -1
  # rt_tol = 1 # min
  mz_tol = mz_tol/1e6 # ppm
  # compound_name="leucine-isoleucine"
  
  H_mass = 1.007825032
  e_mass = 0.000548579
  proton_mass = H_mass - e_mass
  if(is.na(mz_i)){
    
    entry_i = known_rt %>%
      filter(compound == compound_name)
    mz_i = entry_i$mass[1] + proton_mass*ion_mode
    rt_i = entry_i$rt[1]
  } 
  if(is.na(mz_i)){
    stop("check entry name or mz_i")
  }
  
  scanData_filter = scanData %>%
    # filter(retentionTime<rt_i[1])
    filter(abs(retentionTime - rt_i*60) < rt_tol*60)
  
  temp_mzs = mzs[rownames(scanData_filter)]
  temp_inten = intensities[rownames(scanData_filter)]
  
  rn = rownames(scanData_filter)
  
  select = function(mzs, intens){
    scan_filter = mzs >= mz_i*(1-mz_tol) & mzs <= mz_i*(1+mz_tol)
    sum(intens[scan_filter], na.rm = T)
  }
  results = mapply(select, temp_mzs, temp_inten) %>%
    split(scanData_filter$fileIdx)
  
  if(normalized){
    results2 = sapply(results, function(x){
      y = x
      y[is.na(y) | y < 100]=0
      if(any(y>0)){
        y = y/as.vector(sqrt(y%*%y))
      }
      y
    })
  }
  
  if(isTRUE(spline)){
    names_RT_mapping = scanData$retentionTime
    names(names_RT_mapping) = rownames(scanData)
    
    results = lapply(results, function(inten){
      x = names_RT_mapping[names(inten)]
      y = inten
      spline(x,y,n=101)$y
    })
  }
  
  list(results)
  
}
## get_TIC_inten_ls ####
get_TIC_inten_ls = function(scanData){
  # inten_ls = sapply(1:length(fns), function(x){
  #   y = chrs@.Data[[x]]@intensity
  #   y/as.vector(sqrt(y%*%y))
  # })
  
  scanData_ms1 = scanData %>%
    filter(msLevel == 1)
  
  inten_ls = split(scanData_ms1$totIonCurrent, scanData_ms1$fileIdx)
  
  TIC_inten_ls = sapply(inten_ls, function(x){
    y = x/as.vector(sqrt(x%*%x))
  })
  
}

## get_dtw_dist_mtrx ####
get_dtw_dist_mtrx = function(inten_ls){
  
  dataset_dist_mtrx = matrix(0,
                             nrow = length(fns),
                             ncol = length(fns)
  )
  rownames(dataset_dist_mtrx) = colnames(dataset_dist_mtrx) = fns
  for(i in 1:length(fns)){
    for(j in 1:length(fns)){
      if(i<=j){
        next
      }
      dataset_dist_mtrx[i,j] = dtw(inten_ls[[i]],inten_ls[[j]])$distance
      # dataset_dist_mtrx_sbd[i,j] = proxy::dist(inten_ls[[i]], inten_ls[[j]], method = "SBD")
    }
  }
  
  dataset_dist_mtrx = dataset_dist_mtrx+t(dataset_dist_mtrx)
  
}

## get_sbd_dist_mtrx ####
get_sbd_dist_mtrx = function(inten_ls){
  sbd_dist_mtrx = as.matrix(proxy::dist(inten_ls, inten_ls, method = "SBD"))
  rownames(sbd_dist_mtrx) = colnames(sbd_dist_mtrx) = fns
  
  sbd_dist_mtrx
}



## get_sample_group_dist ####
get_sample_group_dist = function(dist_mtrx){
  get_pair_dist = function(dist_mtrx, mat_col_names1, mat_col_names2){
    temp_mat = dist_mtrx[mat_col_names1, mat_col_names2]
    if(identical(mat_col_names1,mat_col_names2)){
      lower_tri_mat = temp_mat[lower.tri(temp_mat)]
      mean_dist = mean(lower_tri_mat)
      sd_dist = sd(lower_tri_mat)
    } else {
      mean_dist = mean(temp_mat)
      sd_dist = sd(temp_mat)
    }
    return(c(mean_dist, sd_dist))
  }
  
  group_dist = matrix(0, nrow = length(name_mapping), ncol = length(name_mapping))
  colnames(group_dist) = rownames(group_dist) = names(name_mapping)
  group_dist_sd = group_dist
  for(i in 1:length(name_mapping)){
    for(j in 1:length(name_mapping)){
      temp_group_dist = get_pair_dist(dist_mtrx, # dist_mtrx_sbd
                                      name_mapping[[i]], name_mapping[[j]])
      group_dist[i,j]=temp_group_dist[1]
      
    }
  }
  group_dist
}
# Main code ####
## File info ####
{
  mzXML <- dir(full.names = TRUE, recursive = T, pattern = ".*.mzXML")
  fns = sub(".mzXML","", basename(mzXML))
  print(mzXML)
  
  name_mapping = gsub("\\.\\/|\\.mzXML","", mzXML) %>%
    strsplit(split = "\\/") %>% data.frame() %>% t() %>% as.data.frame()
  name_mapping = split(name_mapping$V2,name_mapping$V1)
}
## Load files ####
{
  temp_files = list()
  for(i in 1:length(name_mapping)){
    
    foldername = names(name_mapping)[i]
    print(i)
    if(exists(paste(foldername,".RData"))){
      load(paste(foldername,".RData"))
    } else {
      pd <- data.frame(sample_name = name_mapping[[1]],
                       stringsAsFactors = FALSE)
      raw_data <- readMSData(files = mzXML, pdata = new("NAnnotatedDataFrame", pd),
                             mode = "onDisk")
      
      raw_data_ms1 = raw_data %>% filterMsLevel(1)
      
      mzs <- mz(raw_data_ms1)
      intensities <- intensity(raw_data_ms1)
      
      scanData = raw_data_ms1@featureData@data
      
      temp_file = list(scanData = scanData,
                       mzs = mzs,
                       intensities = intensities)
      
      temp_files[[foldername]] = temp_file
      
    }
  }
}

## TIC ####
{
  TIC_inten_ls = get_TIC_inten_ls(scanData)
  
  dist_mtrx_dtw = get_dtw_dist_mtrx(TIC_inten_ls)
  pheatmap(dist_mtrx_dtw)
  
  # group distance calculation ###
  group_dtw = get_sample_group_dist(dist_mtrx_dtw) 
  pheatmap(group_dtw)
  
  # Not very useful
  # dist_mtrx_sbd = get_sbd_dist_mtrx(TIC_inten_ls)
  # pheatmap(dist_mtrx_sbd)
  # group_sbd = get_sample_group_dist(dist_mtrx_sbd) 
  # pheatmap(group_sbd)
}


## Individual peaks ####
{
  known_rt = read_csv("updated_rt_0312.csv")
  
  inten_ls = Chrom(mz_i=NA, rt_i=NA, 
                   ion_mode = -1,
                   rt_tol = 1, # min
                   mz_tol = 20, # ppm
                   compound_name="serine",
                   normalized = T,
                   spline = T)
  
  indvidual_SBD_dist = as.matrix(proxy::dist(inten_ls[[1]], method = "SBD"))
  indvidual_cos_dist = as.matrix(proxy::dist(inten_ls[[1]], method = "cosine"))
  rownames(indvidual_SBD_dist) = colnames(indvidual_SBD_dist) = fns
  rownames(indvidual_cos_dist) = colnames(indvidual_cos_dist) = fns
  
  
  pheatmap(indvidual_SBD_dist)
  pheatmap(indvidual_cos_dist)
  
  
  
  
  profvis::profvis({
    summary_inten_ls = mapply(Chrom, known_rt_IV$mass, known_rt_IV$rt,
                              ion_mode = -1)
    
    summary_SBD_dist = lapply(summary_inten_ls, dist, method = "SBD")
    summary_cos_dist = lapply(summary_inten_ls, dist, method = "cosine")
    temp = Chrom(known_rt_IV$mass[i], known_rt_IV$rt[i])
    indvidual_SBD_dist = as.matrix(proxy::dist(temp[[1]], method = "SBD"))
    indvidual_cos_dist = as.matrix(proxy::dist(temp[[1]], method = "cosine"))
  
  })
  
  
  known_rt_IV = known_rt %>%
    filter(reliability == 'IV') %>%
    filter(!is.na(mass) & !is.na(rt))
  

  names(summary_inten_ls) = known_rt_IV$compound
  
}



# Deprecated ####
## normalize_inten_ls ###

normalize_inten_ls = function(inten_ls){
  sapply(inten_ls, function(x){
    y = x
    y[is.na(y) | y < 100]=0
    z = y/as.vector(sqrt(y%*%y))
    z
  })
}
## get_inten_ls ###
# get_inten_ls = function(raw_data, mz_i=NA, rt_i=NA, 
#                         ion_mode = -1, 
#                         rt_tol = 1, # min
#                         mz_tol = 20, # ppm
#                         compound_name=NA){
#   
#   if(is.na(mz_i) & is.na(compound_name)){
#     return(NULL)
#     warning("both mz and compound_name is empty")
#   }
#   
#   H_mass = 1.007825032
#   e_mass = 0.000548579
#   proton_mass = H_mass - e_mass
#   if(is.na(mz_i)){
#     
#     entry_i = known_rt %>%
#       filter(compound == compound_name)
#     mz_i = entry_i$mass[1] + proton_mass*ion_mode
#     rt_i = entry_i$rt[1]
#   } 
#   
#   sub_rd = raw_data_ms1 %>%
#     filterRt((rt_i + c(-rt_tol,rt_tol))*60) %>%
#     filterMz(mz_i * c(1-mz_tol/1e6,1+mz_tol/1e6))
#   sub_chr = chromatogram(sub_rd)
#   
#   test = sub_rd@featureData@data
#   
#   # MSnbase::plot(sub_chr)
#   
# }
# 
# inten_ls = get_inten_ls(raw_data_ms1,
#                         mz_i=NA, rt_i=NA,
#                         ion_mode = -1,
#                         rt_tol = 1, # min
#                         mz_tol = 20, # ppm
#                         compound_name="alanine")