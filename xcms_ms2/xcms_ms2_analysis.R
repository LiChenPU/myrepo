
# Main ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("xcms_ms2_functions.R")
source("~/myrepo/mz_calculator/mz_calculator_function.R")

# MS2 analysis ####
{
  # load MS2 library #
  library_dir = "C:/Users/lc8/Documents/myrepo/xcms_ms2/library/"
  library_files_name = list.files(path = library_dir, 
                                  pattern = "MS2.*rds", 
                                  recursive = T)
  library_files_path = paste0(library_dir, library_files_name)
  library_files_all = lapply(library_files_path, readRDS)
  names(library_files_all) = sub(".rds","", basename(library_files_name))
  
  ion_mode = 1
  library_files = library_files_all[grepl("pos", names(library_files_all))]
  rm(library_files_all)
  # load experiment ms2 files
  setwd("C:/Users/lc8/Desktop/RNA-Xiaoyang")
  
  load(list.files(pattern = "EXPMS2.RData"))
  
  dir.create("MS2plots", showWarnings=F)
  setwd("./MS2plots")
 
}

## Print all experimental MS2 spectra ####
{
  plotsMS2Spectra_ls = list()
  show_mz_formula = "ms"
  top_n_peaks = 15
  
  precalculated_formula = NA
  if(show_mz_formula == "formula"){
    spec_list = list()
    for(exp_i in 1:length(expMS2Spectra_ls)){
      # for(exp_i in 4:4){
      expMS2Spectra = expMS2Spectra_ls[[exp_i]]
      selectLargestTIC = sapply(expMS2Spectra, tic)
      expMS2Spectra_select = spec2mzIntensity(expMS2Spectra[[which.max(selectLargestTIC)]], top_n_peaks = 20)
      expMS2Spectra_select["label"] = names(expMS2Spectra_ls[exp_i])
      spec_list[[exp_i]] = expMS2Spectra_select
    }
    spec_list = lapply(unlist(expMS2Spectra_ls), spec2mzIntensity, top_n_peaks = 20)
    test_rawdata = bind_rows(spec_list) %>% mutate(ID = 1:nrow(.))
    
    
    result_formula = mz_calculator(test_rawdata, 
                                   expand_formula_to_library("C11H22N6O4"),
                                   connect_depth = 9,
                                   ion_mode = ion_mode)
    result_formula = result_formula %>%
      filter(!is.na(formula), ILP_result>0.6) %>%
      mutate(mz = round(mz, 4)) %>%
      dplyr::select(ID, formula)
    precalculated_formula = merge(test_rawdata, result_formula) %>% mutate(mz = round(mz, 4))
    
  }
  
  
  for(i in 1:length(expMS2Spectra_ls)){
    expMS2Spectra = expMS2Spectra_ls[[i]]
    fig_ls = list()
    for(j in 1:length(expMS2Spectra)){
      temp_spec = expMS2Spectra[[j]]
      fig_ls[[j]] = plot_MS2_spec(MS2Spectra = temp_spec, 
                                  top_n_peaks = top_n_peaks,
                                  show_mz_formula = show_mz_formula, 
                                  exp_inten_cutoff = 500,
                                  ion_mode = ion_mode,
                                  precalculated_formula = NA)
      
    }
    plotsMS2Spectra_ls[[i]] = fig_ls
  }
  names(plotsMS2Spectra_ls) = names(expMS2Spectra_ls)
  # Print out plots
  print_MS2_spec(plotsMS2Spectra_ls, 
                 nrow = 2,
                 ncol = 1, 
                 outputFileName = paste("all_MS2",show_mz_formula,"top",top_n_peaks,sep = "_"))
  
}

## Search library for spectra similar to experimental spectra ####
{
  exp_i = 4
  library_i = 2
  

  
  
  for(exp_i in 1:length(expMS2Spectra_ls)){
    # for(exp_i in 4:4){
    expMS2Spectra = expMS2Spectra_ls[[exp_i]]
    selectLargestTIC = sapply(expMS2Spectra, tic)
    expMS2Spectra_select = spec2mzIntensity(expMS2Spectra[[which.max(selectLargestTIC)]], top_n_peaks = 10)
    
    MS2_similar_result_ls = list()
    for(library_i in 1:length(library_files)){
      # for(library_i in 1:1){
      temp_MS2_library = library_files[[library_i]]
      temp_MS2_library_spectra = lapply(temp_MS2_library, "[[", "spectrum")
      
      MS2_similar_result = search_MS2_library(MS2_library_spectra = temp_MS2_library_spectra,
                                              target_spectrum_df = expMS2Spectra_select,
                                              top_n_spectra = 50)
      
      MS2_similar_result["smiles"] = unlist(lapply(temp_MS2_library, "[[", "SMILES"))[MS2_similar_result$id]
      MS2_similar_result = MS2_similar_result%>%
        arrange(desc(score)) %>%
        distinct(smiles, .keep_all=T)
      MS2_similar_result["library"] = library_i
      MS2_similar_result_ls[[library_i]] = MS2_similar_result
    }
    
    MS2_similar_result = bind_rows(MS2_similar_result_ls) %>%
      arrange(desc(score)) %>% 
      distinct(smiles, .keep_all=T) %>% 
      top_n(12, score) %>%
      filter(score!=0)
    
    
    # plot library MS2
    fig_ls = list()
    show_mz_formula = "mz"
    for(i in 1:nrow(MS2_similar_result)){
      # for(i in 1:10){
      temp_spec = library_files[[MS2_similar_result$library[i]]][[MS2_similar_result$id[i]]]
      ms2plot = plot_MS2_spec(temp_spec,
                              show_mz_formula = show_mz_formula, 
                              top_n_peaks = 10,
                              ion_mode = 1)
      p2 = ms2plot
      my_SMILES2structure(temp_spec$SMILES)
      p1 = recordPlot()
      fig_ls[[i]] = plot_grid(p1, p2,rel_widths = c(1,1.5), rel_heights = c(.6, 1))
    }
    
    outputName = names(expMS2Spectra_ls[exp_i])
    ml <- marrangeGrob(fig_ls, nrow=4, ncol=1, 
                       top = outputName)
    
    pdf(file = paste0(outputName,show_mz_formula, ".pdf"), width = 6.5, height = 10.5, onefile = TRUE)
    print(ml)
    dev.off()
  }
}





## calculate dot product distance between spectra ####
testMS2Spectra = unlist(expMS2Spectra_ls)
spec_df_ls = lapply(testMS2Spectra, spec2mzIntensity, top_n_peaks = 20)
spec_DP_matrix = matrix(1, nrow = length(spec_df_ls), ncol = length(spec_df_ls))
for(i in 1:length(spec_df_ls)){
  for(j in 1:length(spec_df_ls)){
    spec_merge_df = try(mergeMzIntensity(spec_df_ls[[i]], spec_df_ls[[j]], ppmTol = 10E-6), silent = T)
    if(inherits(spec_merge_df, "try-error")){
      spec_merge_df = mergeMzIntensity_backup(spec_df_ls[[i]], spec_df_ls[[j]], ppmTol = 10E-6)
    }
    spec_DP_matrix[i,j] = DotProduct(spec_merge_df[,2], spec_merge_df[,3])
  }
}

# Plot heatmap for spectra distance
{
  spec_DP_matrix_sym = 0.5 *(spec_DP_matrix + t(spec_DP_matrix))
  pheatmap::pheatmap(spec_DP_matrix_sym, clustering_method = "average")
  
  hc.cls = factor(rep(names(expMS2Spectra_ls), lapply(expMS2Spectra_ls, length)))
  annotation <- data.frame(class = hc.cls)
  rownames(annotation) <- unlist(lapply(expMS2Spectra_ls, names))
  colnames(spec_DP_matrix_sym) = rownames(spec_DP_matrix_sym) = unlist(lapply(expMS2Spectra_ls, names))
  # colnames(spec_DP_matrix_sym) = rownames(spec_DP_matrix_sym) = names(unlist(expMS2Spectra_ls))
  
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
  
  
  specDist_heatmap = pheatmap::pheatmap(spec_DP_matrix_sym, 
                     annotation = annotation, 
                     cellwidth = 6, 
                     cellheight = 6,
                     fontsize = 6,
                     # fontsize_row = 8, 
                     # clustering_distance_rows = smplDist,
                     # clustering_distance_cols = smplDist, 
                     clustering_method = "average",
                     # border_color = border.col, 
                     # cluster_rows = T,
                     # cluster_cols = T,
                     # scale = scaleOpt, 
                     # color = colors,
                     annotation_colors = ann_colors)
  
  pdf(file ="specDist_heatmap.pdf", width = 10.5, height = 6.5, onefile = TRUE)
  print(specDist_heatmap)
  dev.off()
  
}

# Benchmark
# mbm = microbenchmark(
#   "m1" = mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6),
#   "m2" = mergeMzIntensity_deprec(spec1_df, spec2_df, ppmTol = 10E-6),
#   "m3" ={
#     error_message=1
#     test_try = try(mergeMzIntensity(MoNA_MS2_pos_spectra[[i]], spec_df_test, ppmTol = 10E-6), silent = T)
#     if(inherits(test_try, "try-error")){
#       error_message=error_message+1
#       test_try = mergeMzIntensity_deprec(MoNA_MS2_pos_spectra[[i]], spec_df_test, ppmTol = 10E-6)
#     }
#   }
# )
# 
# autoplot(mbm)




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

# evaluate function efficiency

# targetMS2Spectra[[row.names(spec_DP_matrix_cls)[1]]]
#   
# check_identical = matrix(nrow=105, ncol=20)
# i=7
# system.time(
# for(i in 1:length(spec_df_ls)){
#   for(j in 1:20){
#     spec1_df = spec_df_ls[[j]]
#     spec2_df = spec_df_ls[[i]]
#     spec_df1 = mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6)
#     # spec_df2 = mergeMzIntensity2(spec1_df, spec2_df, ppmTol = 10E-6)
#     
#     
#     check_identical[i,j] = (all(identical(spec_df1[,1], spec_df2[,1]),
#                                 identical(spec_df1[,2], spec_df2[,2]),
#                                 identical(spec_df1[,3], spec_df2[,3])))
#   }
# }
# )
# 
# 
# profvis::profvis({for(repeating in 1:100){
#   spec_df1 = mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6)
#   spec_df2 = mergeMzIntensity2(spec1_df, spec2_df, ppmTol = 10E-6)
# }})
# mbm = microbenchmark(
#   "m1" = mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6),
#   "m2" = mergeMzIntensity2(spec1_df, spec2_df, ppmTol = 10E-6)
# )
# 
# autoplot(mbm)