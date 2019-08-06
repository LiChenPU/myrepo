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

# predict formula
my_pred_formula=function(mz = df$mz, inten = df$inten, parent_formula = "C99H100N15O15S3P3", N_rule = F, ppm=15, db_max=8){
  
  C_range = 0:(elem_num_query(parent_formula, "C"))
  H_range = 0:(elem_num_query(parent_formula, "H"))
  N_range = 0:(elem_num_query(parent_formula, "N"))
  O_range = 0:(elem_num_query(parent_formula, "O"))
  P_range = 0:(elem_num_query(parent_formula, "P"))
  S_range = 0:(elem_num_query(parent_formula, "S"))
  
  predict_formula = character(length(mz))
  #i=16
  for(i in 1:length(mz)){
    temp = mz_formula(mz[i], N_rule = N_rule, 
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



# Main ####
setwd("C:/study/data/exactive/190731 Melanie young old mice MS2/")
list.files()
peak_list = read.csv("select_peak_list.csv", stringsAsFactors = F)


mzXML <- dir(full.names = TRUE, recursive = F, pattern = ".mzXML")
pd <- data.frame(sample_name = sub(basename(mzXML), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE),
                 stringsAsFactors = FALSE)
raw_data <- readMSData(files = mzXML, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk")
mzs_all <- mz(raw_data)
mzs_by_file <- split(mzs_all, f = fromFile(raw_data))
intens_all = intensity(raw_data)
intens_by_file <- split(intens_all, f = fromFile(raw_data))

spec_all = xcms::spectra(raw_data)
spec = split(spec_all, f = fromFile(raw_data))

# raw = filterFile(raw_data, file = 2)
# mzs = mzs_by_file[[2]]
# intens = intens_by_file[[2]]

scanData = raw_data@featureData@data
MS2ScanData = scanData %>%
  filter(msLevel==2) %>%
  filter(totIonCurrent>1E5) 
  # filter(precursorIntensity!=0)
unique(MS2ScanData$precursorMZ)

fileNames = raw_data@phenoData@data$sample_name

# Examine all MS2 data in one file
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
}


ml_list = list()
for(j in 1:nrow(peak_list)){
  targetMz = peak_list$medMz[j]
  targetRt = peak_list$medRt[j]
  targetMzError = 0.001
  targetRtError = 0.3
  
  
  
  targetMS2Scans = MS2ScanData %>%
    filter(abs(precursorMZ - targetMz) < targetMzError) %>%
    filter(abs(retentionTime - targetRt*60) < targetRtError*60) %>%
    # precursorIntensity is not a good indicator, 
    # because when precursor is fully fragmented, it gives 0 precursorIntensity
    arrange(-totIonCurrent) %>%
    distinct(fileIdx, .keep_all=T)
  
  if(nrow(targetMS2Scans)==0){next}
  
  targetMS2Spectra = spec_all[targetMS2Scans$spectrum]
  fig_ls = list()
  for(i in 1:length(targetMS2Spectra)){
    
    temp_spec = targetMS2Spectra[[i]]
    
    temp_mzs = round(mz(temp_spec),5)
    temp_intens = intensity(temp_spec)
    
    df = as.data.frame(cbind(mz=temp_mzs, inten=temp_intens))
    df = df %>%
      filter(inten > max(inten)*.05 | inten >5E3)
    df["pred_formula"] = my_pred_formula(df$mz, df$inten)
    
    
    ms2Plot = ggplot(df, aes(x=mz, y=inten, ymax = inten, ymin = 0)) +
      geom_linerange() + 
      ggtitle(fileNames[temp_spec@fromFile]) + 
      # geom_text_repel(aes(label=mz),colour='red') +
      geom_text_repel(aes(label=pred_formula),colour='blue') +
      xlim(50, ceiling(max(df$mz)/50)*50)
    # scale_y_continuous(trans = 'dblog', breaks = c(0,10^4,10^5,10^6,10^7,10^8,10^9),
    #                    labels = fancy_scientific)
    fig_ls[[length(fig_ls)+1]] = ms2Plot
  }
  
  
  ml <- marrangeGrob(fig_ls, nrow=3, ncol=2, 
                     top = paste(j,targetMz, targetRt, sep="_"))
  ml_list[[length(ml_list)+1]] = ml
}

pdf("all_ms2_formula.pdf",onefile = TRUE)
print(ml_list)
dev.off()



# Legacy code ####
# Peak picking - ms is not accurate
# mfp <- MatchedFilterParam(binSize = .6)
# xod_all <- findChromPeaks(raw_data, param = mfp)
# xod <- filterFile(xod_all, file = 2)
# plotChromPeaks(xod)
# peakData = xod@msFeatureData$chromPeaks