
# Main ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("xcms_ms2_functions.R")

## Read files ####
{
  work_dir = "C:/study/data/exactive/Yeast_unknown/190402 MS2 yeast 12 13C 50D2O histidine"
  setwd(work_dir)
  
  list.files()
  peak_list = read.csv("select_peak_list.csv", stringsAsFactors = F)
  
  
  mzXML <- dir(full.names = TRUE, recursive = F, pattern = ".mzXML")
  fns = sub(".mzXML","", basename(mzXML))
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
    filter(totIonCurrent>5E3) 
  unique(MS2ScanData$precursorMZ)
  
  MS2ScanData = updatePrecursorIntensity(MS2ScanData,
                                         spec_all,
                                         targetMzError = 10E-6)
}



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
  
  MS2_overview = list(fig_TIC_basepeak, fig_TIC_precursor, fig_RT_precursor)
}



## Extract MS2 of interests ####
{
  # filter MS2 based on peak_list
  # for each peak in peak_list, store a list of MS2 from all sapmles within mz and rt error
  expMS2Spectra_ls = filter_MS2_Spec(MS2ScanData, 
                                     peak_list,
                                     spec_all,
                                     targetMzError = 10E-6,
                                     targetRtError = 0.3)
  
  
  save(fns, MS2_overview, peak_list, expMS2Spectra_ls, file = paste(basename(getwd()), "EXPMS2.RData",sep="_"))
  
}
