
library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(dplyr)

setwd("C:/Users/lc8/Desktop/sample mzXML data")
list.files()


mzXML <- dir(full.names = TRUE, recursive = F, pattern = ".mzXML")
pd <- data.frame(sample_name = sub(basename(mzXML), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE),
                 stringsAsFactors = FALSE)

raw_data <- readMSData(files = mzXML, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk")


mfp <- MatchedFilterParam(binSize = 6)
xod <- findChromPeaks(raw_data, param = mfp)
xod_2 <- filterFile(xod, file = 2)

MS2 = chromPeakSpectra(xod_2, msLevel = 2L, expandRt = 0, expandMz = 0,
                 ppm = 0, method = c("all", "closest_rt", "closest_mz", "signal"),
                 skipFilled = FALSE, return.type = c("Spectra", "list"))
test = MS2@listData

precursor_ls = lapply(test, precursorMz)
test_ls = lapply(test, mz)


