library(mzR)

setwd("C:/study/data/exactive/190731 Melanie young old mice MS2/test")
list.files()


mzXML <- dir(full.names = TRUE, recursive = F, pattern = ".mzXML")


testMzXML <- openMSfile(mzXML[1])
hd <- header(testMzXML)

MS2ScanData = hd %>%
  filter(msLevel==2) %>%
  filter(totIonCurrent>1E5) %>%
  filter(precursorMZ!=0)
  # filter(precursorIntensity!=0)
