## https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms.html
library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(dplyr)

setwd("C:/Users/lc8/Desktop/sample mzXML data")
list.files()

  
# mzXML <- dir(full.names = TRUE, recursive = F, pattern = ".mzXML")
# pd <- data.frame(sample_name = sub(basename(mzXML), pattern = ".mzXML",
#                                    replacement = "", fixed = TRUE),
#                  sample_group = c(rep("WT", 3), rep("lowOD", 2)),
#                  stringsAsFactors = FALSE)
# 
# raw_data <- readMSData(files = mzXML, pdata = new("NAnnotatedDataFrame", pd),
#                        mode = "onDisk")

cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
            recursive = TRUE)
## Create a phenodata data.frame
pd <- data.frame(sample_name = sub(basename(cdfs), pattern = ".CDF",
                                   replacement = "", fixed = TRUE),
                 sample_group = c(rep("KO", 6), rep("WT", 6)),
                 stringsAsFactors = FALSE)
raw_data <- readMSData(files = cdfs, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk")

head(rtime(raw_data))
mzs <- mz(raw_data)
head(mzs)
mzs_by_file <- split(mzs, f = fromFile(raw_data))
length(mzs_by_file)

bpis <- chromatogram(raw_data, aggregationFun = "max")
group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
names(group_colors) <- c("WT", "KO")

## Plot all chromatograms.
plot(bpis, col = group_colors[raw_data$sample_group])


tc <- split(tic(raw_data), f = fromFile(raw_data))
boxplot(tc, col = group_colors[raw_data$sample_group],
        ylab = "intensity", main = "Total ion current")


## Bin the BPC
bpis_bin <- bin(bpis, binSize = 2)

## Calculate correlation on the log2 transformed base peak intensities
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- raw_data$sample_name

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = raw_data$sample_group)
rownames(ann) <- raw_data$sample_name

## Perform the cluster analysis
pheatmap(cormat, annotation = ann,
         annotation_color = list(group = group_colors))

## Define the rt and m/z range of the peak area
rtr <- c(2700, 2900)
mzr <- c(334.9, 335.1)
## extract the chromatogram
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
plot(chr_raw, col = group_colors[chr_raw$sample_group])

raw_data %>%
  filterRt(rt = rtr) %>%
  filterMz(mz = mzr) %>%
  plot(type = "XIC")



xchr <- findChromPeaks(chr_raw, param = CentWaveParam(snthresh = 1))

featureSpectra(raw_data)
