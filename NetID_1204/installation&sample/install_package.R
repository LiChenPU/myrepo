## Install R-tools ##
## https://cran.r-project.org/bin/windows/Rtools/
## Recommended stallation in the default path

## Run these lines ##
## It may take 5-10 min total ##
install.packages("devtools")
install.packages("BiocManager")

## close R and run the rest of codes ##
install.packages("DT")
install.packages("rstudioapi")
install.packages("shiny")
install.packages("shinyjs")
install.packages("igraph")
install.packages("reactlog")
install.packages("ShinyTester")
install.packages("shinythemes")
install.packages("visNetwork")
install.packages("dplyr")
install.packages("RColorBrewer")
install.packages("stringr")
install.packages("enviPat")
BiocManager::install("ChemmineR")
BiocManager::install("ChemmineOB")
devtools::install_github("LiChenPU/Formula_manipulation")

## close R again and all package should be in place ##
