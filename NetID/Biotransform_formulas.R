library(janitor)
# install.packages("readxl")
library(readxl)
library(readr)
library(dplyr)
library(enviPat)
library(fitdistrplus)

# Read files ####
{
  biotransform_formulas_ls = list()
  for(work_dir in c("Lin_Yeast_Pos", "Wenyun_Yeast_neg", "Wenyun_yeast_pos")){
    
    
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
    # work_dir = "Lin_Yeast_Pos"
    print(work_dir)
    setwd(work_dir)
    
    
    
    biotransform_formulas_ls[[length(biotransform_formulas_ls)+1]] = read.csv("output_new_formulas.csv", stringsAsFactors = F) %>%
      mutate(data = work_dir)
    
    
  }
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}


{
  all = bind_rows(biotransform_formulas_ls)
  
  new = all %>%
    filter(!In_HMDB)
  
  new_inten = new %>%
    filter(sig > log10(5e4))
}
