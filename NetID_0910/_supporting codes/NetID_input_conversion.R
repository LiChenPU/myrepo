# library
{
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
}
# Convert WL data format 
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  template_file = "NetID_input_conversion_template.csv"
  template = read_csv(template_file)
  
  work_dir = "WL_liver_neg"
  setwd("..")
  setwd(work_dir)
  
  filename_wl = "pks_liver_neg_buff_2020-02-21.xlsx"
  WL = readxl::read_xlsx(filename_wl,
                         # sheet = "Sheet1",
                         guess_max = 1e6
  ) %>%
    dplyr::rename(medRt = rt,
                  medMz = mz) %>%
    # dplyr::rename(Formula = Formula...32,
    #               Feature = Feature...33,
    #               Background = Background...25) %>%
    filter(T)
  
  raw_data = template %>%
    merge(WL, all = T) %>%
    mutate(groupId = Index,
           pseudo_sample = 10^sig) %>%
    filter(is.na(Background)) %>%
    dplyr::select(colnames(template))
  
  write_csv(raw_data, "raw_data.csv", na="")
}