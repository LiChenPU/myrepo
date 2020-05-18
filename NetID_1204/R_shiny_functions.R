# function ####

## show_peak_list ####
show_peak_list = function(ilp_nodes, mz_interest, ion_form, mz_ppm){
  
  # ilp_nodes = dt$ilp_nodes
  # ion_form = "M"
  # mz_interest = 282.25578
  # mz_ppm = 5
  # inten_lb = 3
  # inten_ub= 10
  # mz_lb = 0 
  # mz_ub = 1500
  # rt_lb = 0
  # rt_ub = 30
  
  ilp_nodes_filter = ilp_nodes %>%
    filter(ilp_result > 1e-6) %>%
    dplyr::select(node_id, medMz, medRt, log10_inten, class, formula, ppm_error) %>%
    dplyr::rename(peak_id = node_id)
    
  if(mz_interest != 0){
    if(ion_form == "M"){mz_adjust = mz_interest}
    if(ion_form == "M+H"){mz_adjust = mz_interest - 1.007276}
    if(ion_form == "M-H"){mz_adjust = mz_interest + 1.007276}
    
    ilp_nodes_filter = ilp_nodes_filter %>%
      filter(abs(medMz - mz_adjust) < medMz * mz_ppm * 1e-6)
  }
  
  return(ilp_nodes_filter)
}

# show_peak_list = function(ilp_nodes, mz_interest, ion_form, mz_ppm,
#                           inten_lb = 3, inten_ub= 10, 
#                           mz_lb = 0, mz_ub = 1500, 
#                           rt_lb = 0, rt_ub = 30){
#   
#   # ilp_nodes = dt$ilp_nodes
#   # ion_form = "M"
#   # mz_interest = 282.25578
#   # mz_ppm = 5
#   # inten_lb = 3
#   # inten_ub= 10
#   # mz_lb = 0 
#   # mz_ub = 1500
#   # rt_lb = 0
#   # rt_ub = 30
#   
#   ilp_nodes_filter = ilp_nodes %>%
#     filter(ilp_result > 1e-6) %>%
#     dplyr::select(node_id, medMz, medRt, log10_inten, class, formula, ppm_error) %>%
#     dplyr::rename(peak_id = node_id)
#   
#   if(mz_interest != 0){
#     if(ion_form == "M"){mz_adjust = mz_interest}
#     if(ion_form == "M+H"){mz_adjust = mz_interest - 1.007276}
#     if(ion_form == "M-H"){mz_adjust = mz_interest + 1.007276}
#     
#     ilp_nodes_filter = ilp_nodes_filter %>%
#       filter(abs(medMz - mz_adjust) < medMz * mz_ppm * 1e-6)
#   }
#   
#   ilp_nodes_filter = ilp_nodes_filter %>%
#     filter(medMz < mz_ub, medMz > mz_lb) %>%
#     filter(medRt < rt_ub, medRt > rt_lb) %>%
#     filter(log10_inten < inten_ub, log10_inten > inten_lb)
#   
#   return(ilp_nodes_filter)
# }

