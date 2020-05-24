# !diagnostics off

# Import library ####

{
  # devtools::install_github("LiChenPU/Formula_manipulation")
  library(lc8)
  library(enviPat)
  library(dplyr)
  library(tidyr)
  # library(fitdistrplus)
  library(slam)
  # library(cplexAPI)
  library(readr)
  library(stringi)
  library(pracma)
  library(igraph)
  library(cplexAPI)
  library(readxl)
  library(stringr)

  
  # setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  options(scipen=999) # Turn off scientific expression
}

# Fucntions ####
# Function for parsing#### 
## read_raw_data for reading ElMAVEN output ####
read_raw_data = function(filename){
  raw_data = read_csv(filename) 
  if("groupId" %in% colnames(raw_data)){
    raw_data = raw_data %>%
      dplyr::rename(id = groupId)
  }
  return(raw_data)
}


## read_manual library ####
read_manual_library = function(manual_library_file){
  if(!file.exists(manual_library_file)){
    # warning("No manual library file found in data folder.")
    return(NULL)
  }
  data("isotopes")
  manual_library = read.csv(manual_library_file, stringsAsFactors = F)
  check_formula = check_chemform(isotopes, manual_library$formula) 
  if(any(check_formula$warning)){
    stop(paste("Check manual library for formula error:", 
               paste(check_formula$new_formula[check_formula$warning], collapse = ", ")))
  }
  manual_library = manual_library %>%
    mutate(formula = check_formula$new_formula) %>%
    mutate(formula = my_calculate_formula(formula, "C1"),
           formula = my_calculate_formula(formula, "C-1")) %>%
    mutate(mass = formula_mz(formula),
           rdbe = formula_rdbe(formula)) %>%
    mutate(note = as.character(note))
  return(manual_library)
}
## Read_rule_table - for Connect_rules ####
Read_rule_table = function(rule_table_file, extend_rule = F){
  data("isotopes")
  Connect_rules = read.csv(rule_table_file,stringsAsFactors = F)
  if(nrow(Connect_rules) == 0){return(Connect_rules)}
  for(i in 1: nrow(Connect_rules)){
    Connect_rules$formula[i] = check_chemform(isotopes,Connect_rules$formula[i])$new_formula
    Connect_rules$formula[i] = my_calculate_formula(Connect_rules$formula[i], "C1")
    Connect_rules$formula[i] = my_calculate_formula(Connect_rules$formula[i], "C1", -1)
    Connect_rules$mass[i] = formula_mz(Connect_rules$formula[i])
  }
  
  # if(extend_rule){
  #   extend_rules = list()
  #   i=3
  #   for(i in 1: nrow(Connect_rules)){
  #     if(Connect_rules$allow_rep[i] <= 1){next}
  #     for(j in 2:Connect_rules$allow_rep[i]){
  #       temp_rule = Connect_rules[i,]
  #       temp_rule$name = paste(temp_rule$name, "x", j, sep="")
  #       temp_rule$formula = my_calculate_formula(temp_rule$formula, temp_rule$formula, 
  #                                                sign = j-1)
  #       temp_rule$mass = temp_rule$mass * j
  #       temp_rule$rdbe = temp_rule$rdbe * j
  #       
  #       extend_rules[[length(extend_rules)+1]] = temp_rule
  #     }
  #   }
  #   
  #   Connect_rules = rbind(Connect_rules, bind_rows(extend_rules))
  #   Connect_rules = Connect_rules %>%
  #     arrange(mass) %>%
  #     dplyr::select(-allow_rep)
  # }
  return(Connect_rules)
}

## Cohort_Info - Data name and cohorts ####
Cohort_Info = function(Mset, first_sample_col_num = 15)
{
  raw = Mset$Raw_data
  all_names=colnames(raw)[first_sample_col_num:ncol(raw)]
  
  if(length(grep("blank|blk", all_names, ignore.case = T))!=0){
    sample_names=all_names[-grep("blank|blk", all_names, ignore.case = T)]
  } else {
    sample_names=all_names
  }
  blank_names=all_names[grep("blank|blk", all_names, ignore.case = T)]
  sample_cohort=stri_replace_last_regex(sample_names,'[:punct:]?[:alnum:]+', '')
  if(length(Mset$Cohort$sample_cohort) != length(Mset$Cohort$sample_names))
  {print("Warning! cohort number does not match sample number.")}
  
  return(list("sample_names"=sample_names,"blank_names"=blank_names, "sample_cohort"=sample_cohort))
}

## Peak_cleanup - Clean up duplicate peaks from peak picking ####
Peak_cleanup = function(Mset, 
                        mz_tol=5/10^6, 
                        rt_tol=0.1,
                        inten_cutoff=500, 
                        high_blank_cutoff = 2,
                        first_sample_col_num = 15
)
{
  
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ion_mode = Mset$Global_parameter$mode
  raw = Mset$Raw_data %>%
    mutate(medMz = medMz - (H_mass-e_mass)*ion_mode)
  
  
  ##Group MS groups
  {
    s = raw[with(raw, order(medMz, medRt)),]
    
    mzs = s$medMz
    count = 1
    MZ_group = rep(1,(length(mzs)))
    for(i in 2:length(mzs)){
      if(mzs[i]-mzs[i-1]>mzs[i-1]*mz_tol){
        count = count+1
      }
      MZ_group[i]=count
    }
    s["MZ_group"]=MZ_group
    }
  
  ##Group RT similar groups based on MS groups
  {
    s2 = s[with(s, order(MZ_group, medRt)),]
    
    rts = s2$medRt
    
    MZRT_group = rep(1,(length(rts)))
    MZ_group = s2$MZ_group
    
    count = 1
    for(i in 2:length(rts)){
      if(MZ_group[i]!=MZ_group[i-1] | rts[i]-rts[i-1]>rt_tol){
        count = count+1
      }
      MZRT_group[i]=count
    }
    s2["MZRT_group"] = MZRT_group
    
  }
  
  # Take median of the mz and rt for peaks with same MZRTgroup
  {
    s3 = s2 %>% 
      arrange(MZRT_group)
    ncol_raw = ncol(raw)
    MZRT_group = s3$MZRT_group
    medMz = s3$medMz
    medRt = s3$medRt
    
    k_max=k_min=1
    while (k_max <= length(MZRT_group)){
      k_min = k_max
      while (MZRT_group[k_min] == MZRT_group[k_max]){
        k_max = k_max+1
        if(k_max > length(MZRT_group)){break}
      }
      if(k_max-k_min ==1){next}
      medMz[k_min:(k_max-1)]=median(medMz[k_min:(k_max-1)], na.rm = T)
      medRt[k_min:(k_max-1)]=median(medRt[k_min:(k_max-1)], na.rm = T)
      temp = s3[k_min:(k_max-1),first_sample_col_num:ncol_raw]
      temp[1,] = apply(temp, 2, function(x){
        if(any(!is.na(x))){
          return(max(x, na.rm = T))
        } else {
          return(NA)
        }
      })
      s3[k_min:(k_max-1), first_sample_col_num:ncol_raw] = temp[1,]
    }
    
    s3$medMz = medMz
    s3$medRt = medRt
  }
  
  #intermediate files, replace below detection number to random small number
  {
    s4=s3 %>%
      mutate(mean_inten = rowMeans(.[,Mset$Cohort$sample_names], na.rm=T)) %>%
      filter(mean_inten > inten_cutoff)
    # s4[,4:ncol(s4)][s4[,4:ncol(s4)]<inten_cutoff]=sample(1:inten_cutoff, 
    #                                                         size=sum(s4[,4:ncol(s4)]<inten_cutoff), 
    #                                                         replace=T)
  }
  
  # Remove high blank
  {
    if(!identical(high_blank_cutoff, F) & length(Mset$Cohort$blank_names) > 0){
      if(identical(high_blank_cutoff, T)){high_blank_cutoff = 2}
      
      s5 = s4 %>%
        mutate(high_blank = rowMeans(s4[,Mset$Cohort$sample_names]) < 
                 rowMeans(s4[,Mset$Cohort$blank_names]) * high_blank_cutoff) %>%
        filter(!high_blank) %>%
        dplyr::select(-"high_blank")
    } else{
      s5 = s4
    }
  }
  
  duplicated = s4 %>%
    filter(MZRT_group %in% unique(.[["MZRT_group"]][duplicated(.[["MZRT_group"]])]))
  
  s6 = s5 %>%
    distinct(MZRT_group, .keep_all=T) %>%
    arrange(id) %>%
    dplyr::rename(Input_id = id) %>%
    mutate(id = 1:nrow(.)) %>%
    dplyr::select(-c("MZ_group", "MZRT_group")) %>%
    # mutate(mean_inten = rowMeans(.[,Mset$Cohort$sample_names], na.rm=T)) %>%
    mutate(log10_inten = log10(mean_inten))
  
  return(s6)
}

## New_function 1204 ####
## Initiate_nodeset ####
Initiate_nodeset = function(Mset){
  NodeSet = apply(Mset$Data, 1, function(x){
    list(
      mz = as.numeric(x["medMz"]),
      RT = as.numeric(x["medRt"]),
      inten = as.numeric(x["log10_inten"]),
      sample_inten = x[Mset$Cohort$sample_names],
      blank_inten = x[Mset$Cohort$blank_names],
      formula = list(),
      MS2 = NULL
    )
  })
  names(NodeSet) = 1:nrow(Mset$Data)
  return(NodeSet)
}
## Initiate_edgeset ####
Initiate_edgeset = function(Mset, NodeSet, mz_tol_abs = 0, mz_tol_ppm = 10, 
                            rt_tol_bio = Inf, rt_tol_nonbio = 0.2){
  mz_tol_ppm = mz_tol_ppm/1e6
  
  temp_mz_list = NodeSet %>% sapply("[[","mz") %>% sort()
  temp_RT_list = NodeSet %>% sapply("[[","RT") 
  temp_id_list = names(temp_mz_list)
  merge_nrow = length(temp_mz_list)
  
  temp_rules = Mset$Empirical_rules %>% arrange(mass)
  timer=Sys.time()
  {
    edge_ls = list()
    for (k in 1:nrow(Mset$Empirical_rules)){
      temp_fg=temp_rules$mass[k]
      temp_deltaRT = ifelse(temp_rules$category[k] == "Biotransform", rt_tol_bio, rt_tol_nonbio)
      
      # Find the i,j combination gives potential transformation 
      ## Memeory efficient, but may be slower
      ## matrix calculation could be a faster & simpler approach 
      temp_edge_list = list()
      i=j=j_pos=1
      while(i<=merge_nrow){
        
        while(TRUE){
          j=j+1
          if(j>merge_nrow){break}
          temp_ms = temp_mz_list[j]-temp_mz_list[i]
          
          mass_tol = max(temp_mz_list[j]*mz_tol_ppm,mz_tol_abs)
          
          if(temp_ms < (temp_fg - mass_tol)){
            j_pos = j # locate the last j that has smaller ms
            next
          }
          
          # Criteria to entry
          if(abs(temp_ms-temp_fg)<mass_tol){
            delta_RT = temp_RT_list[names(temp_mz_list[j])] - temp_RT_list[names(temp_mz_list[i])]
            if(abs(delta_RT) < temp_deltaRT){
              temp_edge_list[[length(temp_edge_list)+1]]= list(node1=temp_id_list[i],
                                                        node2=temp_id_list[j], 
                                                        mass_dif=(temp_ms-temp_fg)/temp_mz_list[j]*1E6)
            }
          }
          if(temp_ms> (temp_fg + mass_tol)){break}
        }
        i = i + 1
        j = j_pos - 1
      }
 
      edge_ls[[k]]= bind_rows(temp_edge_list) %>%
        mutate(category = temp_rules$category[k],
               linktype = temp_rules$formula[k],
               direction = temp_rules$direction[k],
               rdbe = temp_rules$rdbe[k])
      print(paste(temp_rules$category[k], temp_rules$formula[k], nrow(edge_ls[[k]]),"found."))
    }
  }
  
  print(Sys.time()-timer)
  edge_list = bind_rows(edge_ls) %>%
    filter(!linktype == "") %>% # remove data-data isomer connection
    mutate(node1 = as.numeric(node1),
           node2 = as.numeric(node2)) %>%
    # mutate(edge_id = 1:nrow(.)) %>%
    filter(T)
  
  EdgeSet = apply(edge_list, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = as.vector(x["category"]),
      linktype = as.vector(x["linktype"]),
      direction = as.numeric(x["direction"])
    )
  })
  names(EdgeSet) = 1:length(EdgeSet)
  
  return(EdgeSet)
}
## Initiate_libraryset ####
Initiate_libraryset = function(Mset){
  Metabolites_HMDB = Mset$Library_HMDB %>%
    # mutate(category = "Metabolite") %>%
    dplyr::rename(note = accession) %>%
    mutate(origin = "Library_HMDB")
  
  Metabolites_known = Mset$Library_known %>%
    mutate(category = "Metabolite") %>%
    dplyr::rename(note = HMDB) %>%
    mutate(origin = "Library_known") %>%
    mutate(rt = .[,eval(Mset$Global_parameter$LC_method)])

  Adducts = Mset$Empirical_rules %>%
    filter(category == "Adduct") %>%
    # mutate(category = "Artifact") %>%
    dplyr::select(-direction) %>%
    mutate(origin = "Empirical_rules")
  
  Manual = NULL
  if(!is.null(Mset$Manual_library)){
    Manual = Mset$Manual_library %>%
      mutate(origin = "Manual_library")
  }

  
  # Remove entries in HMDB that are adducts 
  Metabolites_HMDB = Metabolites_HMDB %>%
    filter(!formula %in% Adducts$formula)
  
  LibrarySet = bind_rows(Metabolites_HMDB, Metabolites_known, Manual, Adducts) %>%
    distinct(SMILES, formula, .keep_all = T) %>%
    mutate(library_id = (1+nrow(Mset$Data)) : (nrow(Mset$Data)+nrow(.))) %>%
    # group_by(SMILES) %>%
    # filter(n()>1) %>%
    filter(T)

  return(LibrarySet)
}

## Peak_grouping ####
Peak_grouping = function(NodeSet, RT_cutoff = 0.2, inten_cutoff = 1e4)
{
  node_mass = sapply(NodeSet, "[[", "mz")
  node_RT = sapply(NodeSet, "[[", "RT") 
  node_inten = sapply(NodeSet, "[[", "inten")
  temp_id = as.numeric(names(node_RT)) # numeric

  peak_group_ls = list()
  
  for(i in 1:length(node_RT)){
    if(node_inten[i] < log10(inten_cutoff)){next}
    RT_min = node_RT[i] - RT_cutoff
    RT_max = node_RT[i] + RT_cutoff
    partner_id = temp_id[which(node_RT <= RT_max & node_RT >= RT_min)]
    peak_group_ls[[length(peak_group_ls)+1]] = list(node1 = rep(temp_id[i], length(partner_id)),
                                                    node2 = partner_id)
    
  }
  
  peak_group = bind_rows(peak_group_ls) %>%
    mutate(mass1 = node_mass[node1],
           mass2 = node_mass[node2]) %>%
    mutate(RT1 = node_RT[node1],
           RT2 = node_RT[node2]) %>%
    mutate(inten1 = node_inten[node1],
           inten2 = node_inten[node2]) %>%
    filter(T)
  
  return(peak_group)
}

## Ring_artifact_connection ####
Ring_artifact_connection = function(peak_group,
                                        ppm_range_lb = 50, ppm_range_ub = 1000, ring_fold = 50, inten_threshold = 1e6){
  
  ring_artifact = peak_group %>%
    filter(inten1 > log10(inten_threshold)) %>%
    mutate(mz_dif = mass2 - mass1) %>%
    mutate(ppm_mz_dif = mz_dif / mass2 * 1e6) %>%
    filter(abs(ppm_mz_dif) < ppm_range_ub & abs(ppm_mz_dif) > ppm_range_lb) %>%
    mutate(inten_ratio = inten2 - inten1) %>%
    filter(inten_ratio < log10(1/ring_fold))
  
  
  EdgeSet_ring_artifact = apply(ring_artifact, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Ring_artifact",
      linktype = "Ring_artifact",
      direction = 1
    )
  })
  
  return(EdgeSet_ring_artifact)
}

## Oligomer_multicharge_connection ####
Oligomer_multicharge_connection = function(peak_group, ppm_tol = 10){
  oligomer = peak_group %>% 
    mutate(mz_ratio12 = mass1/mass2,
           mz_ratio21 = mass2/mass1) %>% 
    filter(mz_ratio12 > 1.5 | mz_ratio21 > 1.5) %>%
    mutate(mz_ratio12_dif = mz_ratio12 - round(mz_ratio12),
           mz_ratio21_dif = mz_ratio21 - round(mz_ratio21)) %>%
    filter(abs(mz_ratio12_dif) < (ppm_tol / 1e6) | abs(mz_ratio21_dif) < (ppm_tol / 1e6)) %>%
    mutate(direction = ifelse(mz_ratio12 > mz_ratio21, -1, 1),
           ratio = round(pmax(mz_ratio12, mz_ratio21))) 
  
  oligomer1 = oligomer %>%
    filter(direction == 1)
  oligomer2 = oligomer %>%
    filter(direction == -1) %>%
    dplyr::rename(node1 = node2, node2 = node1)
  oligomer = bind_rows(oligomer1, oligomer2) %>%
    distinct(node1, node2, .keep_all = T)
  
  EdgeSet_oligomer = apply(oligomer, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Oligomer",
      linktype = as.character(x["ratio"]),
      direction = 1
    )
  })
  
  EdgeSet_multicharge = apply(oligomer, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Multicharge",
      linktype = as.character(x["ratio"]),
      direction = -1
    )
  })
  
  EdgeSet_oligomer_multicharge = bind_rows(EdgeSet_oligomer, EdgeSet_multicharge)
  
  return(EdgeSet_oligomer_multicharge)
}

## Multicharge_isotope_connection ####
Multicharge_isotope_connection = function(EdgeSet, EdgeSet_oligomer_multicharge){
  
  EdgeSet_natural_abundance = bind_rows(EdgeSet) %>%
    filter(category == "Natural_abundance")
  
  EdgeSet_multicharge = EdgeSet_oligomer_multicharge %>%
    filter(category == "Multicharge") %>%
    mutate(linktype = as.numeric(linktype)) %>%
    dplyr::select(node1, node2, linktype)
  
  edge_multicharge_isotope = EdgeSet_natural_abundance %>%
    filter(node1 %in% EdgeSet_multicharge$node2, node2 %in% EdgeSet_multicharge$node2) %>% 
    merge(EdgeSet_multicharge, by.x = "node1", by.y="node2", suffixes = c("",".parent")) %>%
    merge(EdgeSet_multicharge, by.x = "node2", by.y="node2", suffixes = c("",".isotope")) %>%
    filter(linktype.parent == linktype.isotope) %>%
    mutate(linktype = mapply(my_calculate_formula, linktype, linktype, 
                             sign = -(linktype.parent-1)/linktype.parent)) %>%
    mutate(node1 = node1.parent,
           node2 = node1.isotope,
           category = "Multicharge_isotope") %>%
    dplyr::select(colnames(EdgeSet_natural_abundance))
  
  return(edge_multicharge_isotope)
}
## Heterodimer_connection ####
Heterodimer_connection = function(peak_group, NodeSet, ppm_tol = 10, inten_threshold = 1e5){
  
  peak_group_ls = peak_group %>%
    split(.$node2)
  
  node_inten = sapply(NodeSet, "[[", "inten")
  
  hetero_dimer_ls = list()
  for(i in 1: length(peak_group_ls)){
    temp_e = peak_group_ls[[i]]
    temp_matrix = outer(temp_e$mass1, temp_e$mass1, FUN = "+")  # mz_node1 and mz_dif are the same.
    temp_matrix = (temp_matrix - temp_e$mass2[1])/temp_e$mass2[1] * 10^6
    temp_index = which(abs(temp_matrix) < ppm_tol, arr.ind = T)
    if(dim(temp_index)[1]>0){
      temp_ppm = temp_matrix[temp_index]
      temp_node_1 = temp_e$node1[temp_index[,1]]
      linktype = temp_e$node1[temp_index[,2]]
      temp_df = data.frame(node1 = temp_node_1, linktype = linktype, node2 = temp_e$node2[1], mass_dif = temp_ppm)
      hetero_dimer_ls[[length(hetero_dimer_ls)+1]] = temp_df
    }
  }
  
  
  
  hetero_dimer_df1 = bind_rows(hetero_dimer_ls) %>% 
    mutate(inten1 = node_inten[node1], 
           inten2 = node_inten[node2]) %>%
    filter(inten1 > log10(inten_threshold)) %>% # retain only high intensity as node1
    filter(inten1 > inten2) %>% # Heterodimer parent has higher intensity
    filter(node1 != linktype) %>% # remove homo-dimer
    arrange(node1)
  
  hetero_dimer_df2 = hetero_dimer_df1 %>%
    mutate(temp = linktype,
           linktype = node1,
           node1 = temp) %>%
    dplyr::select(-temp)
  
  hetero_dimer_df = bind_rows(hetero_dimer_df1, hetero_dimer_df2) %>%
    distinct(node1, linktype, .keep_all =T)
  
  EdgeSet_heterodimer = apply(hetero_dimer_df, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Heterodimer",
      linktype = as.character(x["linktype"]),
      direction = 1
    )
  })
  
  return(EdgeSet_heterodimer)
}
 
## Experiment_MS2_fragment_connection ####
Experiment_MS2_fragment_connection = function(peak_group, NodeSet, ppm_tol = 10, inten_threshold = 1e5){
  
  # It requires the parent node has msr MS2, and build connection between nodes
  # See comparison to Library_MS2_fragment_connection
  
  node_MS2_logic = sapply(NodeSet, function(x){!is.null(x$MS2)})
  node_MS2 = (1:length(NodeSet))[node_MS2_logic]
  
  if(length(node_MS2)==0){
    return(NULL)
  }
  
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ion_mode = Mset$Global_parameter$mode
  
  peak_group_MS2 = peak_group %>%
    filter(node2 %in% node_MS2 & inten2 > log10(inten_threshold),
           node1 != node2,
           mass1 < mass2) %>%
    split(.$node2)
  
  experiment_MS2_fragment_ls = list()
  for(i in 1:length(peak_group_MS2)){
    
    node2 = as.numeric(names(peak_group_MS2)[i])

    mz_parent = NodeSet[[node2]]$mz + (H_mass-e_mass)*ion_mode
    
    MS2 = NodeSet[[node2]]$MS2
    MS2_mz = unlist(MS2[,1])
    
    temp_e = peak_group_MS2[[i]]
    node1_mz = temp_e$mass1 + (H_mass-e_mass)*ion_mode
    
    temp_matrix = outer(node1_mz, MS2_mz, FUN = "-") 
    temp_index = which(abs(temp_matrix) < ppm_tol * mz_parent / 1e6, arr.ind = T)
    
    # print(i)
    # print(temp_index)
    
    if(dim(temp_index)[1]>0){
      temp_dif = temp_matrix[temp_index]
      temp_node_1 = temp_e$node1[temp_index[,1]]
      temp_node_2 = node2
      
      temp_linktype = mz_parent - MS2_mz[temp_index[,2]]
      
      temp_df = data.frame(node1 = temp_node_1, linktype = "Experiment_MS2_fragment", 
                           node2 = temp_node_2, mass_dif = temp_dif)
      experiment_MS2_fragment_ls[[length(experiment_MS2_fragment_ls)+1]] = temp_df
    }
    
  }
  
  experiment_MS2_fragment_df = bind_rows(experiment_MS2_fragment_ls)
  
  EdgeSet_experiment_MS2_fragment = apply(experiment_MS2_fragment_df, 1, function(x){
    list(
      node1 = as.numeric(x["node1"]),
      node2 = as.numeric(x["node2"]),
      category = "Experiment_MS2_fragment",
      linktype = as.character(x["linktype"]),
      mass_dif = as.numeric(x["mass_dif"]),
      direction = -1
    )
  })
  
  return(EdgeSet_experiment_MS2_fragment)
}
## Library_MS2_fragment_connection ####
Library_MS2_fragment_connection = function(peak_group, FormulaSet, MS2_library,
                                           inten_threshold = 1e5,
                                           ppm_tol = 10, abs_tol = 1e-4){
  # It only requires the parent node has potential database MS2
  
  FormulaSet_df = bind_rows(FormulaSet)
  MS2_library_external_id = sapply(MS2_library, "[[", 'external_id')
  
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  ion_mode = Mset$Global_parameter$mode
  
  FormulaSet_df_MS2 = FormulaSet_df %>%
    # filter(steps == 0) %>%
    filter(steps == 0, transform == "") %>%
    merge(LibrarySet %>%
            dplyr::select(library_id, note, mass), 
          by.x = "parent_id",
          by.y = "library_id") %>%
    mutate(contain_lib_MS2 = note %in% MS2_library_external_id) %>%
    filter(contain_lib_MS2)
  
  peak_group_MS2 = peak_group %>%
    filter(node2 %in% FormulaSet_df_MS2$node_id & inten2 > log10(inten_threshold), #Only look at high-inten peaks
           inten2 > inten1, # Fragment should have lower inten comparaed to parent
           mass2 > mass1) %>%
    split(.$node2)
  
  Library_MS2_fragment = list()
  for(i in 1:length(peak_group_MS2)){
    id_node1 = peak_group_MS2[[i]]$node1
    mz_node1 = sapply(id_node1, function(x){NodeSet[[x]]$mz + (H_mass-e_mass)*ion_mode})
    
    id_node2 = peak_group_MS2[[i]]$node2[1]
    FormulaSet_df_MS2_filter = FormulaSet_df_MS2 %>%
      filter(node_id == id_node2)
    
    node2 = as.numeric(names(peak_group_MS2)[i])
    mz_parent = NodeSet[[node2]]$mz + (H_mass-e_mass)*ion_mode
    
      for(j in 1:nrow(FormulaSet_df_MS2_filter)){
        match_MS2_IDs = which(FormulaSet_df_MS2_filter$note[j] == MS2_library_external_id)
        if(length(match_MS2_IDs) == 0){next}
        match_MS2 = bind_rows(lapply(MS2_library[match_MS2_IDs], "[[","spectrum"))
        
        temp_matrix = outer(mz_node1, match_MS2$mz, FUN = "-") 
        mz_tol = max(mz_node1 * ppm_tol * 1e-6, abs_tol)
        temp_index = which(abs(temp_matrix) < mz_tol, arr.ind = T)
        
        formula_node1 = match_MS2$formula[temp_index[,2]]
        # Some MS2 library contains NULL formulas
        # Clean MS2 library to make sure each fragment has formula
        if(any(is.null(formula_node1))){next}
        
        length_formula_node1 = length(formula_node1)
        
        formula_node2 = FormulaSet_df_MS2_filter$formula[j]
        
        temp_node1 = id_node1[temp_index[,1]]
        
        
        Library_MS2_fragment[[length(Library_MS2_fragment)+1]] = list(node1 = temp_node1,
                                                              node2 = rep(id_node2,length_formula_node1),
                                                              formula1 = formula_node1,
                                                              formula2 = rep(formula_node2,length_formula_node1))
        
      }
    
  }
    
  
  fragment_df = bind_rows(Library_MS2_fragment) %>%
    # group_by(node1, node2) %>%
    # filter(n()>1) %>%
    distinct() %>% 
    # filter(!is.na(formula1)) %>%
    # mutate(rdbe = formula_rdbe(formula1)) %>%
    mutate(linktype = "Library_MS2_fragment",
           category = "Library_MS2_fragment",
           direction = -1)
}

## merge_edgeset ####
merge_edgeset = function(EdgeSet, ...){
  
  EdgeSet_all_df = bind_rows(EdgeSet, ...) %>%
    mutate(edge_id = 1:nrow(.))
  return(EdgeSet_all_df)
}
### Expand_library ####
expand_library = function(lib, rule, direction, category){
  
  if(nrow(rule) == 0){
    return(NULL)
  }
  # initial_lib_adduct_1 = expand_library(lib_adduct, rule_1, direction = 1, category = "Artifact")
  mz_lib = lib$mass
  mz_rule = rule$mass
  mass_exp = outer(mz_lib, mz_rule * direction, FUN = "+")
  
  rdbe_lib = lib$rdbe
  rdbe_rule = rule$rdbe
  rdbe_exp = outer(rdbe_lib, rdbe_rule * direction, FUN = "+")
  
  formula_lib = lib$formula
  formula_rule = rule$formula
  formula_exp = my_calculate_formula(formula_lib, formula_rule, sign = direction)
  
  parent_lib = lib$formula
  transformation_rule = rule$formula
  
  expansion = merge(parent_lib, transformation_rule, all=T) %>%
    dplyr::rename(parent_formula = x,
           transform = y) %>%
    mutate(parent_formula = as.character(parent_formula),
           transform = as.character(transform),
           mass = as.vector(mass_exp),
           rdbe = as.vector(rdbe_exp),
           formula = as.vector(formula_exp),
           direction = direction,
           parent_id = rep(lib$node_id, nrow(rule)),
           category = category) %>%
    dplyr::select(formula, mass, rdbe, category, parent_id, parent_formula, transform, direction)
  return(expansion)
}
### match_library ####
match_library = function(lib, sf, record_ppm_tol, record_RT_tol, current_step, NodeSet){
  lib_mass = lib$mass
  length_lib = length(lib_mass)
  node_mass = sapply(NodeSet, "[[", "mz") %>% sort()
  node_RT = sapply(NodeSet, "[[", "RT")[names(node_mass)]
  temp_id = as.numeric(names(node_mass)) # numeric
  
  i=i_min=i_max=1
  while(i <= length(node_mass)){
    mass_tol = node_mass[i]*record_ppm_tol
    # Move i_min to a position larger than lower threshold
    while(lib_mass[i_min] < node_mass[i] - mass_tol & i_min < length_lib){
      i_min = i_min + 1
    }
    # if i_min's position larger than upper threhold, then move up i
    if(lib_mass[i_min] > node_mass[i] + mass_tol){
      i = i+1
      next
    }
    # Arriving here, means i_min reaches maximum or/and i_min's position is larger than lower threhold
    i_max = i_min
    
    # Move i_max to a position that is below upper threhold
    while(lib_mass[i_max] - node_mass[i] < mass_tol & i_max < length_lib){
      i_max = i_max + 1
    }
    i_max = i_max - 1
    
    # if there is no overlap of between i_min above lower threhold and i_max below upper threhold
    # it means i_min is maximum and i_max is maximum - 1, then break
    if(i_min > i_max){break}
    
    # Otherwise, record
    candidate = lib[i_min:i_max,]
    if(record_RT_tol < 999){
      candidate_filter = abs(node_RT[i] - node_RT[as.character(candidate$parent_id)]) < record_RT_tol
      candidate = candidate[candidate_filter, ]
      
      if(nrow(candidate) == 0){
        i = i + 1
        next
      }
      # We may implement another filter to break a cycle a->b->a, but it is not necessary
      # Actually, keeping both propagation direction makes it easier to trace parents
    }
    
    candidate["node_id"] = temp_id[i]
    candidate["steps"] = current_step
    
    sf[[temp_id[i]]] = bind_rows(sf[[temp_id[i]]], candidate)
    
    i = i+1
  }
  return(sf)
}
### match_library_slow ####
match_library_propagation = function(lib, sf, record_ppm_tol, record_RT_tol, current_step, NodeSet){
  
  node_mass = sapply(NodeSet, "[[", "mz")
  node_RT = sapply(NodeSet, "[[", "RT")
  
  lib = lib %>%
    mutate(parent_mz = node_mass[parent_id],
           parent_RT = node_RT[parent_id])
  
  for(i in 1:length(NodeSet)){
    temp_node_mz = NodeSet[[i]]$mz
    temp_node_RT = NodeSet[[i]]$RT
    
    lib_filter = lib %>%
      filter(abs(parent_RT - temp_node_RT) < record_RT_tol) %>%
      filter(abs(mass - temp_node_mz) < pmin(parent_mz, temp_node_mz) * record_ppm_tol)
    
    if(nrow(lib_filter) == 0){
      next
    }
    
    adding = lib_filter %>%
      dplyr::select(-parent_RT, -parent_mz) %>%
      mutate(node_id = i) %>%
      mutate(steps = current_step)
    
    sf[[i]] = bind_rows(sf[[i]], adding)
  }
  return(sf)
  
}
## Initialize formula pool ####
Initilize_empty_formulaset = function(NodeSet){
  # sf_str = data.frame(
  #   formula = as.character(),
  #   mass = as.numeric(), 
  #   rdbe = as.numeric(), 
  #   category = as.character(), 
  #   parent_id = as.numeric(), 
  #   parent_formula = as.character(),
  #   transform = as.character(),
  #   direction = as.integer(),
  #   stringsAsFactors = F
  # )
  node_mass = sapply(NodeSet, "[[", "mz")
  sf = lapply(1:length(NodeSet), function(x) {
    data.frame(
      node_id = x,
      formula = "Unknown",
      mass = node_mass[x], 
      rdbe = 0, 
      category = "Unknown", 
      parent_id = x, 
      parent_formula = "Unknown",
      transform = "",
      direction = 0,
      steps = 0,
      stringsAsFactors = F
    )
  })
  return(sf)
}

## Match_library_formulaset ####
Match_library_formulaset = function(FormulaSet, Mset, NodeSet, LibrarySet, 
                                    expand = T,
                                    ppm_tol = 5e-6){
  ## initialize
  seed_library = LibrarySet %>% 
    mutate(node_id = library_id,
           formula = formula,
           mass = mass,
           parent_id = library_id,
           parent_formula = formula,
           transform = "",
           direction = 1,
           rdbe = rdbe,
           category = category
    ) %>%
    dplyr::select(node_id, formula, mass, rdbe, parent_id, parent_formula, transform, direction, category) %>%
    filter(T)
  
  lib_met = seed_library %>% filter(category == "Metabolite")
  lib_adduct = seed_library %>%
    filter(category == "Adduct") %>%
    filter(T)
  lib_manual = seed_library %>%
    filter(category == "Manual")
  ## Expansion of starting formula/structures ####
  if(expand){
    initial_rule = Mset$Empirical_rules %>%
      filter(category %in% c("Biotransform", "Adduct"))

    rule_1 = initial_rule %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,1))
    rule_2 = initial_rule %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,-1))

    lib_known_id = LibrarySet %>%
      filter(origin == "Library_known") %>%
      pull(library_id)
    lib_known = seed_library %>%
      filter(node_id %in% lib_known_id) %>%
      filter(category == "Metabolite")
      
    lib_known_1 = expand_library(lib_known, rule_1, direction = 1, category = "Metabolite")
    lib_known_2 = expand_library(lib_known, rule_2, direction = -1, category = "Metabolite")
    lib_met = bind_rows(lib_known_1, lib_known_2, lib_met)

    rule_1 = initial_rule %>% filter(category == "Adduct") %>% filter(direction %in% c(0,1))
    rule_2 = initial_rule %>% filter(category == "Adduct") %>% filter(direction %in% c(0,-1))

    lib_adduct_1 = expand_library(lib_adduct, rule_1, direction = 1, category = "Adduct")
    lib_adduct_2 = expand_library(lib_adduct, rule_2, direction = -1, category = "Adduct")

    lib_adduct = bind_rows(lib_adduct, lib_adduct_1, lib_adduct_2)
  }

  initial_lib1 = bind_rows(lib_adduct, lib_manual) %>%
    filter(!str_detect(formula, "((?<!H)-)|(-(?!1))"))
    # filter((!grepl("-|NA", formula)) | grepl("H-(1|2)", formula)) 
  
  initial_lib2 = bind_rows(lib_met) %>%
    filter(!grepl("-|NA", formula)) 
  
  initial_lib = bind_rows(initial_lib1, initial_lib2) %>%
    arrange(mass)
  
  sf = FormulaSet
  sf = match_library(lib = initial_lib, sf, 
                     record_ppm_tol = ppm_tol, 
                     record_RT_tol = Inf, 
                     current_step = 0, 
                     NodeSet)
  
  return(sf)

}


## Check_sys_error ####
Check_sys_error = function(NodeSet, FormulaSet, LibrarySet,
                           RT_match = T){
  Library_known = LibrarySet %>%
    filter(origin == "Library_known")
  
  node_RT = sapply(NodeSet, "[[", "RT")
  library_RT = LibrarySet$rt
  names(library_RT) = LibrarySet$library_id
  
  node_mass = sapply(NodeSet, "[[", "mz")
  Library_known_msr = bind_rows(FormulaSet) %>% 
    filter(parent_id %in% Library_known$library_id & transform == "") %>%
    mutate(msr_mass = node_mass[node_id], 
           mass_dif = mass - msr_mass,
           ppm_mass_dif = mass_dif/mass * 1e6) %>%
    filter(quantile(ppm_mass_dif, 0.1)<ppm_mass_dif,
           quantile(ppm_mass_dif, 0.9)>ppm_mass_dif) %>%
    mutate(msr_rt = node_RT[as.character(node_id)],
           lib_rt = library_RT[as.character(parent_id)]) %>%
    mutate(rt_match = abs(msr_rt - lib_rt) < 1) %>%
    arrange(-rt_match) %>%
    distinct(node_id, .keep_all = T) %>%
    # distinct(formula, .keep_all = T) %>%
    filter(T)
  
  lsq_result = lm(Library_known_msr$mass_dif~Library_known_msr$msr_mass)
    
  if(RT_match){
    Library_known_msr = Library_known_msr %>%
      filter(rt_match)
    lsq_result = lm(Library_known_msr$mass_dif~Library_known_msr$msr_mass)
  }
  
  
  ppm_adjust = as.numeric(lsq_result$coefficients[2] * 10^6)
  abs_adjust = as.numeric(lsq_result$coefficients[1])
  
  
  #
  # plot(Library_known_msr$msr_mass, Library_known_msr$mass_dif)
  ## Normal test 
  # shapiro.test(Library_known_msr$mass_dif)
  # shapiro.test(Library_known_msr$ppm_mass_dif)
  fitdistData = fitdistrplus::fitdist(Library_known_msr$ppm_mass_dif, "norm")
  # if(fitdistData$estimate["mean"] > 10*fitdistData$sd["mean"]){
    plot(fitdistData)
    print(fitdistData)
    
  # }
  
  return(list(ppm_adjust = ppm_adjust, 
              abs_adjust = abs_adjust, 
              fitdistData = fitdistData))
  # if((ppm_adjust + abs_adjust / 250 * 1e6) > 0.5)
  # return(NULL)
}

### propagate_ring_artifact ####
propagate_ring_artifact = function(new_nodes_df, sf, EdgeSet_ring_artifact, NodeSet, current_step){
  node_mass = sapply(NodeSet, "[[", "mz")
  
  node1 = EdgeSet_ring_artifact$node1
  node2 = EdgeSet_ring_artifact$node2
  node1_node2_mapping = bind_rows(EdgeSet_ring_artifact) %>%
    dplyr::select(-c("category", "linktype", "direction"))
  
  new_nodes_df_ring_artifact = new_nodes_df %>% 
    filter(node_id %in% node1) %>% 
    merge(node1_node2_mapping, by.x="node_id", by.y="node1") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Ring_artifact",
           transform = "Ring_artifact", 
           direction = 1,
           steps = current_step,
           formula = paste0("Ring_artifact_", formula)) %>%
    mutate(node_id = node2, 
           mass = node_mass[node_id]) %>%
    dplyr::select(-node2)
  
  # for(i in unique(new_nodes_df_ring_artifact$node_id)){
  #   sf[[i]] = bind_rows(sf[[i]], new_nodes_df_ring_artifact[new_nodes_df_ring_artifact$node_id == i, ])
  # }
  
  return(new_nodes_df_ring_artifact)
}
### propagate_oligomer ####
propagate_oligomer = function(new_nodes_df, sf, EdgeSet_oligomer, NodeSet, current_step){
  
  node_mass = sapply(NodeSet, "[[", "mz")
  
  node1_node2_mapping = bind_rows(EdgeSet_oligomer) %>%
    dplyr::select(-c("category","direction")) %>%
    mutate(linktype = as.numeric(linktype))

  new_nodes_df1 = new_nodes_df %>% 
    filter(node_id %in% node1_node2_mapping$node1) %>% 
    dplyr::select(-c("category")) %>%
    merge(node1_node2_mapping, by.x="node_id", by.y="node1") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Oligomer",
           direction = 1,
           steps = current_step,
           transform = as.character(linktype)) %>%
    mutate(node_id = node2, 
           mass = mass * linktype,
           formula = as.character(mapply(my_calculate_formula, formula, formula, sign = linktype - 1)),
           rdbe = rdbe * linktype) %>%
    dplyr::select(-c("node2","linktype"))
  
  new_nodes_df_oligomer = new_nodes_df1

  return(new_nodes_df_oligomer)
}
### propagate_multicharge ####
propagate_multicharge = function(new_nodes_df, sf, EdgeSet_oligomer, NodeSet, current_step){
  
  node_mass = sapply(NodeSet, "[[", "mz")
  
  node1_node2_mapping = bind_rows(EdgeSet_oligomer) %>%
    dplyr::select(-c("category","direction")) %>%
    mutate(linktype = as.numeric(linktype))
  
  new_nodes_df2 = new_nodes_df %>% 
    filter(node_id %in% node1_node2_mapping$node2) %>% 
    dplyr::select(-c("category")) %>%
    merge(node1_node2_mapping, by.x="node_id", by.y="node2") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Multicharge",
           direction = -1,
           steps = current_step,
           transform = as.character(linktype)) %>%
    mutate(node_id = node1, 
           mass = mass / linktype,
           formula = as.character(mapply(my_calculate_formula, formula, formula, sign = -(linktype-1)/linktype)),
           rdbe = rdbe / linktype) %>%
    dplyr::select(-c("node1","linktype"))
  
  new_nodes_df_oligomer = new_nodes_df2
  
  return(new_nodes_df_oligomer)
}

### propagate_heterodimer ####
propagate_heterodimer = function(new_nodes_df, sf, EdgeSet_heterodimer, NodeSet, 
                                 current_step, propagation_ppm_threshold, temp_propagation_category){
  node_mass = sapply(NodeSet, "[[", "mz")
  node1_node2_mapping = bind_rows(EdgeSet_heterodimer) %>%
    dplyr::select(-c("category", "direction"))

  new_nodes_df_heterodimer = new_nodes_df %>% 
    filter(node_id %in% node1_node2_mapping$node1) %>% 
    merge(node1_node2_mapping, by.x="node_id", by.y="node1") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Heterodimer",
           direction = 1,
           steps = current_step,
           transform = linktype) %>%
    dplyr::select(-c("linktype")) 
  
  if(nrow(new_nodes_df_heterodimer) == 0){
    return(NULL)
  }
  
  heterodimer_ls = list()
  
  for(i in 1:nrow(new_nodes_df_heterodimer)){
    temp = new_nodes_df_heterodimer[i,]
    transform = sf[[as.numeric(temp$transform)]] %>%
      distinct(formula, .keep_all = T) %>%
      filter(steps <= 0.01 | (steps >= 1  & steps <= 1.01)) %>%
      filter(category %in% temp_propagation_category) %>%
      filter(abs(node_mass[as.character(node_id)] - mass) < propagation_ppm_threshold * mass) # propagate from accurate formulas
      
    if(nrow(transform) == 0){next}

    heterodimer_ls[[length(heterodimer_ls)+1]] = list(parent_id = rep(temp$parent_id, length(transform$formula)),
                                                      transform = rep(temp$transform, length(transform$formula)),
                                                      node2 = rep(temp$node2, length(transform$formula)),
                                                      parent_formula = rep(temp$parent_formula, length(transform$formula)),
                                                      formula = my_calculate_formula(temp$formula, transform$formula),
                                                      rdbe = temp$rdbe + transform$rdbe,
                                                      mass = temp$mass + transform$mass)
    
  }
  
  heterodimer = bind_rows(heterodimer_ls) %>%
    merge(new_nodes_df_heterodimer %>% dplyr::select(-c("formula", "mass", "rdbe")), all.x = T) %>%
    mutate(node_id = node2) %>%
    dplyr::select(-node2) %>% 
    distinct(formula, node_id, transform, parent_formula, parent_id, .keep_all = T)

  return(heterodimer)
}
### propagate_experiment_MS2_fragment ####
propagate_experiment_MS2_fragment = function(new_nodes_df, sf, 
                                               EdgeSet_experiment_MS2_fragment, 
                                               NodeSet, current_step){
  
  
  if(nrow(EdgeSet_experiment_MS2_fragment)==0){
    return(NULL)
  }
  node_mass = sapply(NodeSet, "[[", "mz")
  
  node1 = EdgeSet_experiment_MS2_fragment$node1
  node2 = EdgeSet_experiment_MS2_fragment$node2
  node1_node2_mapping = bind_rows(EdgeSet_experiment_MS2_fragment) %>%
    dplyr::select(-c("category", "linktype", "direction"))
  
  new_nodes_df_MS2_fragment = new_nodes_df %>% 
    filter(node_id %in% node2) %>% 
    merge(node1_node2_mapping, by.x="node_id", by.y="node2") %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Experiment_MS2_fragment",
           transform = "Experiment_MS2_fragment", 
           direction = -1,
           steps = current_step,
           formula = paste0("MS2_fragment_", formula)) %>%
    mutate(node_id = node1, 
           mass = node_mass[node_id] - mass_dif) %>%
    dplyr::select(-node1, -mass_dif)
  
  return(new_nodes_df_MS2_fragment)
}

### propagate_library_MS2_fragment ####
propagate_library_MS2_fragment = function(new_nodes_df, sf, EdgeSet_library_MS2_fragment_df, 
                                          NodeSet, current_step){
  node_mass = sapply(NodeSet, "[[", "mz")
  
  node1 = EdgeSet_library_MS2_fragment_df$node1
  node2 = EdgeSet_library_MS2_fragment_df$node2
  node1_node2_mapping = EdgeSet_library_MS2_fragment_df %>%
    dplyr::select(-c("category", "linktype", "direction")) %>%
    dplyr::rename(formula = formula2,
                  node_id = node2)
  
  new_nodes_df_library_MS2_fragment = new_nodes_df %>% 
    filter(node_id %in% node2) %>%
    merge(node1_node2_mapping)
  if(nrow(new_nodes_df_library_MS2_fragment) == 0){
    return(NULL)
  }
  # if(current_step == 0.03){browser()}
  
  new_nodes_df_library_MS2_fragment = new_nodes_df_library_MS2_fragment%>% 
    # merge(node1_node2_mapping) %>%
    mutate(parent_id = node_id,
           parent_formula = formula,
           category = "Library_MS2_fragment",
           transform = "Library_MS2_fragment", 
           direction = -1,
           steps = current_step,
           formula = formula1) %>%
    mutate(node_id = node1, 
           mass = formula_mz(formula)) %>%
    dplyr::select(-node1, -formula1,)
  
  return(new_nodes_df_library_MS2_fragment)
}
## Propagate_formulaset ####
Propagate_formulaset = function(Mset, 
                                NodeSet,
                                FormulaSet,
                                EdgeSet_all_df,
                                biotransform_step = 5,
                                artifact_step = 5,
                                propagation_ppm_threshold = 5e-6,
                                propagation_abs_threshold = 2e-4,
                                record_RT_tol = 0.1,
                                record_ppm_tol = 5e-6)
{
  # biotransform_step = 2
  # artifact_step = 3
  # propagation_ppm_threshold = 5e-6
  # propagation_abs_threshold = 2e-4
  # record_RT_tol = 0.1
  # record_ppm_tol = 5e-6

  sf = FormulaSet
  empirical_rules = Mset$Empirical_rules
  node_mass = sapply(NodeSet, "[[", "mz") %>% sort()
  node_RT = sapply(NodeSet, "[[", "RT")[names(node_mass)]
  temp_id = as.numeric(names(node_mass)) # numeric

  
  propagation_rule = Mset$Global_parameter$propagation_rule
  propagation_rule_ls = list()
  for(i in colnames(propagation_rule)){
    propagation_rule_ls[[i]] = colnames(propagation_rule)[propagation_rule[,i]]
  }
  
  ## Expansion 
  timer = Sys.time()
  step_count = 0
  while(step_count < biotransform_step){
    all_nodes_df = bind_rows(sf)
    # Handle artifacts
    sub_step = 0
    while(sub_step < 0.01 * artifact_step){

      node_step = step_count + sub_step
      sub_step = sub_step+0.01
      current_step = step_count + sub_step
      
      all_nodes_df = bind_rows(sf)
      new_nodes_df = all_nodes_df %>%
        distinct(node_id, formula, category, .keep_all = T) %>%
        filter(category != "Unknown") %>%
        filter(abs(node_mass[as.character(node_id)] - mass) < 
                 pmax(propagation_ppm_threshold * mass, propagation_abs_threshold)) # propagate from accurate formulas
      
      if(nrow(new_nodes_df)==0){break}
      
      lib_artifact_ls = list()
      # to Adduct ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Adduct"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step) 
        
        if(nrow(temp_new_nodes_df) != 0){
          temp_rule = empirical_rules %>% filter(category == "Adduct")
          
          temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = 1, category = "Adduct")
          
          lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
        }
      }
      
      # to Fragment ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Fragment"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step) 
        
        if(nrow(temp_new_nodes_df) != 0){
          temp_rule = empirical_rules %>% filter(category == "Fragment")
          
          temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = -1, category = "Fragment")
          
          lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
        }
        
      }
      
      # to Natural_abundance ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Natural_abundance"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step) 
        
        if(nrow(temp_new_nodes_df) != 0){
        
          temp_rule = empirical_rules %>% filter(category == "Natural_abundance") %>% filter(direction == 1)
          temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = 1, category = "Natural_abundance")
          lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
          
          # [10]B has a direction of -1
          temp_rule = empirical_rules %>% filter(category == "Natural_abundance") %>% filter(direction == -1)
          temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = -1, category = "Natural_abundance")
          lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
        }
      }
      
      # to Radical ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Radical"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step)
        
        if(nrow(temp_new_nodes_df) != 0){
          temp_rule = empirical_rules %>% filter(category == "Radical")
          
          temp_propagation = expand_library(temp_new_nodes_df, temp_rule, direction = -1, category = "Radical")
          
          lib_artifact_ls[[length(lib_artifact_ls)+1]] = temp_propagation
        }
      }
      
      # to Heterodimer ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Heterodimer"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step) 
        
        temp_edgeset = EdgeSet_all_df %>%
          filter(category == "Heterodimer")
        
        temp_propagation_category = propagation_rule_ls[["Heterodimer"]]
        heterodimer = propagate_heterodimer(temp_new_nodes_df, sf, 
                                            temp_edgeset, NodeSet, current_step, 
                                            propagation_ppm_threshold, temp_propagation_category)
      }
      
      # to Oligomer ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Oligomer"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step) 
        
        temp_edgeset = EdgeSet_all_df %>%
          filter(category == "Oligomer")
        
        oligomer = propagate_oligomer(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }

      # to Multicharge ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Multicharge"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step)
        
        temp_edgeset = EdgeSet_all_df %>%
          filter(category == "Multicharge")
        
        multicharge = propagate_multicharge(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }
      
      # to Library_MS2_fragment ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Library_MS2_fragment"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step) 
        
        temp_edgeset = EdgeSet_all_df %>%
          filter(category == "Library_MS2_fragment")
        
        library_MS2_fragment = propagate_library_MS2_fragment(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }
      
      # to Experiment_MS2_fragment ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Experiment_MS2_fragment"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step) 
        
        temp_edgeset = EdgeSet_all_df %>%
          filter(category == "Experiment_MS2_fragment")
        
        experiment_MS2_fragment = propagate_experiment_MS2_fragment(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }
      
      # to Ring_artifact ##
      {
        temp_new_nodes_df = new_nodes_df %>%
          filter(category %in% propagation_rule_ls[["Ring_artifact"]]) %>%
          distinct(node_id, formula, .keep_all = T) %>%
          filter(steps==node_step) 
        
        temp_edgeset = EdgeSet_all_df %>%
          filter(category == "Ring_artifact")
        
        ring_artifact = propagate_ring_artifact(temp_new_nodes_df, sf, temp_edgeset, NodeSet, current_step)
      }
      
      lib_artifact = bind_rows(lib_artifact_ls) %>%
        filter(!str_detect(formula, "((?<!H)-)|(-(?!1))")) %>% # any "-" not preceded by H, or not followed by 1 is removed
        # mutate(RT = node_RT[as.character(parent_id)]) %>%
        arrange(mass)
      
      sf_add_mass_filter = bind_rows(oligomer, multicharge, heterodimer, library_MS2_fragment) %>%
        mutate(msr_mass = node_mass[as.character(node_id)],
               mass_dif = abs(msr_mass-mass)) %>%
        mutate(mass_retain = mass_dif / mass <= record_ppm_tol) %>%
        filter(mass_retain) %>%
        dplyr::select(-msr_mass, -mass_dif, -mass_retain, -formula1, -formula2, -edge_id)
      
      sf_add = bind_rows(sf_add_mass_filter, ring_artifact, experiment_MS2_fragment) %>%
        dplyr::select(-formula1, -formula2, -mass_dif, -edge_id)
      
      sf = match_library(lib_artifact,
                         sf,
                         record_ppm_tol,
                         record_RT_tol,
                         current_step,
                         NodeSet)
      
      for(i in unique(sf_add$node_id)){
        sf[[i]] = bind_rows(sf[[i]], sf_add[sf_add$node_id == i, ])
      }
      
      print(paste("Step", current_step, "elapsed="))
      print((Sys.time()-timer))
    }
    
    all_nodes_df = bind_rows(sf)
    new_nodes_df = all_nodes_df %>%
      filter(category == "Metabolite") %>% # only metabolites go to biotransformation, also garantee it is not filtered.
      distinct(node_id, formula, .keep_all=T) %>% # This garantee only new formulas will go to next propagation
      filter(steps == step_count) %>% # only formulas generated from the step go to next propagation
      filter(rdbe > -1) %>% # filter out formula has ring and double bind less than -1
      filter(abs(node_mass[as.character(node_id)] - mass) < pmax(propagation_ppm_threshold * mass, 
                                                                 propagation_abs_threshold)) # propagate from accurate formulas
    
    step_count = step_count+1
    current_step = step_count
    
    if(nrow(new_nodes_df)==0){break}
    
    rule_1 = empirical_rules %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,1))
    rule_2 = empirical_rules %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,-1))
    
    new_nodes_df = new_nodes_df %>%
      mutate(parent_id = node_id)
    lib_1 = expand_library(new_nodes_df, rule_1, direction = 1, category = "Metabolite")
    lib_2 = expand_library(new_nodes_df, rule_2, direction = -1, category = "Metabolite")
    
    lib_met = bind_rows(lib_1, lib_2) %>%
      filter(!grepl("-|NA", formula)) %>% 
      arrange(mass)
    
    sf = match_library(lib_met,
                       sf,
                       record_ppm_tol,
                       record_RT_tol=Inf,
                       current_step,
                       NodeSet)
    
    print(paste("Step",current_step,"elapsed="))
    print(Sys.time()-timer)
  }
  
  return(sf)
}

## Score_formulaset ####
Score_formulaset = function(FormulaSet,
                            database_match = 0.5, 
                            MS2_match = 1,
                            MS2_match_cutoff = 0.8,
                            MS2_similarity = 0.5,
                            MS2_similarity_cutoff = 0.5,
                            rt_match = 1, 
                            known_rt_tol = 0.5,
                            manual_match = 1,
                            bio_decay = 1,
                            artifact_decay = -0.5){
  
  # database_match = 0.5
  # MS2_match = 1
  # MS2_match_cutoff = 0.8
  # MS2_similarity = 0.5
  # MS2_similarity_cutoff = 0.5
  # rt_match = 1
  # known_rt_tol = 0.5
  # manual_match = 1
  # bio_decay = 1
  # artifact_decay = -0.5
  
  FormulaSet_df = bind_rows(FormulaSet) %>%
    mutate(formula_set_id = 1:nrow(.))
  
  # Score HMDB and known adduct match
  {
    FormulaSet_df = FormulaSet_df %>%
      mutate(database_prior = case_when(
        category == "Manual" ~ manual_match,
        category == "Unknown" ~ 0,
        steps == 0 & transform == "" ~ database_match,
        
        TRUE ~ 0 # Everything else
      ))
  }
  
  # Score MS2 match and similarity
  node_MS2_logic = sapply(NodeSet, function(x){
    !is.null(x$MS2)
  })
  if(MS2_match != 0 & any(node_MS2_logic)){
    
    node_MS2 = (1:length(NodeSet))[node_MS2_logic]
    
    MS2_library_external_id = sapply(MS2_library, "[[", 'external_id')
 
    FormulaSet_df_MS2 = FormulaSet_df %>%
      filter(steps == 0) %>%
      # filter(steps == 0, transform == "") %>%
      merge(LibrarySet %>%
              dplyr::select(library_id, note, mass), 
            by.x = "parent_id",
            by.y = "library_id") %>%
      mutate(contain_msr_MS2 = node_id %in% node_MS2) %>%
      mutate(contain_lib_MS2 = note %in% MS2_library_external_id) %>%
      filter(contain_msr_MS2, contain_lib_MS2)
    
    fwd_dp = rev_dp = rep(0, nrow(FormulaSet_df_MS2))
    
    for(i in 1:nrow(FormulaSet_df_MS2)){
      
      H_mass = 1.007825032
      e_mass = 0.000548579
      ion_mode = Mset$Global_parameter$mode
      
      node_id = FormulaSet_df_MS2$node_id[i]
      mz1 = FormulaSet_df_MS2$mass.x[i] + (H_mass-e_mass) * ion_mode
      MS2_msr = NodeSet[[node_id]]$MS2
      colnames(MS2_msr)=c("mz", "inten")
      if(nrow(MS2_msr)==1){
        next
      }
      
      HMDB_id = FormulaSet_df_MS2$note[i]
      mz2 = FormulaSet_df_MS2$mass.y[i] + (H_mass-e_mass) * ion_mode
      lib_MS2_id = which(MS2_library_external_id == HMDB_id)
      MS2_lib = MS2_library[lib_MS2_id]
   
      # skip if all MS2_lib spectra are one-row
      if(all(sapply(MS2_lib, function(x){nrow(x$spectrum)}) == 1)){
        next
      }
      
      if(abs(mz1 - mz2) < max(1e-3, mz2*10e-6)){
        temp_mz_parent = (mz1 + mz2)/2
      } else {
        temp_mz_parent = -Inf
      }
      # Fwd dot product
      spec_score = 0
      for(j in 1:length(MS2_lib)){
        
        spec1_df = as.data.frame(MS2_msr)
        spec2_df = MS2_lib[[j]]$spectrum
        if(nrow(spec2_df)<=1){
          next
        }
        
        spec_merge_df = try(mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3), silent = T)
        if(inherits(spec_merge_df, "try-error")){
          spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3)
        }
        spec_score[j] = Score_merge_MS2(spec_merge_df, mz_parent = temp_mz_parent)
        ## warning or even error occur here if one spectrum have same mz or close mz
        if(is.na(spec_score[j])){
          spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3)
          spec_score[j] = Score_merge_MS2(spec_merge_df, mz_parent = temp_mz_parent)
        }
      }
      # print(spec_score)
      fwd_dp[i] = max(spec_score, na.rm = T)
      
      # Rev dot product
      spec_score = 0
      for(j in 1:length(MS2_lib)){
        spec1_df = as.data.frame(MS2_msr)
        spec2_df = MS2_lib[[j]]$spectrum
        
        spec1_df[,1] = mz1 - spec1_df[,1]
        spec2_df[,1] = mz2 - spec2_df[,1]
        
        spec1_df = spec1_df[nrow(spec1_df):1,]
        spec2_df = spec2_df[nrow(spec2_df):1,]
        
        spec_merge_df = try(mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3), silent = T)
        if(inherits(spec_merge_df, "try-error")){
          spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3)
        }
        spec_score[j] = Score_merge_MS2(spec_merge_df, mz_parent = 0)
        
        ## warning or even error occur here if one spectrum have same mz or close mz
        if(is.na(spec_score[j])){
          spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3)
          spec_score[j] = Score_merge_MS2(spec_merge_df, mz_parent = 0)
        }
      }
      rev_dp[i] = max(spec_score)
    }
    
    # How to combine fwd score and rev score into a MS2 score
    # For exact matched MS2, it takes only fwd spectrum score
    # For transformed MS2, it takes the higher one of the fwd and rev spectrum score
    FormulaSet_df_MS2_ = FormulaSet_df_MS2 %>%
      mutate(fwd_score = fwd_dp,
             rev_score = rev_dp,
             higher_fwd_rev = pmax(fwd_score, rev_score)) %>%
      mutate(MS2score_prior = case_when(
        transform=="" & fwd_score>MS2_match_cutoff ~ fwd_score*MS2_match,
        transform!="" & higher_fwd_rev>MS2_similarity_cutoff ~ higher_fwd_rev*MS2_similarity,
        TRUE ~ 0 # The rest
        ))%>%
      dplyr::select(formula_set_id, MS2score_prior)
    
    FormulaSet_df = FormulaSet_df %>%
      merge(FormulaSet_df_MS2_, all.x = T) %>%
      mutate(MS2score_prior = ifelse(is.na(MS2score_prior), 0, MS2score_prior))
  }
  
  # Score known RT match
  {
    node_RT = sapply(NodeSet, "[[", "RT")
    library_RT = LibrarySet$rt
    names(library_RT) = LibrarySet$library_id
    
    FormulaSet_df = FormulaSet_df %>%
      mutate(node_rt = node_RT[node_id], 
             library_rt = library_RT[as.character(parent_id)]) %>%
      mutate(known_rt_prior = ifelse(abs(node_rt - library_rt) < known_rt_tol & 
                                       steps == 0 & transform == "", rt_match-abs(node_rt - library_rt), 0)) %>%
      mutate(known_rt_prior = replace_na(known_rt_prior, 0)) %>%
      dplyr::select(-node_rt, -library_rt)
  }
  
  # Penalize formula based on PO ratio
  {
    temp_formula = FormulaSet_df$formula
    
    formula_O = str_extract_all(temp_formula, "(?<=O)([:digit:]|\\.)+")
    formula_O = sapply(formula_O, function(x){sum(as.numeric(x))})
    formula_P = str_extract_all(temp_formula, "(?<=P)([:digit:]|\\.)+")
    formula_P = sapply(formula_P, function(x){sum(as.numeric(x))})
    formula_Si = str_extract_all(temp_formula, "(?<=Si)([:digit:]|\\.)+")
    formula_Si = sapply(formula_Si, function(x){sum(as.numeric(x))})
    
    FormulaSet_df = FormulaSet_df %>%
      mutate(empirical_POratio_prior = ifelse(formula_O >= 3*formula_P & 
                                                formula_O >= 1.01*formula_Si, 0, -10)) # make sure when no O and no Si, it is true
  }
  
  # Penalize formula based on RDBE rule
  {
    FormulaSet_df = FormulaSet_df %>%
      mutate(empirical_RDBE_prior = ifelse(rdbe >= 0 | category != "Metabolite", 0, -10))
  }
  
  
  # Prior formula scores 
  {
    prior_formula_score = FormulaSet_df %>% 
      dplyr::select(ends_with("_prior")) %>% 
      rowSums(na.rm = T)
    
    FormulaSet_df = FormulaSet_df %>%
      mutate(prior_score = prior_formula_score)
  }
  
  return(FormulaSet_df)
}


## initiate_ilp_nodes ####
initiate_ilp_nodes = function(FormulaSet_df,
                              artifact_decay = 1){

  # artifact_decay = 1
  
  ilp_nodes = FormulaSet_df %>%
    mutate(class = ifelse(category %in% c("Metabolite", "Unknown"), category, "Artifact")) %>%
    distinct(node_id, formula, class, .keep_all=T) %>%
    arrange(node_id) %>%
    mutate(ilp_node_id = 1:nrow(.)) %>%
    dplyr::select(ilp_node_id, node_id, class, formula) %>%
    filter(T)
  
  # Unknowns 
  {
    FormulaSet_unknowns = FormulaSet_df %>% 
      filter(category == "Unknown")
    ilp_nodes_unknowns = ilp_nodes %>%
      filter(class == "Unknown") %>%
      merge(FormulaSet_unknowns) %>%
      mutate(score_prior_propagation = prior_score)
  }
  
  # Metabolites
  {
    FormulaSet_metabolites = FormulaSet_df %>% 
      filter(category == "Metabolite")
    ilp_nodes_metabolites = ilp_nodes %>%
      filter(class == "Metabolite") %>%
      merge(FormulaSet_metabolites) %>%
      mutate(score_prior_propagation = prior_score) %>%
      arrange(-score_prior_propagation, steps) %>%
      distinct(node_id, formula, .keep_all = T)
  }
  
  # Artifact steps 0
  {
    FormulaSet_artifacts0 = FormulaSet_df %>% 
      filter(!category %in% c("Metabolite", "Unknown", "Heterodimer")) %>%
      filter(steps == 0)
    ilp_nodes_artifacts0 = ilp_nodes %>%
      filter(class == "Artifact") %>%
      merge(FormulaSet_artifacts0, all.y = T) %>%
      mutate(score_prior_propagation = prior_score) %>%
      arrange(-score_prior_propagation, steps) %>%
      distinct(node_id, formula, .keep_all = T)
  }
  
  ilp_nodes2 = bind_rows(ilp_nodes_unknowns, ilp_nodes_metabolites, ilp_nodes_artifacts0)
  
  # Artifacts
  {
    
    FormulaSet_artifacts = FormulaSet_df %>% 
      # filter(!category %in% c("Metabolite", "Unknown", "Heterodimer"))
      filter(!category %in% c("Metabolite", "Unknown")) %>%
      filter(steps != 0)
    
    floating_error = 1e-6
    ilp_nodes_artifacts = ilp_nodes2
    unique_step = unique(FormulaSet_artifacts$steps)[1]
    for(unique_step in sort(unique(FormulaSet_artifacts$steps))){
      temp_ilp_nodes = ilp_nodes_artifacts %>%
        filter(abs(steps - (unique_step - 0.01)) < floating_error) %>%
        arrange(-score_prior_propagation, category) %>%
        dplyr::select(node_id, formula, score_prior_propagation) %>%
        distinct(node_id, formula, .keep_all = T) %>%
        dplyr::rename(parent_id = node_id, 
                      parent_formula = formula)
        
      
      temp = FormulaSet_artifacts %>%
        filter(steps == unique_step) %>%
        merge(temp_ilp_nodes) %>%
        mutate(score_prior_propagation = case_when(
          category %in% c("Library_MS2_fragment", "Ring_artifact", "Natural_abundance") ~ 0,
          score_prior_propagation <= 0 ~ score_prior_propagation + prior_score,  
          score_prior_propagation < artifact_decay ~ prior_score,
          score_prior_propagation >= artifact_decay ~ score_prior_propagation + prior_score - artifact_decay
        )) %>%
        mutate(class = "Artifact")
      
      ilp_nodes_artifacts = bind_rows(ilp_nodes_artifacts, temp)
    }

    ilp_nodes_artifacts2 = ilp_nodes_artifacts %>%
      filter(is.na(ilp_node_id)) %>%
      dplyr::select(-ilp_node_id) %>%
      merge(ilp_nodes)


  }
  
  ilp_nodes_result = bind_rows(ilp_nodes2, ilp_nodes_artifacts2) %>%
    arrange(node_id, -score_prior_propagation, steps, category) %>%
    distinct(ilp_node_id, formula, class, .keep_all = T) %>%
    # group_by(ilp_node_id) %>%
    # filter(n()>1)
    arrange(ilp_node_id)
  
  return(ilp_nodes_result)
  
}

## score_ilp_nodes ####
score_ilp_nodes = function(CplexSet, mass_dist_gamma_rate = mass_dist_gamma_rate, 
                           formula_score = 0, metabolite_score = 0.1, 
                           isotope_Cl_penalty = -1, unassigned_penalty = -0.5){
  
  ilp_nodes = CplexSet$ilp_nodes
  ilp_edges = CplexSet$ilp_edges
  
  node_mass = sapply(NodeSet, "[[", "mz")
  # Formula and metabolite score
  {
    ilp_nodes = ilp_nodes %>%
      mutate(msr_mass = node_mass[node_id],
             ppm_error = (mass - msr_mass) / mass * 1e6) %>%
      mutate(score_formula = formula_score) %>%
      mutate(score_formula = ifelse(category == "Metabolite", score_formula + metabolite_score, score_formula))
  }
  

  # Score mass accuracy
  {
    ilp_nodes = ilp_nodes %>%
      mutate(score_mass = dgamma(abs(ppm_error), 1, scale = mass_dist_gamma_rate) * mass_dist_gamma_rate) %>%
      mutate(score_mass = log10(score_mass)) %>%
      filter(T)
  }

  # Score penalties for isotopes
  if(isotope_Cl_penalty != 0){
    node_inten = sapply(NodeSet,"[[","inten")
    ilp_edges_Cl = ilp_edges %>% filter(grepl("Cl-1", linktype))
    ilp_nodes_Cl = ilp_nodes %>%
      filter(str_detect(formula, "(?<!\\])Cl")) %>% # Find Cl but not Cl isotope formulas
      mutate(score_Cl_isotope = ifelse(ilp_node_id %in% ilp_edges_Cl$ilp_nodes1 &
                                         node_inten[node_id] > 4.5, 
                                       0, isotope_Cl_penalty)) %>%
      dplyr::select(ilp_node_id, score_Cl_isotope)
    
    # isotope_Cr_penalty = isotope_Cl_penalty
    # ilp_edges_Cr = ilp_edges %>% filter(grepl("Cr-1", linktype))
    # ilp_nodes_Cr = ilp_nodes %>%
    #   filter(str_detect(formula, "(?<!\\])Cr")) %>%
    #   mutate(score_Cr_isotope = ifelse(ilp_node_id %in% ilp_edges_Cr$ilp_nodes1 &
    #                                      node_inten[node_id] > 5, 
    #                                    0, isotope_Cr_penalty)) %>%
    #   dplyr::select(ilp_node_id, score_Cr_isotope)
    ilp_nodes = ilp_nodes %>%
      merge(ilp_nodes_Cl, all.x = T) %>%
      mutate(score_Cl_isotope = ifelse(is.na(score_Cl_isotope), 0, score_Cl_isotope))
  }
  
  # Score penalties for unassigned peaks
  if(unassigned_penalty != 0){
    nodes_inten = sapply(NodeSet, "[[", "inten")
    ilp_nodes_unassigned = ilp_nodes %>%
      filter(category == "Unknown") %>%
      mutate(inten = nodes_inten[node_id]) %>% 
      mutate(score_unassigned = ifelse(inten < 5, 0, unassigned_penalty * (inten-5))) %>%
      dplyr::select(ilp_node_id, score_unassigned)
    ilp_nodes = ilp_nodes %>%
      merge(ilp_nodes_unassigned, all.x = T) %>%
      mutate(score_unassigned = ifelse(is.na(score_unassigned), 0, score_unassigned))
  }
  
  
  
  cplex_score_node = ilp_nodes %>% 
    dplyr::select(starts_with("score_")) %>% 
    rowSums()
  
  ilp_nodes = ilp_nodes %>%
    mutate(cplex_score = cplex_score_node)
  
  return(ilp_nodes)
}

## initiate_ilp_edges ####
initiate_ilp_edges = function(EdgeSet_all_df, CplexSet, Exclude = "Biotransform"){
  
  ilp_nodes = CplexSet$ilp_nodes %>%
    arrange(ilp_node_id) 
  
  EdgeSet_df = EdgeSet_all_df %>%
    filter(category != "Heterodimer") %>%
    filter(!category %in% Exclude)
  
  ilp_nodes_assigned = ilp_nodes %>% filter(category != "Unknown")
  node_id_ilp_node_mapping = split(ilp_nodes_assigned$ilp_node_id, ilp_nodes_assigned$node_id)
  
  ilp_nodes_met = ilp_nodes %>% filter(category == "Metabolite") 
  node_id_ilp_node_mapping_met = split(ilp_nodes_met$ilp_node_id, ilp_nodes_met$node_id)
  
  ilp_nodes_nonmet = ilp_nodes %>% filter(!category %in% c("Metabolite","Unknown"))
  node_id_ilp_node_mapping_nonmet = split(ilp_nodes_nonmet$ilp_node_id, ilp_nodes_nonmet$node_id)
  
  ilp_nodes_radical = ilp_nodes %>% filter(category == "Radical")
  node_id_ilp_node_mapping_radical = split(ilp_nodes_radical$ilp_node_id, ilp_nodes_radical$node_id)
  
  ilp_nodes_formula = ilp_nodes %>%
    pull(formula)
  
  match_matrix_index_ls = list()
  
  for(i in 1:nrow(EdgeSet_df)){
    # print(i)
    edge_id = EdgeSet_df$edge_id[i]
    category = EdgeSet_df$category[i]
    
    if(category == "Biotransform"){
      # Only allow two metabolite ilp_node to form biotransform connection
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping_met[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping_met[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      
      transform = EdgeSet_df$linktype[i]
      formula_transform = my_calculate_formula(formula1, transform)
    #------------------------------------------------------------------------#
    } else if(category == "Fragment") {
      # node1 cannot be metabolite, and node2 can be anything
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping_nonmet[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      
      transform = EdgeSet_df$linktype[i]
      formula_transform = my_calculate_formula(formula1, transform)
    #------------------------------------------------------------------------#
    } else if(category == "Oligomer") {
      # Both can be anything for simplicity
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping_nonmet[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      
      fold = as.numeric(EdgeSet_df$linktype[i])
      formula_transform = mapply(my_calculate_formula, formula1, formula1, fold-1)
    #------------------------------------------------------------------------#
    } else if(category == "Multicharge") {
      # Both can be anything for simplicity
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping_nonmet[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      
      fold = as.numeric(EdgeSet_df$linktype[i])
      formula_transform = mapply(my_calculate_formula, formula1, formula1, fold-1)
      #------------------------------------------------------------------------#
    } else if(category == "Experiment_MS2_fragment"){
      # node1 cannot be metabolite, and node2 to be anything
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping_nonmet[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      
      formula_transform = sub("MS2_fragment_", "", formula1)
    } else if(category == "Library_MS2_fragment"){
      # node2 has to be metabolite because of HMDB library matching
      
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping_nonmet[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping_met[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      
      formula1_position = which(formula1 == EdgeSet_df$formula1[i])
      formula_transform = rep("", length(formula1))
      formula_transform[formula1_position] = EdgeSet_df$formula2[i]
    } else if(category == "Radical"){
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping_radical[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      transform = EdgeSet_df$linktype[i]
      formula_transform = my_calculate_formula(formula1, transform)
    } else {
      # if not above situation, then artifact
      # node2 has to be artifact, node1 can be anything
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping_nonmet[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      
      if(any(category == c("Adduct", "Natural_abundance", "Multicharge_isotope"))){
        # Direction is not needed as it is always node1 + linktype = node2
        transform = EdgeSet_df$linktype[i]
        formula_transform = my_calculate_formula(formula1, transform)
      } else if(category == c("Ring_artifact")){
        formula_transform = paste0("Ring_artifact_", formula1)
      } 
    }
    # Matrix can parallel calculate all formulas of node1
    # And compares with all formulas of node2

    formula_match_matrix = matrix(TRUE, length(formula_transform), length(formula2))
    for(j in 1:length(formula_transform)){
      formula_match_matrix[j,] = formula_transform[j] == formula2
    }
    
    formula_match_matrix_index = which(formula_match_matrix == T, arr.ind = TRUE) 
    if(dim(formula_match_matrix_index)[1] == 0){next}
    

    match_matrix_index_ls[[length(match_matrix_index_ls)+1]] = list(edge_id = rep(edge_id, dim(formula_match_matrix_index)[1]),
                                      ilp_nodes1 = ilp_nodes1[formula_match_matrix_index[, 1]],
                                      ilp_nodes2 = ilp_nodes2[formula_match_matrix_index[, 2]])
  }
  
  
  
  ilp_edges = bind_rows(match_matrix_index_ls) %>%
    merge(EdgeSet_df, all.x = T) %>%
    mutate(formula1 = ilp_nodes_formula[ilp_nodes1],
           formula2 = ilp_nodes_formula[ilp_nodes2]) %>%
    mutate(ilp_edge_id = 1:nrow(.))
  
  return(ilp_edges)
}

## initiate_heterodimer_ilp_edges ####
initiate_heterodimer_ilp_edges = function(EdgeSet_all_df, CplexSet, NodeSet){
  
  node_inten = sapply(NodeSet, "[[", "inten")
  EdgeSet_df = EdgeSet_all_df %>%
    filter(category == "Heterodimer") %>%
    mutate(linktype = as.numeric(linktype)) %>%
    filter(node_inten[node1] > node_inten[linktype]) # node1 and linktype are interchangeable in heterodimer

  if(nrow(EdgeSet_df) == 0){
    return(NULL)
  }
  
  ilp_nodes = CplexSet$ilp_nodes %>%
    arrange(ilp_node_id) 

  ilp_nodes_assigned = ilp_nodes %>% filter(category != "Unknown")
  node_id_ilp_node_mapping = split(ilp_nodes_assigned$ilp_node_id, ilp_nodes_assigned$node_id)
  
  ilp_nodes_nonmet = ilp_nodes %>% 
    filter(!category %in% c("Metabolite","Unknown")) %>%
    filter(rdbe %% 1 == 0) %>% # Remove heterodimer which may be radicals
    filter(!grepl("\\[", formula)) # Remove natural abundance in heterodimer
    
  node_id_ilp_node_mapping_nonmet = split(ilp_nodes_nonmet$ilp_node_id, ilp_nodes_nonmet$node_id)
  
  ilp_nodes_formula = ilp_nodes %>%
    pull(formula)

  match_matrix_index_ls = list()
  
  for(i in 1:nrow(EdgeSet_df)){
    edge_id = EdgeSet_df$edge_id[i]
    
    # node1 and node_link can be anything, node2 has to be artifact
    node1 = EdgeSet_df$node1[i]
    ilp_nodes1 = node_id_ilp_node_mapping[[as.character(node1)]]
    formula1 = ilp_nodes_formula[ilp_nodes1]
    if(length(formula1) == 0){next}
    
    node2 = EdgeSet_df$node2[i]
    ilp_nodes2 = node_id_ilp_node_mapping_nonmet[[as.character(node2)]]
    formula2 = ilp_nodes_formula[ilp_nodes2]
    if(length(formula2) == 0){next}
    
    node_link = EdgeSet_df$linktype[i]
    ilp_nodes_link = node_id_ilp_node_mapping[[as.character(node_link)]]
    formula_link = ilp_nodes_formula[ilp_nodes_link]
    if(length(formula_link) == 0){next}
    
    # if(length(formula_link)>1 & length(formula1)>1)break
    transform = formula_link
    formula_transform_all = my_calculate_formula(formula1, transform)
    
    for(k in 1:length(transform)){
      if(length(transform) == 1){
        formula_transform = formula_transform_all
      } else {
        formula_transform = formula_transform_all[, k]      
      }
      
      formula_match_matrix = matrix(TRUE, length(formula_transform), length(formula2))
      for(j in 1:length(formula_transform)){
        formula_match_matrix[j,] = formula_transform[j] == formula2
      }
      
      formula_match_matrix_index = which(formula_match_matrix == T, arr.ind = TRUE) 
      
      if(dim(formula_match_matrix_index)[1] == 0){next}
      
      match_matrix_index_ls[[length(match_matrix_index_ls)+1]] = list(edge_id = rep(edge_id, dim(formula_match_matrix_index)[1]),
                                                                      ilp_nodes1 = ilp_nodes1[formula_match_matrix_index[, 1]],
                                                                      ilp_nodes2 = ilp_nodes2[formula_match_matrix_index[, 2]],
                                                                      ilp_nodes_link = rep(ilp_nodes_link[k], dim(formula_match_matrix_index)[1]))
    }
  }
  heterodimer_ilp_edges = bind_rows(match_matrix_index_ls) %>%
    merge(EdgeSet_df, all.x = T) %>%
    mutate(formula1 = ilp_nodes_formula[ilp_nodes1],
           formula2 = ilp_nodes_formula[ilp_nodes2],
           formula_link = ilp_nodes_formula[ilp_nodes_link])
  return(heterodimer_ilp_edges)
}



## score_ilp_edges ####
score_ilp_edges = function(CplexSet, NodeSet,
                           rule_score_biotransform = 0, rule_score_adduct = 0.5, 
                           rule_score_oligomer = 0.5, rule_score_natural_abundance = 1,
                           rule_score_fragment = 0.3, rule_score_ring_artifact = 2,
                           rule_score_radical = 0.2, rule_score_multicharge = 0.5, rule_score_multicharge_isotope = 1,
                           rule_score_experiment_MS2_fragment = 1, rule_score_library_MS2_fragment = 0.3,
                           rt_penalty_artifact_cutoff = 0.05, rt_penalty_artifact_ratio = 5,
                           inten_score_isotope = 1, 
                           MS2_score_similarity = 1, MS2_similarity_cutoff = 0.3,
                           MS2_score_experiment_fragment = 0.5){

  # rule_score_biotransform = 0
  # rule_score_adduct = 0.5
  # rule_score_oligomer = 0.5
  # rule_score_natural_abundance = 1
  # rule_score_fragment = 0.3
  # rule_score_ring_artifact = 2
  # rule_score_radical = 0.2
  # rule_score_experiment_MS2_fragment = 1
  # rule_score_library_MS2_fragment = 0.3
  # rt_penalty_artifact_cutoff = 0.05
  # rt_penalty_artifact_ratio = 5
  # inten_score_isotope = 1
  # MS2_score_similarity = 1
  # MS2_similarity_cutoff = 0.3
  # MS2_score_experiment_fragment = 0.5
  
  
  ilp_edges = CplexSet$ilp_edges
    
  # Score rule category 
  {
    ilp_edges = ilp_edges %>%
      mutate(score_category = case_when(
        category == "Biotransform" ~ rule_score_biotransform,
        category == "Ring_artifact" ~ rule_score_ring_artifact,
        category == "Experiment_MS2_fragment" ~ rule_score_experiment_MS2_fragment,
        category == "Library_MS2_fragment" ~ 0, # special handling below
        category == "Natural_abundance" ~ rule_score_natural_abundance,
        category == "Oligomer" ~ rule_score_oligomer,
        category == "Multicharge" ~ rule_score_multicharge, 
        category == "Fragment" ~ rule_score_fragment,
        category == "Adduct" ~ rule_score_adduct, 
        category == "Radical" ~ rule_score_radical,
        category == "Multicharge_isotope" ~ rule_score_multicharge_isotope
        # category != "Biotransform" ~ rule_score_artifact, 
      )) %>%
      filter(T)
  }
  
  # Penalty for large RT dif
  # This should be replaced by chormatogram correlation   
  if(rt_penalty_artifact_ratio != 0){
    node_rt = sapply(NodeSet, "[[", "RT")
    
    ilp_edges_rt = ilp_edges %>%
      mutate(rt_node1 = node_rt[node1],
             rt_node2 = node_rt[node2]) %>%
      mutate(rt_dif = rt_node2 - rt_node1) %>%
      mutate(score_rt_artifact = case_when(
        category %in% c("Biotransform") ~ 0,
        abs(rt_dif) < rt_penalty_artifact_cutoff ~ 0,
        T ~ (rt_penalty_artifact_cutoff - abs(rt_dif)) * rt_penalty_artifact_ratio)) %>%
      dplyr::select(ilp_edge_id, score_rt_artifact) %>%
      filter(T)
    ilp_edges = ilp_edges %>%
      merge(ilp_edges_rt, all.x = T) %>%
      mutate(score_rt_artifact = ifelse(is.na(score_rt_artifact), 0, score_rt_artifact))
  }
  
  # Score library_MS2_fragment
  # Describe if given any peak related to a HMDB metabolite, 
  # look for the metabolite's MS2, to see if any fragment is presented in MS1
  # If so, strengthen the formula pair of the fragment and parent
  if(rule_score_library_MS2_fragment != 0){
    ilp_edges_library_MS2_fragment = ilp_edges %>%
      filter(category == "Library_MS2_fragment") %>%
      group_by(node2) %>%
      mutate(n = n()) %>%
      mutate(score_library_MS2_fragment = rule_score_library_MS2_fragment/n) %>%
      merge(ilp_edges) %>%
      dplyr::select(ilp_edge_id, score_library_MS2_fragment) %>%
      filter(T)
    ilp_edges = ilp_edges %>%
      merge(ilp_edges_library_MS2_fragment, all.x = T) %>%
      mutate(score_library_MS2_fragment = ifelse(is.na(score_library_MS2_fragment), 0, score_library_MS2_fragment))
  }
  # Score isotope intensity 
  if(inten_score_isotope != 0){
    node_inten = sapply(NodeSet, "[[", "inten")
    ilp_edges_isotope = ilp_edges %>%
      filter(category == "Natural_abundance") %>%
      mutate(inten1 = node_inten[node1],
             inten2 = node_inten[node2]) %>%
      mutate(inten_ratio_measured = inten2 - inten1, 
             inten_ratio_calculated = log10(mapply(isotopic_abundance, .$formula1, .$linktype)),
             measured_calculated_ratio = 10^(inten_ratio_measured-inten_ratio_calculated)) %>%
      mutate(p_obs = dnorm(measured_calculated_ratio, 1, 0.2+10^(3-pmin(inten1, inten2))),
             p_theory = dnorm(1, 1, 0.2+10^(3-pmin(inten1, inten2)))) %>%
      mutate(score_inten_isotope = log10(p_obs/p_theory+1e-10) + inten_score_isotope) %>%
      dplyr::select(ilp_edge_id, score_inten_isotope)
    
    ilp_edges = ilp_edges %>%
      merge(ilp_edges_isotope, all.x = T) %>%
      mutate(score_inten_isotope = ifelse(is.na(score_inten_isotope), 0, score_inten_isotope))
    # hist(ilp_edges_isotope %>% filter(score_inten_isotope<1 & score_inten_isotope>-8 ) %>% pull(score_inten_isotope))
  }
  
  
  # Score Experiment MS2 similarity
  # It compares any connection in edge list where both experiment MS2 exist
  # If their MS2 are similar, then strengthen such connection
  node_MS2_logic = sapply(NodeSet, function(x){
    !is.null(x$MS2)
  })
  if(MS2_score_similarity != 0 & any(node_MS2_logic)){

    node_MS2 = (1:length(NodeSet))[node_MS2_logic]
    
    ilp_edges_MS2_similarity = ilp_edges %>%
      distinct(node1, node2, .keep_all = T) %>%
      filter(node1 %in% node_MS2, node2 %in% node_MS2)
    
    fwd_dp = rev_dp = rep(0, nrow(ilp_edges_MS2_similarity))
    for(i in 1:nrow(ilp_edges_MS2_similarity)){
      
      H_mass = 1.007825032
      e_mass = 0.000548579
      ion_mode = Mset$Global_parameter$mode

      node1 = ilp_edges_MS2_similarity$node1[i]
      node2 = ilp_edges_MS2_similarity$node2[i]
      mz1 = NodeSet[[node1]]$mz + (H_mass - e_mass) * ion_mode
      mz2 = NodeSet[[node2]]$mz + (H_mass - e_mass) * ion_mode
      
      spec1_df = as.data.frame(NodeSet[[node1]]$MS2)
      spec2_df = as.data.frame(NodeSet[[node2]]$MS2)
      
      if(nrow(spec1_df)==1 | nrow(spec2_df)==1){
        next
      }
      
      if(abs(mz1 - mz2) < max(1e-3, mz2*10e-6)){
        temp_mz_parent = (mz1 + mz2)/2
      } else {
        temp_mz_parent = -Inf
      }
      
      colnames(spec1_df)=colnames(spec2_df)=c("mz", "inten")
      
      # Fwd dot product
      spec_merge_df = try(mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3), silent = T)
      if(inherits(spec_merge_df, "try-error")){
        spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3)
      }
      fwd_dp[i] = Score_merge_MS2(spec_merge_df, mz_parent = temp_mz_parent)
      ## warning or even error occur here if one spectrum have same mz or close mz
      if(is.na(fwd_dp[i])){
        spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3)
        fwd_dp[i] = Score_merge_MS2(spec_merge_df, mz_parent = temp_mz_parent)
      }
      
      # Rev dot product
      spec1_df[,1] = mz1 - spec1_df[,1]
      spec2_df[,1] = mz2 - spec2_df[,1]
      spec1_df = spec1_df[nrow(spec1_df):1,]
      spec2_df = spec2_df[nrow(spec2_df):1,]
      spec_merge_df = try(mergeMzIntensity(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3), silent = T)
      if(inherits(spec_merge_df, "try-error")){
        spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3)
      }
      rev_dp[i] = Score_merge_MS2(spec_merge_df, mz_parent = 0)
      if(is.na(rev_dp[i])){
        spec_merge_df = mergeMzIntensity_backup(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3)
        rev_dp[i] = Score_merge_MS2(spec_merge_df, mz_parent = temp_mz_parent)
      }
    }

    ilp_edges_MS2_similarity_ = ilp_edges_MS2_similarity %>%
      mutate(fwd_score = fwd_dp,
             rev_score = rev_dp,
             higher_fwd_rev = pmax(fwd_score, rev_score)) %>%
      filter(higher_fwd_rev>MS2_similarity_cutoff) %>%
      group_by(node1) %>% mutate(n1 = n()) %>% ungroup() %>%
      group_by(node2) %>% mutate(n2 = n()) %>% ungroup() %>%
      mutate(higher_n = pmax(n1, n2)) %>%
      mutate(score_MS2_similarity =  higher_fwd_rev*MS2_score_similarity/higher_n)%>%
      dplyr::select(node1, node2, score_MS2_similarity)
    
    ilp_edges = ilp_edges %>%
      merge(ilp_edges_MS2_similarity_, all.x = T) %>%
      mutate(score_MS2_similarity = ifelse(is.na(score_MS2_similarity), 0, score_MS2_similarity))
  }
  
  # Score MS2 fragment confidence to adducts/heterodimer/oligomer/fragment etc.
  # It 
  if(MS2_score_experiment_fragment != 0){
    ilp_edges_experiment_MS2_fragment = ilp_edges %>%
      filter(category == "Experiment_MS2_fragment") %>%
      distinct(node1, node2) %>%
      merge(ilp_edges) %>%
      filter(category != "Experiment_MS2_fragment") %>%
      # Because biotransform is unlikely to have same RT and show fragment in a MS2 spectrum
      filter(category != "Biotransform") %>% 
      mutate(score_experiment_MS2_fragment = MS2_score_experiment_fragment) %>%
      dplyr::select(ilp_edge_id, score_experiment_MS2_fragment) %>%
      filter(T)
    
    ilp_edges = ilp_edges %>%
      merge(ilp_edges_experiment_MS2_fragment, all.x = T) %>%
      mutate(score_experiment_MS2_fragment = ifelse(is.na(score_experiment_MS2_fragment), 
                                           0, 
                                           score_experiment_MS2_fragment))
  }
  
  
  cplex_score_edge = ilp_edges %>% 
    dplyr::select(starts_with("score_")) %>% 
    rowSums(na.rm=T)
  
  ilp_edges = ilp_edges %>%
    mutate(cplex_score = cplex_score_edge)
  
  return(ilp_edges)
}

## score_heterodimer_ilp_edges ####
score_heterodimer_ilp_edges = function(CplexSet, rule_score_heterodimer = 1,
                                       MS2_score_experiment_fragment = 0.5){
  
  heterodimer_ilp_edges = CplexSet$heterodimer_ilp_edges %>%
    mutate(ilp_edge_id = 1:nrow(.) + nrow(CplexSet$ilp_edges)) %>%
    mutate(score_category = case_when(
      category == "Heterodimer" ~ rule_score_heterodimer,
      category != "Heterodimer" ~ 0
    )) %>%
    filter(T)
  
  if(MS2_score_experiment_fragment != 0){
    heterodimer_ilp_edges_experiment_MS2_fragment = CplexSet$ilp_edges %>%
      filter(category == "Experiment_MS2_fragment") %>%
      distinct(node1, node2) %>%
      merge(heterodimer_ilp_edges) %>%
      mutate(score_experiment_MS2_fragment = MS2_score_experiment_fragment) %>%
      dplyr::select(ilp_edge_id, score_experiment_MS2_fragment) %>%
      filter(T)
    
    heterodimer_ilp_edges = heterodimer_ilp_edges %>%
      merge(heterodimer_ilp_edges_experiment_MS2_fragment, all.x = T) %>%
      mutate(score_experiment_MS2_fragment = ifelse(is.na(score_experiment_MS2_fragment), 
                                                    0, 
                                                    score_experiment_MS2_fragment))
    
  }
  
  cplex_score_edge = heterodimer_ilp_edges %>% 
    dplyr::select(starts_with("score_")) %>% 
    rowSums(na.rm=T)
  
  heterodimer_ilp_edges = heterodimer_ilp_edges %>%
    mutate(cplex_score = cplex_score_edge)
  
  return(heterodimer_ilp_edges)
}

## Prepare_CPLEX parameter ####
Initiate_cplexset = function(CplexSet){
  
  ilp_nodes = CplexSet$ilp_nodes %>%
    arrange(ilp_node_id)
  ilp_edges = CplexSet$ilp_edges %>%
    arrange(ilp_edge_id)
  
  ## Core codes
  # Construct constraint matrix 
  # triplet_nodes
  {
    ilp_rows = rep(1,length = nrow(ilp_nodes))
    ilp_nodes.node_id = ilp_nodes$node_id
    for(i in 2:length(ilp_nodes.node_id)){
      if(ilp_nodes.node_id[i] == ilp_nodes.node_id[i-1]){
        ilp_rows[i] = ilp_rows[i-1]
      } else{
        ilp_rows[i] = ilp_rows[i-1] + 1
      }
    }
    
    triplet_node = ilp_nodes %>%
      mutate(ilp_row_id = ilp_rows) %>%
      dplyr::select(ilp_row_id, ilp_node_id) %>%
      transmute(i = ilp_row_id,  ## Caution: row number is discontinous as not all node_id exist
                j = ilp_node_id, 
                v = 1)
  }
  
  # triplet_edges
  # Because triplet_node take max(triplet_node$i) rows, and max(triplet_node$j) columns
  {
    triplet_edge_edge = ilp_edges %>%
      transmute(i = ilp_edge_id + max(triplet_node$i), 
                j = ilp_edge_id + max(triplet_node$j),
                v = 2)
    
    triplet_edge_node1 = ilp_edges %>%
      transmute(i = ilp_edge_id + max(triplet_node$i),
                j = ilp_nodes1,
                v = -1)
    
    triplet_edge_node2 = ilp_edges %>%
      transmute(i = ilp_edge_id + max(triplet_node$i),
                j = ilp_nodes2,
                v = -1)
    
    triplet_edge = bind_rows(triplet_edge_edge, 
                             triplet_edge_node1,
                             triplet_edge_node2)
  }
  
  # triplet_isotope
  # constrain an isotope formula must come with an isotope edge
  {
    ilp_edges_isotope = ilp_edges %>%
      filter(category == "Natural_abundance") %>%
      group_by(edge_id, ilp_nodes1) %>%
      mutate(n1 = n()) %>%
      ungroup() %>%
      group_by(edge_id, ilp_nodes2) %>%
      mutate(n2 = n()) %>%
      ungroup()
    
    ilp_edges_isotope_1 = ilp_edges_isotope %>%
      filter(n1 == 1, n2 == 1)
    
    triplet_isotope_edge_1 = ilp_edges_isotope_1 %>%
      transmute(i = 1:nrow(.) + max(triplet_edge$i), 
                j = ilp_edge_id + max(triplet_node$j),
                v = 1)
    triplet_isotope_node_1 = ilp_edges_isotope_1 %>%
      transmute(i = 1:nrow(.) + max(triplet_edge$i), 
                j = ifelse(direction == 1, ilp_nodes2, ilp_nodes1),
                v = -1)
    
    ## When two ilp_node1 point to same isotope ilp_node2, two lines needs to be combined
    
    ## Try catching bug here
    {
      # test = ilp_nodes %>%
      #   group_by(node_id, formula) %>%
      #   filter(n()>2)
      # Bug1: one edge contains two same ilp_nodes
      ilp_edges_isotope_bug1 = ilp_edges_isotope %>%
        filter(n1 != 1 & n2 != 1)
      # Bug2: the direction of isotope edge causes problem
      ilp_edges_isotope_bug2.1 = ilp_edges_isotope %>%
        filter(n1 == 1 & n2 != 1 & direction != 1)
      ilp_edges_isotope_bug2.2 = ilp_edges_isotope %>%
        filter(n2 == 1 & n1 != 1 & direction != -1)
      if(nrow(ilp_edges_isotope_bug1) != 0 |
         nrow(ilp_edges_isotope_bug2.1) != 0 |
         nrow(ilp_edges_isotope_bug2.2) != 0){
        stop("Bugs in isotope triplex.")
      }
    }
    
    ilp_edges_isotope_2 = ilp_edges_isotope %>%
      filter(n1 != 1 | n2 != 1) 
    
    ilp_edges_isotope_2.1 = ilp_edges_isotope %>%
      filter(n1 != 1) %>%
      group_by(edge_id, ilp_nodes1) %>%
      group_split()
    ilp_edges_isotope_2.2 = ilp_edges_isotope %>%
      filter(n2 != 1) %>%
      group_by(edge_id, ilp_nodes2) %>%
      group_split()

    # It is assumed that isotope_id is counting the number of list
    ilp_edges_isotope_2 = bind_rows(ilp_edges_isotope_2.1, 
                                    ilp_edges_isotope_2.2, 
                                    .id = "isotope_id") %>%
      mutate(isotope_id = as.numeric(isotope_id))
    
    # isotope_id specify which row the triplex should be in
    
    triplet_isotope_edge_2 = ilp_edges_isotope_2 %>%
      transmute(i = isotope_id + max(triplet_isotope_edge_1$i), 
                j = ilp_edge_id + max(triplet_node$j),
                v = 1)
    triplet_isotope_node_2 = ilp_edges_isotope_2 %>%
      transmute(i = isotope_id + max(triplet_isotope_edge_1$i), 
                j = ifelse(direction == 1, ilp_nodes2, ilp_nodes1),
                v = -1) %>%
      distinct()
    
    nrow_triplet_isotope = nrow(ilp_edges_isotope_1) + max(ilp_edges_isotope_2$isotope_id)
    

    triplet_isotope = bind_rows(triplet_isotope_edge_1, 
                                triplet_isotope_node_1,
                                triplet_isotope_edge_2,
                                triplet_isotope_node_2)
  }
  
  # triplet_heterodimer
  {
    heterodimer_ilp_edges = CplexSet$heterodimer_ilp_edges %>%
      filter(category == "Heterodimer") %>%
      mutate(ilp_edge_id = 1:nrow(.))
    
    triplet_edge_edge = heterodimer_ilp_edges %>%
      transmute(i = ilp_edge_id + max(triplet_isotope$i), 
                j = ilp_edge_id + max(triplet_edge$j),
                v = 3)
    
    triplet_edge_node1 = heterodimer_ilp_edges %>%
      transmute(i = ilp_edge_id + max(triplet_isotope$i),
                j = ilp_nodes1,
                v = -1)
    
    triplet_edge_node2 = heterodimer_ilp_edges %>%
      transmute(i = ilp_edge_id + max(triplet_isotope$i),
                j = ilp_nodes2,
                v = -1)
    
    triplet_edge_node_link = heterodimer_ilp_edges %>%
      transmute(i = ilp_edge_id + max(triplet_isotope$i),
                j = ilp_nodes_link,
                v = -1)
    
    triplet_edge_heterodimer = bind_rows(triplet_edge_edge, 
                                         triplet_edge_node1,
                                         triplet_edge_node2,
                                         triplet_edge_node_link)
  }
    
  # triplet that force only one edge exist between two ilp_nodes
  # E.g. adduct vs fragment of NH3, HPO3 adduct and heterodimer, etc
  {
    double_edges = bind_rows(ilp_edges, heterodimer_ilp_edges %>% mutate(linktype = as.character(linktype))) %>%
      filter(category != "Library_MS2_fragment") %>% # allow Library_MS2_fragment to form double edge
      group_by(ilp_nodes1, ilp_nodes2) %>%
      mutate(n_twonodes = n()) %>%
      group_by(ilp_nodes1, ilp_nodes2, category) %>%
      mutate(n_category = n()) %>%
      filter(n_twonodes != n_category) %>% # Mainly to remove heterodimer when ilp_nodes1/2 are the same, but link different
      arrange(ilp_nodes1, ilp_nodes2) %>%
      group_by(ilp_nodes1, ilp_nodes2) %>%
      group_split() %>%
      bind_rows(.id = "double_edges_id") %>%
      mutate(double_edges_id = as.numeric(double_edges_id))
    
    triplet_double_edges_regular = double_edges %>%
      filter(category != "Heterodimer") %>%
      transmute(i = double_edges_id + max(triplet_edge_heterodimer$i), 
                j = ilp_edge_id + max(triplet_node$j),
                v = 1)
    
    triplet_double_edges_heterodimer = double_edges %>%
      filter(category == "Heterodimer") %>%
      transmute(i = double_edges_id + max(triplet_edge_heterodimer$i), 
                j = ilp_edge_id + max(triplet_edge$j),
                v = 1)
    
    triplet_edge_double_edges = bind_rows(triplet_double_edges_regular, 
                                          triplet_double_edges_heterodimer)
    
    nrow_triplet_double_edges = max(double_edges$double_edges_id)
  }
  
  # Generate sparse matrix on left hand side
  triplet_df = rbind(
    triplet_node,
    triplet_edge,
    triplet_isotope,
    triplet_edge_heterodimer,
    triplet_edge_double_edges
  )


  # converts the triplet into matrix
  mat = slam::simple_triplet_matrix(i=triplet_df$i,
                                    j=triplet_df$j,
                                    v=triplet_df$v)
  
  
  #CPLEX solver parameter
  {
    nc <- max(mat$j)
    obj <- c(ilp_nodes$cplex_score, 
             ilp_edges$cplex_score,
             heterodimer_ilp_edges$cplex_score)
    lb <- rep(0, nc)
    ub <- rep(1, nc)
    ctype <- rep("B",nc)
    
    nr <- max(mat$i)
    ## Three parts of constraints:
    ## 1. For each peak, binary sum of formula = 1. Each peak chooses unknown formula or 1 formula from potential formulas for the peak
    ## 2. For each edge, an edge exists only both formula it connects exists
    ## 3. For isotopic peak, an isotopic formula is given only if the isotope edge is chosen.
    ## 4. For heterodimer edge, an edge exists only when both formula and the linktype connection exist
    ## 5. For double edges, where two ilp_nodes are connected by more than 1 edge, force to choose at most 1 edge
    rhs = c(rep(1, max(ilp_rows)), 
            rep(0, nrow(ilp_edges)), 
            rep(0, nrow_triplet_isotope), 
            rep(0, nrow(heterodimer_ilp_edges)),
            rep(1, nrow_triplet_double_edges))
    sense <- c(rep("E", max(ilp_rows)), 
               rep("L", nrow(ilp_edges)),
               rep("E", nrow_triplet_isotope), 
               rep("L", nrow(heterodimer_ilp_edges)),
               rep("L", nrow_triplet_double_edges))
    
    triplet_df = triplet_df %>% arrange(j)
    cnt=as.vector(table(triplet_df$j))
    beg=vector()
    beg[1]=0
    for(i in 2:length(cnt)){beg[i]=beg[i-1]+cnt[i-1]}
    ind=triplet_df$i-1
    val = triplet_df$v
  }
  CPX_MAX = -1
  CPLEX_para = list(nc = nc,
                    nr = nr,
                    CPX_MAX = CPX_MAX,
                    obj = obj,
                    rhs = rhs,
                    sense = sense,
                    beg = beg,
                    cnt = cnt,
                    ind = ind, 
                    val = val,
                    lb = lb,
                    ub = ub,
                    ctype = ctype
  )
  
  return(CPLEX_para)
}
## Run_cplex ####
Run_cplex = function(CplexSet, obj_cplex, relative_gap = 1e-4, total_run_time = 3000){
  # obj_cplex = CplexSet$para$obj
    env <- openEnvCPLEX()
    prob <- initProbCPLEX(env)
    
    nc = CplexSet$para$nc
    nr = CplexSet$para$nr
    CPX_MAX = CplexSet$para$CPX_MAX
    rhs = CplexSet$para$rhs
    sense = CplexSet$para$sense
    beg = CplexSet$para$beg
    cnt = CplexSet$para$cnt
    ind = CplexSet$para$ind
    val = CplexSet$para$val
    lb = CplexSet$para$lb
    ub = CplexSet$para$ub
    ctype = CplexSet$para$ctype
    
    copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj = obj_cplex, rhs, sense,
                      beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
    
    copyColTypeCPLEX(env, prob, ctype)
    
    setDblParmCPLEX(env, CPX_PARAM_TILIM, total_run_time) # total run time
    
    # Conserve memory true
    setIntParmCPLEX(env, CPX_PARAM_MEMORYEMPHASIS, CPX_ON)
    setIntParmCPLEX(env, CPX_PARAM_PROBE, 3)
    # Sets which continuous optimizer will be used to solve the initial relaxation of a MIP.
    setIntParmCPLEX(env, 2025, 3) # 0 is default. one test shows 3 may works better and 4 is ok.
    
    
    # Sets an absolute tolerance on the gap between the best integer objective and the objective of the best node remaining
    setDblParmCPLEX(env, 2009, relative_gap) # Setting relative seems better
    
    
    # Sets an absolute tolerance on the gap between the best integer objective and the objective of the best node remaining
    # setDblParmCPLEX(env, 2008, 1000) 
    
    # setIntParmCPLEX(env, CPX_PARAM_INTSOLLIM, 2)
    # setDefaultParmCPLEX(env)
    # getChgParmCPLEX(env)
    
    
    # Assess parameters
    # getParmNameCPLEX(env, 2025)
    
    
    # Access Relative Objective Gap for a MIP Optimization Description
    # getMIPrelGapCPLEX(env, prob)
    
    
    tictoc::tic()
    # test = basicPresolveCPLEX(env, prob)
    return_code = mipoptCPLEX(env, prob)
    result_solution=solutionCPLEX(env, prob)
    # result_solution_info = solnInfoCPLEX(env, prob)
    
    print(paste(return_codeCPLEX(return_code),"-",
                status_codeCPLEX(env, getStatCPLEX(env, prob)),
                " - OBJ_value =", result_solution$objval))
    tictoc::toc()
    
    # writeProbCPLEX(env, prob, "prob.lp")
    delProbCPLEX(env, prob)
    closeEnvCPLEX(env)
  
  return(list(obj = obj_cplex, result_solution = result_solution))
}
## Test_para_CPLEX ####
Test_para_CPLEX = function(CplexSet, obj_cplex,  para = 0, para2 = 0, 
                           relative_gap=1e-4, total_run_time = 3000){
  
  # for(test_para_CPX_PARAM_PROBE in test_para1){
  for(rep_2 in 1:length(para2)){
    temp_para2 = para2[rep_2]
    # print(temp_para2)
  for(rep_1 in 1:length(para)){
    temp_para = para[rep_1]
    # print(test_para_CPX_PARAM_PROBE)
    print(temp_para)
    
    # obj_cplex = CplexSet$para$obj
    env <- openEnvCPLEX()
    prob <- initProbCPLEX(env)
    
    nc = CplexSet$para$nc
    nr = CplexSet$para$nr
    CPX_MAX = CplexSet$para$CPX_MAX
    rhs = CplexSet$para$rhs
    sense = CplexSet$para$sense
    beg = CplexSet$para$beg
    cnt = CplexSet$para$cnt
    ind = CplexSet$para$ind
    val = CplexSet$para$val
    lb = CplexSet$para$lb
    ub = CplexSet$para$ub
    ctype = CplexSet$para$ctype
    
    
    
    copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj = obj_cplex, rhs, sense,
                      beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
    
    
    copyColTypeCPLEX(env, prob, ctype)
    
    # Conserve memory true
    setIntParmCPLEX(env, CPX_PARAM_MEMORYEMPHASIS, CPX_ON)
    setIntParmCPLEX(env, CPX_PARAM_PROBE, 3)
    # MIP starting algorithm
    # Sets which continuous optimizer will be used to solve the initial relaxation of a MIP.
    setIntParmCPLEX(env, 2025, temp_para) # 3 and 4 seem better
    
    # Set time is Dbl not Int
    setDblParmCPLEX(env, CPX_PARAM_TILIM, total_run_time) # total run time
    # setIntParmCPLEX(env, CPX_PARAM_TUNINGTILIM, 200) # run time for each tuning (each optimizatoin run will test serveral tuning)
    
 
    
    # Sets an absolute tolerance on the gap between the best integer objective and the objective of the best node remaining
    setDblParmCPLEX(env, 2009, relative_gap) # Setting relative seems better
    
    
    
    # MIP subproblem algorithm
    # Decides which continuous optimizer will be used to solve the subproblems in a MIP, after the initial relaxation.
    # setIntParmCPLEX(env, 2026, temp_para) # Seems not much difference, one test shows default best
    
    # MIP variable selection strategy
    # setIntParmCPLEX(env, 2028, temp_para) #  Seems not much difference, one test shows 3 is best, but only 5%
    
    
    
    # setIntParmCPLEX(env, CPX_PARAM_BBINTERVAL, temp_para)
    # setIntParmCPLEX(env, CPX_PARAM_NODESEL, temp_para) # 0:3 No effect
    # setIntParmCPLEX(env, CPX_PARAM_CLIQUES, temp_para) # 0:3 No effect
    # 
    
    # setIntParmCPLEX(env, CPX_PARAM_CLIQUES, temp_para)
    
    # setIntParmCPLEX(env, CPX_PARAM_INTSOLLIM, 2)
    # setIntParmCPLEX(env, CPX_PARAM_PROBE, 2)
    # setDefaultParmCPLEX(env)
    # getChgParmCPLEX(env)
    
    # Assess parameters
    # getParmNameCPLEX(env, 1082)
    
    
    # Access Relative Objective Gap for a MIP Optimization Description
    # getMIPrelGapCPLEX(env, prob)
    
    tictoc::tic()
    # test = basicPresolveCPLEX(env, prob)
    return_code = mipoptCPLEX(env, prob)
    result_solution=solutionCPLEX(env, prob)
    # result_solution_info = solnInfoCPLEX(env, prob)
    
    print(paste(return_codeCPLEX(return_code),"-",
                status_codeCPLEX(env, getStatCPLEX(env, prob)),
                " - OBJ_value =", result_solution$objval))
    tictoc::toc()
    
    # writeProbCPLEX(env, prob, "prob.lp")
    delProbCPLEX(env, prob)
    closeEnvCPLEX(env)
    # }
  }
  }
  return(0)
}
## CPLEX_permutation ####
CPLEX_permutation = function(CplexSet, n_pmt = 5, sd_rel_max = 0.5){
  unknown_formula = CplexSet$formula$unknown_formula
  obj = CplexSet$Init_solution[[1]]$obj
  obj_node = obj[1:nrow(unknown_formula)]
  obj_edge = obj[(nrow(unknown_formula)+1):length(obj)]
  solution_ls = list()
  for(i in 1:n_pmt){
    
    temp_obj_edge = obj_edge + rnorm(length(obj_edge), mean = 0, sd = max(obj_edge) * sd_rel_max)
    temp_obj <- c(obj_node, temp_obj_edge)
    
    result_solution = Run_CPLEX(CplexSet, temp_obj)
    solution_ls[[length(solution_ls)+1]] = result_solution
    
  }
  return(solution_ls)
}




## Add_constraint_CPLEX ####
Add_constraint_CPLEX = function(CplexSet, obj){
  # obj = obj_cplex
  env <- openEnvCPLEX()
  prob <- initProbCPLEX(env)
  
  nc = CplexSet$para$nc
  nr = CplexSet$para$nr
  CPX_MAX = CplexSet$para$CPX_MAX
  rhs = CplexSet$para$rhs
  sense = CplexSet$para$sense
  beg = CplexSet$para$beg
  cnt = CplexSet$para$cnt
  ind = CplexSet$para$ind
  val = CplexSet$para$val
  lb = CplexSet$para$lb
  ub = CplexSet$para$ub
  ctype = CplexSet$para$ctype
  
  copyLpwNamesCPLEX(env, prob, nc, nr, CPX_MAX, obj = obj, rhs, sense,
                    beg, cnt, ind, val, lb, ub, NULL, NULL, NULL)
  copyColTypeCPLEX(env, prob, ctype)
  
  
  # nnz = sum((unknown_formula$ILP_result!=0)==T)
  # matbeg = 0
  # matval = rep(1,nnz)
  # matind = which(unknown_formula$ILP_result!=0)-1
  # addRowsCPLEX(env, prob, ncols=0, nrows=1, nnz=nnz, matbeg=matbeg, matind=matind, matval=matval,
  #              rhs = base::floor(nnz*.99), sense = "L",
  #              cnames = NULL, rnames = NULL)
  # addRowsCPLEX(env, prob, ncols=0, nrows=1, nnz=1, matbeg=0, matind=3909, matval=1,
  #              rhs = 1, sense = "E",
  #              cnames = NULL, rnames = NULL)
  # delRowsCPLEX(env, prob, begin = nr, end = getNumRowsCPLEX(env, prob)-1)
  # getNumRowsCPLEX(env, prob)
  
  # addMIPstartsCPLEX(env, prob, mcnt = 1, nzcnt = nc, beg = 0, varindices = 1:nc,
  #                   values = CplexSet$Init_solution2$CPLEX_x, effortlevel = 1, mipstartname = NULL)
  # 
  
  
  tictoc::tic()
  return_code = mipoptCPLEX(env, prob)
  result_solution=solutionCPLEX(env, prob)
  
  print(paste(return_codeCPLEX(return_code),"-",
              status_codeCPLEX(env, getStatCPLEX(env, prob)),
              " - OBJ_value =", result_solution$objval))
  
  tictoc::toc()
  
  # writeProbCPLEX(env, prob, "prob.lp")
  delProbCPLEX(env, prob)
  closeEnvCPLEX(env)
  return(list(obj = obj, result_solution = result_solution))
}
## Read CPLEX result ####
Read_cplex_result = function(solution){
  CPLEX_all_x=list()
  for(i in 1:length(solution)){
    CPLEX_all_x[[i]] =  solution[[i]]$result_solution$x
  }
  CPLEX_all_x = bind_cols(CPLEX_all_x)
  return(CPLEX_all_x)
}

# Annotation functions -------- ####
## core_annotate ####
core_annotate = function(ilp_nodes, FormulaSet_df, LibrarySet){
  # Make annotation to core peaks (where steps == 0)
  
  core_rank = FormulaSet_df %>%
    # filter(category == "Metabolite") %>%
    filter(steps == 0) %>%
    merge(LibrarySet %>%
            dplyr::select(library_id, name, note, origin, SMILES, status) %>% 
            dplyr::rename(parent_id = library_id)) %>%
    mutate(score_source = case_when(
      origin == "Manual_library" ~ 0.5,
      origin == "Library_known" ~ 0.5,
      status == "quantified" ~ 0.6,
      status == "detected" ~ 0.3,
      status == "expected" ~ 0,
      status == "predicted" ~ -0.5
    )) %>% 
    mutate(score_source = ifelse(transform == "", score_source, score_source-0.7)) %>%
    mutate(rank_score = score_source + prior_score - database_prior) %>%
    arrange(-rank_score) %>% 
    mutate(core_annotate = case_when(
      transform == "" ~ paste(name, parent_formula),
      direction == 1 ~ paste(name, parent_formula, "+", transform),
      direction == -1 ~ paste(name, parent_formula, "-", transform)
    )) %>%
    mutate(class = case_when(
      category == "Unknown" ~ "Unknown",
      category == "Metabolite" ~ "Metabolite", 
      T ~ "Artifact"
    )) %>%
    inner_join(ilp_nodes %>%
                 dplyr::select(node_id, formula, category, ilp_node_id))
  
  return(core_rank)
}
## core_annotate_unique ####
core_annotate_unique = function(core_annotation){
  core_annotation_filter = core_annotation %>%
    distinct(node_id, formula, category, .keep_all = T) %>%
    dplyr::select(node_id, rank_score)
  core_annotation_best = core_annotation %>%
    semi_join(core_annotation_filter)
  core_annotation_unique = core_annotation_best %>%
    distinct(node_id, formula, category, .keep_all = T)
}

## initiate_g_met ####
initiate_g_met = function(ilp_nodes, ilp_edges){
  ilp_nodes_met = ilp_nodes %>%
    filter(class == "Metabolite") %>%
    arrange(-ilp_result)
  ilp_edges_met = ilp_edges %>%
    filter(category == "Biotransform") %>%
    mutate(direction = ifelse(ilp_nodes1 < ilp_nodes2, 1, -1)) %>%
    dplyr::select(ilp_nodes1, ilp_nodes2, everything())
  
  g_met = graph_from_data_frame(ilp_edges_met, 
                                directed = T, 
                                vertices = ilp_nodes_met)
  
}

## initiate_g_nonmet ####
initiate_g_nonmet = function(ilp_nodes, ilp_edges, heterodimer_ilp_edges){
  
  ilp_edges_merge = merge(ilp_edges, heterodimer_ilp_edges, all = T) 
  
  ilp_edges_nonmet = ilp_edges_merge %>%
    filter(category != "Biotransform") %>%
    dplyr::select(ilp_nodes1, ilp_nodes2, everything())
  
  ilp_edges_nonmet_reorder1 = ilp_edges_nonmet %>%
    filter(direction == 1) %>%
    dplyr::rename(from = ilp_nodes1, to = ilp_nodes2)
  
  ilp_edges_nonmet_reorder2 = ilp_edges_nonmet %>%
    filter(direction == -1) %>%
    dplyr::rename(from = ilp_nodes2, to = ilp_nodes1)
  
  ilp_edges_nonmet_directed = bind_rows(ilp_edges_nonmet_reorder1, ilp_edges_nonmet_reorder2)
    
  ilp_edges_weights = ilp_edges_nonmet_directed %>%
    mutate(edge_weight = 2 - ilp_result - 0.1*cplex_score) %>%
    mutate(edge_weight = ifelse(category == "Heterodimer", edge_weight + 1.5, edge_weight)) %>%
    arrange(edge_weight)
  
  g_nonmet = graph_from_data_frame(ilp_edges_weights, 
                                directed = T, 
                                vertices = ilp_nodes)
  
  return(g_nonmet)
}

## initiate_met_dist_mat ####
initiate_met_dist_mat = function(g_met, ilp_nodes, core_annotation_unique){
  
  ilp_nodes_met = ilp_nodes %>%
    filter(class == "Metabolite") %>%
    arrange(-ilp_result)

  met_annotation_unique = core_annotation_unique %>%
    filter(category == "Metabolite") 
  # 2273 core_met * 15342 ilp_metabolites means 267 Mb ~ 35MB
  core_met_names = met_annotation_unique %>% pull(ilp_node_id) %>% as.character()
  
  # Keep order of from and to but treated as undirected graph
  distMatrix <- shortest.paths(g_met, 
                               v=core_met_names,
                               to=V(g_met),
                               mode = "all") 
  
  distMatrix[distMatrix==0] = Inf # So that shortest distance prevent finding itself
  
  return(distMatrix)
}

## initiate_nonmet_dist_mat ####
initiate_nonmet_dist_mat = function(g_nonmet, ilp_nodes, only_ilp_result = T){
  
  ilp_nodes_nonmet = ilp_nodes %>%
    filter(!category %in% c("Metabolite", "Unknown")) %>%
    arrange(-ilp_result)
  
  if(only_ilp_result){
    ilp_nodes_nonmet = ilp_nodes_nonmet %>%
      filter(ilp_result > 0.01)
  }

  target_names = ilp_nodes_nonmet %>% 
    pull(ilp_node_id) %>% as.character()
  
  core_names = ilp_nodes %>%
    filter(category != "Unknown") %>%
    filter(steps %% 1 == 0) %>% 
    pull(ilp_node_id) %>% as.character()
  
  g_edges = igraph::as_data_frame(g_nonmet, "edges") 
  
  # core_names = "61508"
  # target_names = "28"
  
  distMatrix <- shortest.paths(g_nonmet, 
                               v=core_names,
                               to=target_names, 
                               weights = g_edges$edge_weight,
                               mode = c("out"))
  distMatrix[distMatrix==0] = Inf 
  
  # An idea of compression 
  # Calculate sparse matrix on the run (Compression sparse column/row format)
  # Make sure each column can be extracted conviniently 
  
  return(distMatrix)
  
}

## track_annotation_met ####  
track_annotation_met = function(query_ilp_id, 
                                ilp_edges_annotate,
                                g_annotation = g_met, 
                                graph_path_mode = "all", 
                                dist_mat = met_dist_mat, 
                                core_annotation_unique){
  
  # g_annotation = g_met
  # dist_mat = met_dist_mat
  # query_ilp_id = 100000
  # graph_path_mode = "all"
  # ilp_edges_annotate = ilp_edges_annotate_met
  
  core_ilp_id = as.integer(row.names(dist_mat))
  
  if(!any(as.character(query_ilp_id) == colnames(dist_mat))){
    return("Not existed in distance matrix")
  }
  
  if(query_ilp_id %in% core_ilp_id){
    core_annotation_unique_filter = core_annotation_unique$ilp_node_id == query_ilp_id
    temp_annotation = core_annotation_unique$core_annotate[core_annotation_unique_filter]
    
    return(temp_annotation)
  }
  
  query_distMatrix = dist_mat[, as.character(query_ilp_id)]
  query_distMatrix_min = min(query_distMatrix)
  
  if(is.infinite(query_distMatrix_min)){
    return("No edge connections.")
  }
  
  {
    shortest_ilp_nodes = which(query_distMatrix == query_distMatrix_min)
    parent_selected = names(shortest_ilp_nodes)[1]
    paths_connect_ij_nodes = shortest_paths(g_annotation, 
                                            from = parent_selected, 
                                            to = as.character(query_ilp_id), 
                                            mode = graph_path_mode, 
                                            output = "vpath")
    
    
    ilp_node_path = paths_connect_ij_nodes[[1]] %>% unlist() %>% names() %>% as.numeric()
    
    temp_core = core_annotation_unique$core_annotate[core_annotation_unique$ilp_node_id == ilp_node_path[1]]
    
    transform_path = c()
    for(i in 1:(length(ilp_node_path)-1)){
      temp_ilp_edge_filter = ilp_edges_annotate$from == ilp_node_path[i] &
        ilp_edges_annotate$to == ilp_node_path[i+1]
      temp_ilp_edge = ilp_edges_annotate[temp_ilp_edge_filter,]
      
      if(dim(temp_ilp_edge)[1] == 0){
        temp_ilp_edge_filter = ilp_edges_annotate$to == ilp_node_path[i] &
          ilp_edges_annotate$from == ilp_node_path[i+1]
        temp_ilp_edge = ilp_edges_annotate[temp_ilp_edge_filter,]
        temp_ilp_edge$direction = -1
      }
      
      temp_sign = ifelse(temp_ilp_edge$direction != -1, "+", "-")
      temp_linktype = temp_ilp_edge$linktype
      temp_formula = ifelse(temp_ilp_edge$direction != -1, temp_ilp_edge$formula2, temp_ilp_edge$formula1)
      
      transform_path = paste(transform_path, temp_sign, temp_linktype, "->", temp_formula)
      
    }
    temp_annotation = paste0(temp_core, transform_path, collapse = "")
  }
  
  return(temp_annotation)
}


## track_annotation_nonmet ####  
track_annotation_nonmet = function(query_ilp_id, 
                                   ilp_edges_annotate,
                                g_annotation = g_nonmet, 
                                graph_path_mode = "out", 
                                dist_mat = nonmet_dist_mat, 
                                core_annotation_unique){
  
  # g_annotation = g_nonmet
  # graph_path_mode = "out"
  # dist_mat = nonmet_dist_mat
  # query_ilp_id = 287
  
  if(!any(as.character(query_ilp_id) == colnames(dist_mat))){
    return("Not existed in distance matrix")
  }
  
  core_ilp_id = as.integer(row.names(dist_mat))
  if(query_ilp_id %in% core_ilp_id){
    core_annotation_unique_filter = core_annotation_unique$ilp_node_id == query_ilp_id
    temp_annotation = core_annotation_unique$core_annotate[core_annotation_unique_filter]
    return(temp_annotation)
  } 
  
  query_distMatrix = dist_mat[, as.character(query_ilp_id)]
  query_distMatrix_min = min(query_distMatrix)
  
  if(is.infinite(query_distMatrix_min)){
    return("No edge connections.")
  }
  
  ## Network ####
  {
    shortest_ilp_nodes = which(query_distMatrix == query_distMatrix_min)
    parent_selected = names(shortest_ilp_nodes)[1]
    paths_connect_ij_nodes = shortest_paths(g_annotation, 
                                            from = parent_selected, 
                                            to = as.character(query_ilp_id), 
                                            mode = graph_path_mode, 
                                            output = "vpath")
    
    
    ilp_node_path = paths_connect_ij_nodes[[1]] %>% unlist() %>% names() %>% as.numeric()
    
    temp_core = core_annotation_unique$core_annotate[core_annotation_unique$ilp_node_id == ilp_node_path[1]]
    
    transform_path = c()
    for(i in 1:(length(ilp_node_path)-1)){
      temp_ilp_edge_filter = ilp_edges_annotate$from == ilp_node_path[i] &
        ilp_edges_annotate$to == ilp_node_path[i+1]
      temp_ilp_edge = ilp_edges_annotate[temp_ilp_edge_filter,]
      
      if(dim(temp_ilp_edge)[1] == 0){
        temp_ilp_edge_filter = ilp_edges_annotate$to == ilp_node_path[i] &
          ilp_edges_annotate$from == ilp_node_path[i+1]
        temp_ilp_edge = ilp_edges_annotate[temp_ilp_edge_filter,]
        temp_ilp_edge$direction = -1
      }
      
      temp_sign = case_when(
        temp_ilp_edge$category[1] == "Oligomer" ~ "*",
        temp_ilp_edge$category[1] == "Multicharge" ~ "/",
        temp_ilp_edge$category[1] == "Heterodimer" ~ "+ Peak",
        temp_ilp_edge$direction[1] == 1 ~ "+",
        temp_ilp_edge$direction[1] == -1 ~ "-"
      )
      
      temp_linktype = temp_ilp_edge$linktype[1]
      temp_formula = ifelse(temp_ilp_edge$direction[1] != -1, temp_ilp_edge$formula2, temp_ilp_edge$formula1)
      
      transform_path = paste(transform_path, temp_sign, temp_linktype, "->", temp_formula)
    }
    temp_annotation = paste0(temp_core, transform_path, collapse = "")
  }
  
  return(temp_annotation)
}

## Network annotation ####
empty_function = function(){
  {
    # Given a query id (node_id), find relevant ilp_id, distMatrix and closest nodes
    query_id = 996
    
    query_ilp_id = ilp_nodes_met %>% 
      filter(node_id==query_id) %>%
      pull(1)
    
    query_distMatrix = distMatrix[, query_ilp_id]
    
    if(is.null(dim(query_distMatrix))){
      shortest_ilp_nodes = list(which(query_distMatrix == min(query_distMatrix)))
    } else {
      shortest_ilp_nodes = apply(query_distMatrix, 2, function(x){
        which(x == min(x, na.rm=T))
      })
    }
  }
  
  # path
  # Given a pair of parent and query ilp_id, generate a network
  {
    parent_names = names(shortest_ilp_nodes[[1]])
    
    paths_connect_ij_nodes = lapply(parent_names, function(x){
      temp = shortest_paths(g_met, from = x, to = query_ilp_id[1], mode = "all")
      unlist(temp$vpath)
    })
    paths_rownum = unique(as.numeric(unlist(paths_connect_ij_nodes)))
    
    ilp_nodes_met_selected = ilp_nodes_met[paths_rownum,]
    
    ilp_id_selected = ilp_nodes_met$name[paths_rownum]
    ilp_edges_met_selected = ilp_edges_met %>%
      filter(from %in% ilp_id_selected, to %in% ilp_id_selected)
    
    g_met_selected = graph_from_data_frame(ilp_edges_met_selected, directed = F)
    
    plot.igraph(g_met_selected)
  }
}
# MS2 functions ---------- #####
## Read Xi MS2 format file ####
Add_MS2_nodeset = function(MS2_folder, NodeSet){
  
  if(!any(list.files()==MS2_folder)){
    return(NodeSet)
  }
  
  old_path = getwd()
  setwd(MS2_folder)
  
  MS2_filenames = list.files(pattern = ".xlsx")
  
  for(i_filename in 1:length(MS2_filenames)){
    MS2_filename = MS2_filenames[i_filename]
    
    ## Sensitive to format changes in MS2 input files
    WL_MS2 = read_xlsx(MS2_filename)
    IDs = WL_MS2$Comment %>%
      str_sub(start = 4) %>%
      as.numeric()
    
    sheetnames = excel_sheets(MS2_filename)[-c(1,2,3)]
    
    if(length(sheetnames)!=length(IDs)){
      warning(paste("Inconsistent number of peaks and spectra.", MS2_filename))
      next
    }

    MS2_ls = lapply(sheetnames, function(sheetname){
      df = read_xlsx(path=MS2_filename, sheet = sheetname, col_names = F, .name_repair = "minimal")
      if(ncol(df) == 1){
        return(NULL)
      }
      colnames(df) = c("mz", "inten")
      df = as.data.frame(df)
    })
    
    valid_MS2 = !sapply(MS2_ls, is.null)
    
    
    #########################
    
    id_map = Mset$Data$id
    names(id_map) = as.character(Mset$Data$Input_id)
    
    for(i in 1:length(IDs)){
      if(valid_MS2[i]){
        node_id = id_map[as.character(IDs[i])]
        if(is.na(node_id)){next}
        NodeSet[[node_id]]$MS2 = MS2_ls[[i]]
      }
    }
  }
  
  
  setwd(old_path)
  return(NodeSet)
}
## mergeMzIntensity ####
mergeMzIntensity = function(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3){
  
  
  # absTol = 1e-3
  # ppmTol = 10E-6
  # true_merge = merge(as.data.frame(spec1_df),spec2_df, by="mz", all=T)
  # Assuming mz came in sorted
  mz1 = spec1_df[,1]
  mz2 = spec2_df[,1]
  mzs = c(mz1, mz2)
  mzs = sort(mzs)
  inten1 = spec1_df[,2]
  inten2 = spec2_df[,2]
  keeps = c(1)
  count = 1
  MS_group = rep(1,(length(mzs)))
  for(i in 2:length(mzs)){
    # Determine when to treat two mz as same, and when they are different
    # Currently, when two masses are at least ppmTol different, and absTol different, they are considered different.
    if(mzs[i]-mzs[i-1]>mzs[i-1]*ppmTol & mzs[i]-mzs[i-1]>absTol){

      count = count+1
      keeps=c(keeps,i)
    }
    MS_group[i]=count
  }
  
  mz_result = mzs[keeps]
  count_mzs = 2
  count_keep = 1
  i = j = 1
  intens_mat = matrix(0, ncol=3, nrow = length(mz_result))
  while(count_mzs <= length(mzs)){
    if(MS_group[count_mzs] == MS_group[count_mzs-1]){
      intens_mat[count_keep,2] = inten1[i]
      intens_mat[count_keep,3] = inten2[j]
      i = i+1
      j = j+1
      count_mzs = count_mzs+2
    } else if(is.na(mz1[i])){
      intens_mat[count_keep,3] = inten2[j]
      count_mzs = count_mzs+1
      j=j+1
    } else if(mz1[i] == mzs[count_mzs-1]){
      intens_mat[count_keep,2] = inten1[i]
      i=i+1
      count_mzs = count_mzs+1
    } else {
      intens_mat[count_keep,3] = inten2[j]
      j=j+1
      count_mzs = count_mzs+1
    }
    count_keep = count_keep+1
  }
  # Handle boundary
  if(!is.na(mz1[i])){intens_mat[count_keep,2] = inten1[i]}
  if(!is.na(mz2[j])){intens_mat[count_keep,3] = inten2[j]}
  
  intens_mat[,1] = mz_result
  return(intens_mat)
}
mergeMzIntensity_backup = function(spec1_df, spec2_df, ppmTol = 10E-6, absTol = 1e-3){
  # if(identical(spec1_df, spec2_df))
  colnames(spec1_df)[2] = "inten1"
  colnames(spec2_df)[2] = "inten2"
  spec_df = merge(spec1_df, spec2_df, all=T, by = "mz")
  
  MS_group = rep(1,(nrow(spec_df)))
  mzs = spec_df$mz
  keeps = c(1)
  count = 1
  for(i in 2:length(mzs)){
    if(mzs[i]-mzs[i-1]>mzs[i-1]*ppmTol & mzs[i]-mzs[i-1]>absTol){
      count = count+1
      keeps=c(keeps,i)
    }
    MS_group[i]=count
  }
  spec_df[is.na(spec_df)]=0
  
  k_max=k_min=1
  while (k_max <= length(MS_group)){
    k_min = k_max
    while (MS_group[k_min] == MS_group[k_max]){
      k_max = k_max+1
      if(k_max > length(MS_group)){break}
    }
    # spec_df$mz[k_min]=max(spec_df$mz[k_min:(k_max-1)], na.rm = T)
    spec_df$inten1[k_min]=max(spec_df$inten1[k_min:(k_max-1)], na.rm = T)
    spec_df$inten2[k_min]=max(spec_df$inten2[k_min:(k_max-1)], na.rm = T)
  }
  spec_df = spec_df[keeps,]
  return(spec_df)
}
## Score_merge_MS2 ####
Score_merge_MS2 = function(spec_merge_df, mz_parent = 0){
  
  # A more realistic evaluation of two MS2 spectra simialrity is needed.
  # idea1 discount the weight of parent peak matched
  # idea2 add weights to the number of shared fragment
  # idea3 make sqrt of the intensity - 
  # which adds more weight to low-signal shared peaks
  
  # Implementing idea 3
  # a = a/max(a)
  # b = b/max(b)
  # a = sqrt(a)
  # b = sqrt(b)
  # Needs to use experimental data to validate the score cutoff for match and similarity
  # Implementing idea 1
  
  a = spec_merge_df[,2]
  b = spec_merge_df[,3]
  
  # 5e-4 and divide by 5 is arbitrary, needs experimental data to validate.
  mz_parent_position = which(abs(spec_merge_df[,1] - mz_parent) < 5e-4)
  if(length(mz_parent_position) != 0){
    a[mz_parent_position] = a[mz_parent_position]/5
    b[mz_parent_position] = b[mz_parent_position]/5
  }

  spec_score = a%*%b / sqrt(a%*%a * b%*%b)
  
  return(as.numeric(spec_score))
}



## ---------------------- ####  
## Deprecated ####
## read_library ####
# read_library = function(library_file){
#   data(isotopes)
#   hmdb_lib = read_csv(library_file)
#   hmdb_lib$formula = check_chemform(isotopes, hmdb_lib$formula)$new_formula
#   hmdb_lib$formula = sapply(hmdb_lib$formula, my_calculate_formula,"C1")
#   hmdb_lib$formula = sapply(hmdb_lib$formula, my_calculate_formula,"C1",-1)
#   hmdb_lib$mass = formula_mz(hmdb_lib$formula)
#   hmdb_lib["rdbe"]=formula_rdbe(hmdb_lib$formula)
#   return(hmdb_lib)
# }
## Special_fragment_connection ####
# Special_fragment_connection = function(NodeSet, FormulaSet, EdgeSet_MS2_fragment, MS2_library,
#                                        ppm_tol = 10, abs_tol = 1e-4){
#   # It requires the parent node has both msr MS2 and database MS2
#   # See comparison to General_fragment_connection
#   
#   parent_nodes = sapply(EdgeSet_MS2_fragment, "[[", "node2")
#   
#   FormulaSet_df = bind_rows(FormulaSet)
#   
#   MS2_library_external_id = sapply(MS2_library, "[[", 'external_id')
#   
#   FormulaSet_df_MS2 = FormulaSet_df %>%
#     filter(steps == 0) %>%
#     # filter(steps == 0, transform == "") %>%
#     merge(LibrarySet %>%
#             dplyr::select(library_id, note, mass), 
#           by.x = "parent_id",
#           by.y = "library_id") %>%
#     mutate(contain_msr_MS2 = node_id %in% parent_nodes) %>%
#     mutate(contain_lib_MS2 = note %in% MS2_library_external_id) %>%
#     filter(contain_msr_MS2, contain_lib_MS2)
#   
#   EdgeSet_MS2_fragment_df = bind_rows(EdgeSet_MS2_fragment) %>%
#     mutate(category = "Special_fragment")
#   
#   
#   H_mass = 1.00782503224
#   e_mass = 0.00054857990943
#   ion_mode = Mset$Global_parameter$mode
#   
#   Special_fragment = list()
#   for(i in 1:nrow(EdgeSet_MS2_fragment_df)){
#     id_node1 = EdgeSet_MS2_fragment_df$node1[i]
#     mz_node1 = NodeSet[[id_node1]]$mz + (H_mass-e_mass)*ion_mode
#     
#     id_node2 = EdgeSet_MS2_fragment_df$node2[i]
#     
#     FormulaSet_df_MS2_filter = FormulaSet_df_MS2 %>%
#       filter(node_id == id_node2)
#     if(nrow(FormulaSet_df_MS2_filter) == 0){next}
#     
#     for(j in 1:nrow(FormulaSet_df_MS2_filter)){
#       match_MS2_IDs = which(FormulaSet_df_MS2_filter$note[j] == MS2_library_external_id)
#       if(length(match_MS2_IDs) == 0){next}
#       match_MS2 = bind_rows(lapply(MS2_library[match_MS2_IDs], "[[","spectrum"))
#       
#       mz_dif = mz_node1 - match_MS2$mz
#       mz_tol = max(mz_node1 * ppm_tol * 1e-6, abs_tol)
#       mz_dif_position = which(abs(mz_dif) < mz_tol)
#       formula_node1 = unique(match_MS2$formula[mz_dif_position])
#       formula_node1 = formula_node1[!is.na(formula_node1)]
#       
#       length_formula_node1 = length(formula_node1)
#       if(length_formula_node1==0){next}
#       
#       formula_node2 = FormulaSet_df_MS2_filter$formula[j]
#       
#       Special_fragment[[length(Special_fragment)+1]] = list(node1 = rep(id_node1,length_formula_node1),
#                                                             node2 = rep(id_node2,length_formula_node1),
#                                                             formula1 = formula_node1,
#                                                             formula2 = rep(formula_node2,length_formula_node1))
#       
#     }
#   }
#   
#   special_fragment = bind_rows(Special_fragment) %>%
#     merge(EdgeSet_MS2_fragment_df) %>%
#     # group_by(node1, node2) %>%
#     # filter(n()>1) %>%
#     distinct()
#   
#   
#   
# }
## track_annotation ####
# track_annotation = function(ilp_nodes, ilp_edges){
#   
#   ilp_nodes_ls = list()
#   for(i in 1:nrow(ilp_nodes)){
#     # print(i)
#     temp_ilp_node_id = ilp_nodes_ilp$ilp_node_id[i]
#     
#     temp_ilp_edge = ilp_edges %>%
#       filter(ilp_nodes1 == temp_ilp_node_id | ilp_nodes2 == temp_ilp_node_id)
#     
#     temp_parent = temp_ilp_edge %>%
#       filter((ilp_nodes1 == temp_ilp_node_id & direction != 1 )| 
#                (ilp_nodes2 == temp_ilp_node_id & direction != -1)) %>%
#       arrange(-ilp_result, -cplex_score)
#     # temp_parent_high = temp_parent %>% filter(ilp_result > 0.01)
#     # temp_parent_low = temp_parent %>% filter(ilp_result <= 0.01)
#     
#     temp_child = temp_ilp_edge %>%
#       filter((ilp_nodes1 == temp_ilp_node_id & direction != -1 )| 
#                (ilp_nodes2 == temp_ilp_node_id & direction != 1)) %>%
#       arrange(-ilp_result)
#     # temp_child_high = temp_child %>% filter(ilp_result > 0.01)
#     # temp_child_low = temp_child %>% filter(ilp_result <= 0.01)
#     
#     ilp_nodes_ls[[as.character(temp_ilp_node_id)]] = list(parent = temp_parent,
#                                                           child = temp_child)
#   }
#   return(ilp_nodes_ls)
# }

## query_path #####
# query_path = function(query_node_id,
#                       ilp_nodes_ls, 
#                       LibrarySet){
# 
#   # query_node_id = 87
#   trace_ls = list()
#   ilp_id_ls = numeric()
#   
#   # For each ilp_nodes, find its best parent
#   
#   temp_parent = ilp_nodes_ls[[as.character(query_node_id)]]
#   
#     temp_ilp_node = ilp_nodes_ls[[87]] %>%
#       arrange(-ilp_result) %>%
#       slice(1)
#     
#     temp_ilp_edge = ilp_edges_ilp %>%
#       filter((ilp_nodes1 == temp_ilp_node$ilp_node_id & direction != 1 )| 
#              (ilp_nodes2 == temp_ilp_node$ilp_node_id & direction != -1)) %>%
#       slice(1)
#     
#   
#   
#   
#   
#   
#   node_id2 = temp_ilp_edge$node1
#   
#   temp_ilp_node2 = ilp_nodes_ls[[node_id2]] %>%
#     arrange(-ilp_result) %>%
#     filter(ilp_result > 0.01)
#   
#   temp_ilp_node = temp_ilp_node2$ilp_node_id
#   
#   if(temp_ilp_node2$steps == 0){
#     stop()
#   }
#   
#   temp_ilp_edge = ilp_edges_ilp %>%
#     filter((ilp_nodes1 == temp_ilp_node & direction != 1 )| 
#              (ilp_nodes2 == temp_ilp_node & direction != -1)) %>%
#     slice(1)
#   
#   
#   temp = result_ls[[87]]
#   current_step = max(temp$steps)
#     trace_ls = list()
#     id_ls = numeric()
#     while(nrow(temp) != 0){
#       result_ilp_node = temp %>% 
#         filter(ilp_result > 0.01) %>%
#         arrange(-ilp_result, -cplex_score, -score_prior_propagation, steps) %>%
#         slice(1)
#       
#       if(nrow(result_ilp_node)>1){
#         
#         
#       }
#       
#       if(nrow(result_ilp_node)==0){
#         result_ilp_node = temp %>% 
#           filter(steps <= current_step) %>%
#           arrange(-cplex_score, -score_prior_propagation, steps, category, parent_id) %>%
#           slice(1)
#       }
#       
#       parent_id = result_ilp_node$parent_id
#       if(parent_id %in% id_ls){
#         break
#       }
#       
#       if(result_ilp_node$steps == 0){
#         id_ls[[length(id_ls) + 1]] = parent_id
#         # trace_ls[[length(trace_ls) + 1]] = list(LibrarySet$name[LibrarySet$library_id == parent_id],
#         #                                         LibrarySet$formula[LibrarySet$library_id == parent_id])
#         if(result_ilp_node$transform == ""){
#           trace_ls[[length(trace_ls) + 1]] = c(LibrarySet$name[LibrarySet$library_id == parent_id],
#                                                LibrarySet$formula[LibrarySet$library_id == parent_id])
#         } else{
#           if(result_ilp_node$direction == 1){
#             temp_sign = "+"
#           } else {
#             temp_sign = "-"
#           }
#           trace_ls[[length(trace_ls) + 1]] = c(LibrarySet$name[LibrarySet$library_id == parent_id],
#                                                LibrarySet$formula[LibrarySet$library_id == parent_id],
#                                                temp_sign, result_ilp_node$transform)
#         }
#         
#         break
#       }
#       
#       if(result_ilp_node$category %in% c("Metabolite", "Adduct", "Fragment", 
#                                          "Ring_artifact", "Natural_abundance", "Radical",
#                                          "Library_MS2_fragment", "Experiment_MS2_fragment")){
#         temp_sign = ifelse(result_ilp_node$direction == 1, "+", "-")
#         transform = c(temp_sign, result_ilp_node$transform)
#       } else if(result_ilp_node$category %in% c("Heterodimer")){
#         transform = c("+ Peak", result_ilp_node$transform)
#       } else if(result_ilp_node$category %in% c("Multicharge")){
#         transform = c("/", result_ilp_node$transform)
#       } else if(result_ilp_node$category %in% c("Oligomer")){
#         transform = c("*", result_ilp_node$transform)
#       }
#       
#       trace_ls[[length(trace_ls) + 1]] = c(transform, "->", result_ilp_node$formula)
#       id_ls[[length(id_ls) + 1]] = parent_id
#       
#       current_step = result_ilp_node$steps
#       
#       temp = result_ls[[parent_id]] %>%
#         filter(formula == result_ilp_node$parent_formula) 
#     }
#   }
#   
#   
#   # Output formating
#   {
#     paste_combine = c()
#     for(i in length(trace_ls):1){
#       paste_combine = c(paste_combine, trace_ls[[i]])
#     }
#   }
#   
#   paste(paste_combine, collapse = " ")
# }
