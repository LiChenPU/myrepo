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
  # setwd(dirname(rstudioapi::getSourceEditorContext()$path))
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
    mutate(mass = formula_mz(formula),
           rdbe = formula_rdbe(formula))
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
    rename(Input_id = id) %>%
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
      formula = list()
    )
  })
  names(NodeSet) = 1:nrow(Mset$Data)
  return(NodeSet)
}
## Initiate_edgeset ####
Initiate_edgeset = function(Mset, NodeSet, mz_tol_abs = 0, mz_tol_ppm = 10, rt_tol_bio = Inf, rt_tol_nonbio = 0.2){
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
        mass_tol = max(temp_mz_list[i]*mz_tol_ppm,mz_tol_abs)
        while(1){
          j=j+1
          if(j>merge_nrow){break}
          temp_ms = temp_mz_list[j]-temp_mz_list[i]
          
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
    mutate(category = "Metabolite") %>%
    dplyr::rename(note = accession) %>%
    mutate(origin = "Library_HMDB")
  
  Metabolites_known = Mset$Library_known %>%
    mutate(category = "Metabolite") %>%
    dplyr::rename(note = HMDB) %>%
    mutate(origin = "Library_known") %>%
    mutate(rt = .[,eval(Mset$Global_parameter$LC_method)])

  Adducts = Mset$Empirical_rules %>%
    filter(category == "Adduct") %>%
    mutate(category = "Artifact") %>%
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
      node1 = as.vector(x["node1"]),
      node2 = as.vector(x["node2"]),
      category = "Ring_artifact",
      linktype = "Ring_artifact",
      direction = 1
    )
  })
  
  return(EdgeSet_ring_artifact)
}

## Oligomer_connection ####
Oligomer_connection = function(peak_group, ppm_tol = 10){
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
      # linktype = as.vector(paste0("x", x["ratio"])),
      direction = 0
    )
  })
  
  return(EdgeSet_oligomer)
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
    mutate(inten1 = node_inten[node1]) %>%
    filter(inten1 > log10(inten_threshold)) %>% # retain only high intensity as node1
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
      node1 = as.vector(x["node1"]),
      node2 = as.vector(x["node2"]),
      category = "Heterodimer",
      linktype = as.character(x["linktype"]),
      direction = 1
    )
  })
  
  return(EdgeSet_heterodimer)
}
 
## merge_edgeset ####
merge_edgeset = function(EdgeSet, ...){
  
  EdgeSet_all_df = bind_rows(EdgeSet, ...) %>%
    mutate(edge_id = 1:nrow(.))
  return(EdgeSet_all_df)
}
### Expand_library ####
expand_library = function(lib, rule, direction, category){
  
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
### Match_library ####
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
    if(record_RT_tol < 99){
      candidate = candidate %>%
        filter(abs(node_RT[i] - node_RT[as.character(parent_id)]) < record_RT_tol) 
      
      if(nrow(candidate) == 0){
        i = i + 1
        next
      }
      # We may implement another filter to break a cycle a->b->a, but it is not necessary
      # Actually, keeping both propagation direction makes it easier to trace parents
    }
    
    adding = candidate %>%
      mutate(node_id = temp_id[i]) %>%
      mutate(steps = current_step)
    sf[[temp_id[i]]] = bind_rows(sf[[temp_id[i]]], adding)
    
    i = i+1
  }
  return(sf)
}

## Initialize formula pool ####
Initilize_empty_formulaset = function(NodeSet){
  sf_str = data.frame(
    formula = as.character(),
    mass = as.numeric(), 
    rdbe = as.numeric(), 
    category = as.character(), 
    parent_id = as.numeric(), 
    parent_formula = as.character(),
    transform = as.character(),
    direction = as.integer(),
    stringsAsFactors = F
  )
  sf = lapply(1:length(NodeSet), function(x) sf_str)
  return(sf)
}

## Match_library_formulaset ####
Match_library_formulaset = function(FormulaSet, Mset, NodeSet, LibrarySet, 
                                    expand = F,
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
    filter(category == "Artifact") %>%
    filter(T)
  lib_manual = seed_library %>%
    filter(category == "Manual")
  ## Expansion of starting formula/structures ####
  if(expand){
    initial_rule = Mset$Empirical_rules %>%
      filter(category %in% c("Biotransform", "Adduct"))

    rule_1 = initial_rule %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,1))
    rule_2 = initial_rule %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,-1))

    lib_met_1 = expand_library(lib_met, rule_1, direction = 1, category = "Metabolite")
    lib_met_2 = expand_library(lib_met, rule_2, direction = -1, category = "Metabolite")
    lib_met = bind_rows(lib_met_1, lib_met_2, lib_met)

    rule_1 = initial_rule %>% filter(category == "Adduct") %>% filter(direction %in% c(0,1))
    rule_2 = initial_rule %>% filter(category == "Adduct") %>% filter(direction %in% c(0,-1))

    lib_adduct_1 = expand_library(lib_adduct, rule_1, direction = 1, category = "Artifact")
    lib_adduct_2 = expand_library(lib_adduct, rule_2, direction = -1, category = "Artifact")

    lib_adduct = bind_rows(lib_adduct, lib_adduct_1, lib_adduct_2)
  }

  initial_lib = bind_rows(lib_met, lib_adduct, lib_manual) %>%
    # mutate(RT = -1) %>%
    # filter(grepl("(?<!H)-.", formula, perl=T))
    filter(!grepl("-|NA", formula)) %>% # in case a Rb1H-1 is measured, it will get filtered
    arrange(mass)
  
  sf = FormulaSet
  sf = match_library(lib = initial_lib, sf, 
                     record_ppm_tol = ppm_tol, 
                     record_RT_tol = Inf, 
                     current_step = 0, 
                     NodeSet)
  
  return(sf)

}


## propagate_ring_artifact ####
propagate_ring_artifact = function(new_nodes_df, sf, EdgeSet_ring_artifact, NodeSet, current_step){
  node_mass = sapply(NodeSet, "[[", "mz")
  
  node1 = sapply(EdgeSet_ring_artifact, "[[", "node1")
  node2 = sapply(EdgeSet_ring_artifact, "[[", "node2")
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
## propagate_oligomer ####
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
  
  
  new_nodes_df_oligomer = bind_rows(new_nodes_df1, new_nodes_df2)
  

  
  
  return(new_nodes_df_oligomer)
}
## propagate_heterodimer ####
propagate_heterodimer = function(new_nodes_df, sf, EdgeSet_heterodimer, NodeSet, current_step, propagation_ppm_threshold){
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
  
  heterodimer_ls = list()
  
  if(nrow(new_nodes_df_heterodimer) == 0){
    return(NULL)
  }
 
  
  for(i in 1:nrow(new_nodes_df_heterodimer)){
    temp = new_nodes_df_heterodimer[i,]
    transform = sf[[as.numeric(temp$transform)]] %>%
      distinct(formula, .keep_all = T) %>%
      filter(node_mass[node_id] - mass < mass * propagation_ppm_threshold,
             !grepl("\\.|Ring_artifact", formula)) %>%
      # Need to avoid exponential growth because heterodimer entries grows at n^2 rate
      filter(category != "Heterodimer") %>%
      filter(steps < 0.02 | (steps >= 1  & steps <1.02))
    
    if(nrow(transform) == 0){next}
    heterodimer_ls[[length(heterodimer_ls)+1]] = list(parent_id = rep(temp$parent_id, length(transform$formula)),
                                                      transform = rep(temp$transform, length(transform$formula)),
                                                      node2 = rep(temp$node2, length(transform$formula)),
                                                      formula = my_calculate_formula(temp$formula, transform$formula),
                                                      rdbe = temp$rdbe + transform$rdbe,
                                                      mass = temp$mass + transform$mass)
    
  }
  
  heterodimer = bind_rows(heterodimer_ls) %>%
    merge(new_nodes_df_heterodimer %>% dplyr::select(-c("formula", "mass", "rdbe")), all.x = T) %>%
    mutate(node_id = node2) %>%
    dplyr::select(-node2)
  
  # for(i in unique(new_nodes_df_heterodimer$node_id)){
  #   sf[[i]] = bind_rows(sf[[i]], new_nodes_df_heterodimer[new_nodes_df_heterodimer$node_id == i, ])
  # }
  
  return(heterodimer)
}
## Propagate_formulaset ####
Propagate_formulaset = function(Mset, 
                                NodeSet,
                                FormulaSet,
                                biotransform_step = 5,
                                artifact_step = 5,
                                propagation_ppm_threshold = 1e-6,
                                record_RT_tol = 0.1,
                                record_ppm_tol = 5e-6)
{

  # biotransform_step = 3
  # artifact_step = 4
  # propagation_ppm_threshold = 1e-6
  # record_RT_tol = 0.1
  # record_ppm_tol = 5e-6

  sf = FormulaSet
  empirical_rules = Mset$Empirical_rules
  node_mass = sapply(NodeSet, "[[", "mz") %>% sort()
  node_RT = sapply(NodeSet, "[[", "RT")[names(node_mass)]
  temp_id = as.numeric(names(node_mass)) # numeric

  
  ## Expansion 
  timer = Sys.time()
  step_count = 0
  while(step_count < biotransform_step){
    all_nodes_df = bind_rows(sf)
    # Handle artifacts
    sub_step = 0
    while(sub_step < 0.01 * artifact_step){
      all_nodes_df = bind_rows(sf)
      new_nodes_df = all_nodes_df %>%
        distinct(node_id,formula, .keep_all=T) %>% # This garantee only new formulas will go to next propagation
        filter(steps==(step_count + sub_step)) %>% # This garantee only formulas generated from the step go to next propagation
        filter(rdbe > -1) %>% # filter out formula has ring and double bind less than -1
        filter(!grepl("\\.|Ring_artifact", formula)) %>%  # filter formula with decimal point and ring artifacts
        filter(node_mass[as.character(node_id)] - mass < propagation_ppm_threshold * mass) # propagate from accurate formulas
      
      
      if(nrow(new_nodes_df)==0){break}
      
      sub_step = sub_step+0.01
      current_step = step_count + sub_step
      print(paste("Step",step_count + sub_step,"elapsed="))
      print((Sys.time()-timer))
      
      
      rule_1 = empirical_rules %>% filter(category != "Biotransform") %>% filter(direction %in% c(0,1))
      rule_2 = empirical_rules %>% filter(category != "Biotransform") %>% filter(direction %in% c(0,-1))
      
      lib_1 = expand_library(new_nodes_df, rule_1, direction = 1, category = "Artifact")
      lib_2 = expand_library(new_nodes_df, rule_2, direction = -1, category = "Artifact")
      
      lib_adduct = bind_rows(lib_1, lib_2) %>%
        filter(!grepl("-|NA", formula)) %>% # in case a Rb1H-1 is measured
        # mutate(RT = node_RT[as.character(parent_id)]) %>%
        arrange(mass)
      

      ## sf[[11]] to test if [13]C1C15H30O2 is added
      ring_artifact = propagate_ring_artifact(new_nodes_df, sf, EdgeSet_ring_artifact, NodeSet, current_step)
      ## sf[[302]]
      oligomer = propagate_oligomer(new_nodes_df, sf, EdgeSet_oligomer, NodeSet, current_step)
      ## sf[[460]]
      heterodimer = propagate_heterodimer(new_nodes_df, sf, EdgeSet_heterodimer, NodeSet, current_step, propagation_ppm_threshold)
      ## sf[[436]]
      sf_add = bind_rows(ring_artifact, oligomer, heterodimer)
      
      sf = match_library(lib_adduct,
                         sf,
                         record_ppm_tol,
                         record_RT_tol,
                         current_step,
                         NodeSet)
      
      for(i in unique(sf_add$node_id)){
        sf[[i]] = bind_rows(sf[[i]], sf_add[sf_add$node_id == i, ])
      }

    }
    all_nodes_df = bind_rows(sf)
    new_nodes_df = all_nodes_df %>%
      filter(category == "Metabolite") %>% # only metabolites go to biotransformation, also garantee it is not filtered.
      distinct(node_id,formula, .keep_all=T) %>% # This garantee only new formulas will go to next propagation
      filter(steps == step_count) %>% # only formulas generated from the step go to next propagation
      filter(rdbe > -1) %>% # filter out formula has ring and double bind less than -1
      filter(!grepl("\\.|Ring_artifact", formula)) %>%  # filter formula with decimal point and ring artifacts
      filter(node_mass[as.character(node_id)] - mass < propagation_ppm_threshold * mass) # propagate from accurate formulas
    
    step_count = step_count+1
    current_step = step_count
    print(paste("Step",current_step,"elapsed="))
    print(Sys.time()-timer)
    
    if(nrow(new_nodes_df)==0){break}
    
    rule_1 = empirical_rules %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,1))
    rule_2 = empirical_rules %>% filter(category == "Biotransform") %>% filter(direction %in% c(0,-1))
    
    new_nodes_df = new_nodes_df %>%
      mutate(parent_id = node_id)
    lib_1 = expand_library(new_nodes_df, rule_1, direction = 1, category = "Metabolite")
    lib_2 = expand_library(new_nodes_df, rule_2, direction = -1, category = "Metabolite")
    
    lib_met = bind_rows(lib_1, lib_2) %>%
      filter(!grepl("-|NA", formula)) %>% # in case a Rb1H-1 is measured
      # mutate(RT = node_RT[as.character(parent_id)]) %>%
      arrange(mass)
    
    sf = match_library(lib_met,
                       sf,
                       record_ppm_tol,
                       record_RT_tol=Inf,
                       current_step,
                       NodeSet)

  }
  
  return(sf)
}

## Score_formulaset ####
Score_formulaset = function(FormulaSet,
                            database_match = 0.5, 
                            rt_match = 1, 
                            known_rt_tol = 0.5,
                            manual_match = 1,
                            bio_decay = 1,
                            artifact_decay = -0.5){
  

  # database_match = 0.5
  # rt_match = 1
  # known_rt_tol = 0.5
  # manual_match = 1
  # bio_decay = -0.2
  # artifact_decay = -0.1
  
  FormulaSet_df = bind_rows(FormulaSet)
  
  # Score HMDB and known adduct match
  {
    FormulaSet_df = FormulaSet_df %>%
      mutate(database_prior = case_when(
        category == "Manual" ~ manual_match,
        steps == 0 & transform == "" ~ database_match,
        TRUE ~ 0 # Everything else
      ))
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
    # Slow here - may need optimization 
    # Consider use stringr to substring
    formula_P = sapply(temp_formula, elem_num_query, "P")
    formula_Si = sapply(temp_formula, elem_num_query, "Si")
    formula_O = sapply(temp_formula, elem_num_query, "O")
    FormulaSet_df = FormulaSet_df %>%
      mutate(empirical_POratio_prior = ifelse(formula_O >= 3*formula_P & formula_O >= 2*formula_Si, 0, -10))
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
  
  # Score propagation
  {
    summary_ls = list()
    
    initial = FormulaSet_df %>%
      filter(steps == 0) %>%
      mutate(score = prior_score)
    
    summary_ls[[length(summary_ls)+1]] = initial
    
    nonzero = initial %>%
      filter(score > 0) %>%
      dplyr::select(node_id, formula, score) %>%
      distinct(node_id, formula, .keep_all=T)
    
    unique_steps = unique(FormulaSet_df)$steps
    step = 0
    floating_error = 1e-8
    while(any(abs(step - unique_steps) < floating_error)){
      # artifact step
      sub_step = step + 0.01
      sub_nonzero = nonzero 
      while(any(abs(sub_step - unique_steps) < floating_error)){
        temp = FormulaSet_df %>%
          filter(abs(steps - sub_step) < floating_error)
        
        temp = temp %>%
          left_join(sub_nonzero, by = c("parent_id"="node_id", "parent_formula" = "formula")) %>%
          mutate(score = score + artifact_decay + prior_score)
        
        sub_nonzero = temp %>%
          filter(score > 0) %>%
          dplyr::select(node_id, formula, score) %>%
          distinct(node_id, formula, .keep_all=T)
        
        summary_ls[[length(summary_ls)+1]] = temp
        
        sub_step = sub_step + 0.01
      }
      
      step = step + 1
      temp = FormulaSet_df %>%
        filter(abs(steps - step) < floating_error)
      
      temp = temp %>%
        left_join(nonzero, by = c("parent_id"="node_id", "parent_formula" = "formula")) %>%
        mutate(score = score + bio_decay + prior_score)
      
      summary_ls[[length(summary_ls)+1]] = temp
      
      nonzero = temp %>% 
        filter(score > 0) %>%
        dplyr::select(node_id, formula, score) %>%
        distinct(node_id, formula, .keep_all=T)
      
    }
    
    summary = bind_rows(summary_ls) %>%
      mutate(score = ifelse(is.na(score), prior_score, score)) %>%
      mutate(score = ifelse(score < 0 & score > -5, 0, score)) %>%
      dplyr::rename(score_prior_propagation = score)
  }
  
  return(summary)
}


## initiate_ilp_nodes ####
initiate_ilp_nodes = function(FormulaSet_df){
  ilp_nodes = FormulaSet_df %>%
    arrange(-score_prior_propagation) %>% 
    mutate(temp_cat = ifelse(category == "Metabolite", "Metabolite", "Artifact")) %>%
    distinct(node_id, formula, temp_cat, .keep_all=T) %>%
    dplyr::select(-temp_cat) %>%
    arrange(node_id) %>%
    mutate(ilp_node_id = 1:nrow(.)) %>%
    dplyr::select(ilp_node_id, everything()) %>%
    filter(T)
  return(ilp_nodes)
  
}

## score_ilp_nodes ####
score_ilp_nodes = function(ilp_nodes, MassDistsigma = MassDistsigma, 
                           formula_score = 1){
  
  node_mass = sapply(NodeSet, "[[", "mz")
  
  ilp_nodes = ilp_nodes %>%
    mutate(msr_mass = node_mass[node_id],
           ppm_error = (mass - msr_mass) / mass * 1e6) %>%
    mutate(score_formula = formula_score)
  
  # Score mass accuracy
  ilp_nodes = ilp_nodes %>%
    mutate(score_mass = dnorm(ppm_error, 0, MassDistsigma)/dnorm(0, 0, MassDistsigma)) %>%
    mutate(score_mass = score_mass+1e-10) %>%
    mutate(score_mass = log10(score_mass)) %>%
    filter(T)
  
  cplex_score_node = ilp_nodes %>% 
    dplyr::select(starts_with("score_")) %>% 
    rowSums()
  
  ilp_nodes = ilp_nodes %>%
    mutate(cplex_score = cplex_score_node)
  
  return(ilp_nodes)
}

## initiate_ilp_edges ####
initiate_ilp_edges = function(EdgeSet_all_df, CplexSet){
  ilp_nodes = CplexSet$ilp_nodes %>%
    arrange(ilp_node_id) 
  
  EdgeSet_df = EdgeSet_all_df %>%
    filter(category != "Heterodimer")

  node_id_ilp_node_mapping = split(ilp_nodes$ilp_node_id, ilp_nodes$node_id)
  
  ilp_nodes_met = ilp_nodes %>% filter(category == "Metabolite") 
  node_id_ilp_node_mapping_met = split(ilp_nodes_met$ilp_node_id, ilp_nodes_met$node_id)
  
  ilp_nodes_nonmet = ilp_nodes %>% filter(category != "Metabolite") 
  node_id_ilp_node_mapping_nonmet = split(ilp_nodes_nonmet$ilp_node_id, ilp_nodes_nonmet$node_id)
  
  ilp_nodes_formula = ilp_nodes %>%
    pull(formula)
  
  match_matrix_index_ls = list()
  
  for(i in 1:nrow(EdgeSet_df)){
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
    } else if(category == "Fragment") {
      # node1 must be artifact, and node2 can be anything
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
    } else if(category == "Oligomer") {
      # Both can be anything for simplicity
      node1 = EdgeSet_df$node1[i]
      ilp_nodes1 = node_id_ilp_node_mapping[[as.character(node1)]]
      formula1 = ilp_nodes_formula[ilp_nodes1]
      if(length(formula1) == 0){next}
      node2 = EdgeSet_df$node2[i]
      ilp_nodes2 = node_id_ilp_node_mapping[[as.character(node2)]]
      formula2 = ilp_nodes_formula[ilp_nodes2]
      if(length(formula2) == 0){next}
      
      fold = as.numeric(EdgeSet_df$linktype[i])
      formula_transform = mapply(my_calculate_formula, formula1, formula1, fold-1)
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
      
      if(any(category == c("Adduct", "Natural_abundance"))){
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

  node_id_ilp_node_mapping = split(ilp_nodes$ilp_node_id, ilp_nodes$node_id)
  
  ilp_nodes_met = ilp_nodes %>% filter(category == "Metabolite") 
  node_id_ilp_node_mapping_met = split(ilp_nodes_met$ilp_node_id, ilp_nodes_met$node_id)
  
  ilp_nodes_nonmet = ilp_nodes %>% filter(category != "Metabolite") 
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
score_ilp_edges = function(ilp_edges, NodeSet, MassDistsigma = MassDistsigma, 
                           rule_score_biotransform = 0.1, rule_score_artifact = 1, 
                           rule_score_oligomer = 1, rule_score_ring_artifact = 5,
                           inten_score_isotope = 1){

  # rule_score_biotransform = 0
  # rule_score_artifact = 1
  # rule_score_oligomer = 1
  # rule_score_ring_artifact = 5
  # inten_score_isotope = 1
    
  # Score rule category 
  ilp_edges = CplexSet$ilp_edges %>%
    mutate(score_category = case_when(
      category == "Biotransform" ~ rule_score_biotransform,
      category == "Ring_artifact" ~ rule_score_ring_artifact,
      category == "Oligomer" ~ rule_score_oligomer,
      category != "Biotransform" ~ rule_score_artifact, 
    )) %>%
    filter(T)
  
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
      mutate(score_inten_isotope = log10(p_obs/p_theory+1e-10) + inten_score_isotope)
    
    ilp_edges = ilp_edges %>%
      merge(ilp_edges_isotope %>% dplyr::select(ilp_edge_id, score_inten_isotope), all.x = T)
    # hist(ilp_edges_isotope %>% filter(score_inten_isotope<1 & score_inten_isotope>-8 ) %>% pull(score_inten_isotope))
  }
  
  cplex_score_edge = ilp_edges %>% 
    dplyr::select(starts_with("score_")) %>% 
    rowSums(na.rm=T)
  
  ilp_edges = ilp_edges %>%
    mutate(cplex_score = cplex_score_edge)
  
  return(ilp_edges)
}

## score_heterodimer_ilp_edges ####
score_heterodimer_ilp_edges = function(CplexSet, rule_score_heterodimer = 1){
  
  heterodimer_ilp_edges = CplexSet$heterodimer_ilp_edges %>%
    mutate(score_category = case_when(
      category == "Heterodimer" ~ rule_score_heterodimer,
      category != "Heterodimer" ~ 0
    )) %>%
    filter(T)
  
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
  # 
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
    
  # Generate sparse matrix on left hand side
  triplet_df = rbind(
    triplet_node,
    triplet_edge,
    triplet_isotope,
    triplet_edge_heterodimer
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
    ## 1. For each peak, binary sum of formula <= 1. Each peak chooses 0 or 1 formula from potential formulas for the peak
    ## 2. For each edge, an edge exists only both formula it connects exists
    ## 3. For isotopic peak, an isotopic formula is given only if the isotope edge is chosen.
    ## 4. For heterodimer edge, an edge exists only when both formula and the linktype connection exist
    rhs = c(rep(1, max(ilp_rows)), rep(0, nrow(ilp_edges)), rep(0, nrow_triplet_isotope), rep(0, nrow(heterodimer_ilp_edges)))
    sense <- c(rep("L", max(ilp_rows)), rep("L", nrow(ilp_edges)), rep("E", nrow_triplet_isotope), rep("L", nrow(heterodimer_ilp_edges)))
    
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
Run_cplex = function(CplexSet, obj_cplex){
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
  # setIntParmCPLEX(env, CPX_PARAM_INTSOLLIM, 2)
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
  
  return(list(obj = obj_cplex, result_solution = result_solution))
}
## Test_para_CPLEX ####
Test_para_CPLEX = function(CplexSet, obj_cplex,  test_para = -1:4){
  
  # for(test_para_CPX_PARAM_PROBE in test_para1){
  for(temp_para in test_para){
    
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
    
    # Set time is Dbl not Int
    # setDblParmCPLEX(env, CPX_PARAM_TILIM, 1000) # total run time
    # setIntParmCPLEX(env, CPX_PARAM_TUNINGTILIM, 200) # run time for each tuning (each optimizatoin run will test serveral tuning)
    
    
    
    # setIntParmCPLEX(env, CPX_PARAM_BBINTERVAL, temp_para)
    # setIntParmCPLEX(env, CPX_PARAM_NODESEL, temp_para) # 0:3 No effect
    setIntParmCPLEX(env, CPX_PARAM_CLIQUES, temp_para) # 0:3 No effect
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
## query_path ####
query_path = function(query_node_id = 6,
                      result_ls, 
                      LibrarySet){
  
  # query_node_id = 204
  # test1 = result_ls[[503]]
  # test2 = result_ls[[204]]
  temp = result_ls[[query_node_id]]
  current_step = max(temp$steps)
  if(is.null(temp)){
    trace_ls = NULL
    id_ls = NULL
  } else {
    trace_ls = list()
    id_ls = numeric()
    while(nrow(temp) != 0){
      result_ilp_node = temp %>% 
        filter(ilp_result > 0.01) %>%
        arrange(-ilp_result, -cplex_score, -score_prior_propagation, steps, parent_id) %>%
        slice(1)
      
      if(nrow(result_ilp_node)==0){
        result_ilp_node = temp %>% 
          filter(steps <= current_step) %>%
          arrange(-cplex_score, -score_prior_propagation, steps, parent_id) %>%
          slice(1)
      }
      
      parent_id = result_ilp_node$parent_id
      if(parent_id %in% id_ls){
        break
      }
      
      if(result_ilp_node$steps == 0){
        id_ls[[length(id_ls) + 1]] = parent_id
        # trace_ls[[length(trace_ls) + 1]] = list(LibrarySet$name[LibrarySet$library_id == parent_id],
        #                                         LibrarySet$formula[LibrarySet$library_id == parent_id])
        trace_ls[[length(trace_ls) + 1]] = c(LibrarySet$name[LibrarySet$library_id == parent_id],
                                             LibrarySet$formula[LibrarySet$library_id == parent_id])
        break
      }
      
      if(result_ilp_node$category %in% c("Metabolite", "Manual", "Artifact", "Ring_artifact")){
        if(result_ilp_node$direction == 1){
          temp_sign = "+"
        } else {
          temp_sign = "-"
        }
        transform = c(temp_sign, result_ilp_node$transform)
      } else if(result_ilp_node$category %in% c("Heterodimer")){
        transform = c("+ Peak", result_ilp_node$transform)
      } else if(result_ilp_node$category %in% c("Multicharge")){
        transform = c("/", result_ilp_node$transform)
      } else if(result_ilp_node$category %in% c("Oligomer")){
        transform = c("*", result_ilp_node$transform)
      }
      
      trace_ls[[length(trace_ls) + 1]] = c(transform, "->", result_ilp_node$formula)
      id_ls[[length(id_ls) + 1]] = parent_id
      
      current_step = result_ilp_node$steps
      
      temp = result_ls[[parent_id]] %>%
        filter(formula == result_ilp_node$parent_formula) 
      
    }
    
  }
  
  
  # Output formating
  {
    paste_combine = c()
    for(i in length(trace_ls):1){
      paste_combine = c(paste_combine, trace_ls[[i]])
    }
  }
  
  paste(paste_combine, collapse = " ")
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