#显示中文 
#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
#import library
{
  library(readr)
  #install.packages("rcdk")
  library(rcdk)
  library(igraph)
  #install.packages("fitdistrplus")
  library("fitdistrplus")
  library(tidyr)
  library(dplyr)
  library(enviPat)
}


##Read basic data information##
time = data.frame(Sys.time())
colnames(time)="read_data"
{
  
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
  output_csv=T
  
  #data1 records select data from raw table.
  HMDB_node_list = read_csv("HMDB_all_nodes.csv")
  data(isotopes)
  HMDB_node_list$MF=check_chemform(isotopes, HMDB_node_list$MF)$new_formula
  
  #HMDB_edge_list = read_csv("HMDB_edge_list.csv")
  df_raw = read.csv("yeast pos.csv")
  #df_raw = df_raw[df_raw$feature=="Metabolite"|df_raw$feature=="Adduct",]
  df_raw = df_raw[df_raw$feature=="Metabolite"|df_raw$feature=="[]",]
  #df_raw = df_raw[(!is.na(df_raw$Formula))|(!df_raw$high_blank),]
  
  #df_raw = sample_n(df_raw, 4000, replace = F)
  #df_raw = df_raw[1:2000,]
}
#Initialize data
mode = 1
{
  raw = df_raw
  colnames(raw)[1:3]=c("originID", "mz", "RT")
  raw = raw[with(raw, order(mz)),]
  raw["ID"]=1:nrow(raw)
  raw["MF"]=NA
  raw["category"]=2
  raw$category[raw$feature=="Metabolite"]=1
  raw["compound_name"]=NA
  
  # node_list = data.frame(1:nrow(raw),raw$mz, raw$MF, raw$name, raw$category, stringsAsFactors=FALSE)
  
  node_list = raw[,c("ID","mz","RT","MF", "category","compound_name")]
  
  raw_rnum = nrow(raw)
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  node_list$mz = node_list$mz*abs(mode) - (H_mass-e_mass)*mode
  
  #node_list=node_list[1:10000,]
}
#Initialize HMDB database and merge into targeted data
{
  #HMDB_edge_list[,1:2]=HMDB_edge_list[,1:2]+nrow(node_list)
  HMDB_node_list[,1]=HMDB_node_list[,1]
  HMDB_node_list["category"]=0
  HMDB_node_list["RT"]=0
  colnames(HMDB_node_list)[1:2]=c("ID", "mz")
  #HMDB_edge_list["category"]=0
  #HMDB_edge_list$mass_dif=0
  #HMDB_edge_list$edge_massdif_score=1
  
  #merge_node_list = rbind(node_list, HMDB_node_list)
  merge_node_list = HMDB_node_list
  #merge_node_list = merge_node_list[1:10000,]
  merge_nrow = nrow(merge_node_list)
}

##find mass difference corresponding to a functional group, such as "CH2"
time["edge_list"]=data.frame(Sys.time())


#isotope information
{isotope_info = list(
  H_iso=(1.00782503224-1),
  N_iso=(14.0030740048-14),
  O_iso=(15.99491462-16),
  Cl_iso=(34.96885268-35),
  F_iso=(18.99840322-19),
  Br_iso=(78.9183371-79),
  I_iso=(126.904473-127),
  S_iso=(31.97207100-32),
  P_iso=(30.97376163-31),
  Si_iso=(27.97692653-28),
  B_iso=(10.0129370-10),
  Na_iso=(22.9897692809-23),
  K_iso=(38.96370668-39),
  M_ELECTRON=0.00054857990943
)
}

#function group information
{
  fun_group_1 = data.frame(c("H2",
                             "C2H2O", #ACETYL
                             "C10H12N5O6P", #AMP
                             "C16H30O", #Palmityl
                             "C6H10O5", #Glycosyl
                             "C11H17NO9", #Sialic acid
                             "HPO3",
                             "CO", #FORMYL
                             "NH3",
                             #"NH3-O" ,#transaminase
                             #"OH-NH2", #amino group
                             #"C",
                             "CH2O",
                             #"C9H13N3O10P2", #CDP
                             "S",
                             "SO3",
                             "C2H4",
                             "O",
                             "NH",
                             "CH2",
                             "CO2",
                             "H2O"
  ),stringsAsFactors=F)
  
  data(isotopes)
  fun_group_1[,1] = check_chemform(isotopes,fun_group_1[,1])$new_formula
  fun_group_1["mass"] = check_chemform(isotopes,fun_group_1[,1])$monoisotopic_mass
  
  colnames(fun_group_1) = c("fun_group", "mass")
  
  fun_group_1[fun_group_1$fun_group=="NH3-O",2]=get.formula("NH3")@mass-get.formula("O")@mass
  fun_group_1[fun_group_1$fun_group=="OH-NH2",2]=get.formula("OH")@mass-get.formula("NH2")@mass
  fun_group_1=rbind(list("Same", 0),fun_group_1)
}



#Edge list
#mass_tol = 0.001
mass_ppm = 5/10^6
merge_node_list = merge_node_list[with(merge_node_list, order(mz)),]

timer=Sys.time()
{
  edge_ls = list()
  temp_mz_list=merge_node_list$mz
  
  for (k in 1:nrow(fun_group_1)){
    temp_fg=fun_group_1[k,2]
    i=j=1
    temp_edge_list = data.frame(node1=as.numeric(), node2=as.numeric(), linktype=as.numeric(), mass_dif=as.numeric())
    while(i<=merge_nrow){
      temp_ms=0
      mass_tol = temp_mz_list[i]*mass_ppm
      while(1){
        j=j+1
        if(j>merge_nrow){break}
        temp_ms = temp_mz_list[j]-temp_mz_list[i]
        
        if(abs(temp_ms-temp_fg)<mass_tol){
          temp_edge_list[nrow(temp_edge_list)+1,]=list(merge_node_list$ID[i], merge_node_list$ID[j], k, (temp_ms-temp_fg)/temp_mz_list[j]*1E6)
        }
        if((temp_ms-temp_fg)>mass_tol){break}
      }
      i=i+1
      
      while(j>i){
        j=j-1
        temp_ms = temp_mz_list[j]-temp_mz_list[i]
        if((temp_ms-temp_fg)<(-mass_tol)){break}
      }
    }
    edge_ls[[k]]=temp_edge_list
  }
}

print(Sys.time()-timer)
edge_list = bind_rows(edge_ls)
edge_list$linktype=fun_group_1$fun_group[edge_list$linktype]
#merge_edge_list = edge_list[edge_list$edge_massdif_score>0.7,]

# test_j2 = edge_list
# test_ji = edge_list
# test_intersect=intersect(test_ji,test_j2)
# test_dif = setdiff(test_ji,test_intersect)
# test_unique=unique(test_j2)
# 
# edge_list_sub = subset(edge_list, 
#                        (edge_list$node1<=nrow(node_list)|
#                           edge_list$node2<=nrow(node_list))&
#                          edge_list$node1!=edge_list$node2
# )
edge_list_sub = edge_list
edge_list_sub["category"]=1


#Scoring edge based on mass accuracy
{
  edge_mzdif_FIT <- fitdist(as.numeric(edge_list_sub[,4])*10^3, "norm")    
  hist(edge_list_sub$mass_dif)
  
  plot(edge_mzdif_FIT)  
  edge_list_sub["edge_massdif_score"]=dnorm(edge_list_sub[,4]*10^3, edge_mzdif_FIT$estimate[1], edge_mzdif_FIT$estimate[2])
  edge_list_sub["edge_massdif_score"]=edge_list_sub["edge_massdif_score"]/max(edge_list_sub["edge_massdif_score"])
}

# node_list_with_edge = unique(c(edge_list_sub$node1,edge_list_sub$node2))
# node_list_with_edge = node_list_with_edge[node_list_with_edge<=nrow(node_list)]
if(output_csv){
  # write.csv(edge_list_sub,"formula_network_edge_list2.csv",row.names = F)
  # write.csv(raw[,c("ID","originID")], "ID.csv",row.names = F)
}

##Predict formula based on known formula and edgelist 
time["formula_pred"]=Sys.time()

#formula_manipulate is used to process two formula, calling subfunction formula_table and table_to_formula
{
  formula_table = function(X){
    Xformula = get.formula(X)
    temp=Xformula@isotopes
    df = data.frame(as.numeric(temp[,2]))
    rownames(df)=temp[,1]
    colnames(df)=Xformula@string
    return(df)
  }
  table_to_formula = function(result_formula){
    r = as.character()
    for(i in 1:nrow(result_formula)){
      if(result_formula$result[i]!=0){
        r=paste(r, result_formula$Row.names[i],  result_formula$result[i], sep="")
      }
    }
    if(sum(result_formula$result<0)!=0){
      r=paste(r, "Illegal_formula")
    }
    return (r)
  }
  formula_manipulate = function(F1, F2, plusminus) {
    temp_1 = formula_table(F1)
    temp_2 = formula_table(F2)
    result_formula = merge(temp_1,temp_2, all=T,  by="row.names")
    result_formula[is.na(result_formula)]=0
    result_formula["result"]=result_formula[,2]+result_formula[,3]*plusminus
    return (table_to_formula(result_formula))
  }
  ##My_mergefrom&My_subfrom
  #formula1<-c("C10H8O3Cl0")
  #formula2<-c("C2H5O2")
  My_mergefrom=function (formula1, formula2)
  {
    formula2 <- gsub("D", "[2]H", formula2)
    ende2 <- nchar(formula2)
    element2 <- c()
    number2 <- c()
    j <- c(1)
    while (j <= ende2) {
      #browser()
      if (substr(formula2, j, j) == c("[")) {
        b <- j
        while (any(substr(formula2, j, j) == c("]")) != 
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula2, j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element2 <- c(element2, substr(formula2, b, m))
      }
      if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                          "4", "5", "6", "7", "8", "9")) != TRUE) {
        k <- j
        while (any(substr(formula2, j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element2 <- c(element2, substr(formula2, k, m))
      }
      if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                          "4", "5", "6", "7", "8", "9")) == TRUE) {
        k <- j
        while (any(substr(formula2, j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) == 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number2 <- c(number2, as.numeric(substr(formula2, 
                                                k, m)))
      }
      j <- j + 1
    }
    formulas <- c()
    for (i in 1:length(formula1)) {
      #i=1
      formula1[i] <- gsub("D", "[2]H", formula1[i])
      ende1 <- nchar(formula1[i])
      element1 <- c()
      number1 <- c()
      j <- c(1)
      while (j <= ende1) {
        if (substr(formula1[i], j, j) == c("[")) {
          b <- j
          while (any(substr(formula1[i], j, j) == c("]")) != 
                 TRUE) {
            j <- c(j + 1)
          }
          k <- j
          while (any(substr(formula1[i], j, j) == c("0", 
                                                    "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          element1 <- c(element1, substr(formula1[i], 
                                         b, m))
        }
        if (any(substr(formula1[i], j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) != 
            TRUE) {
          k <- j
          while (any(substr(formula1[i], j, j) == c("0", 
                                                    "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          element1 <- c(element1, substr(formula1[i], 
                                         k, m))
        }
        if (any(substr(formula1[i], j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) == 
            TRUE) {
          k <- j
          while (any(substr(formula1[i], j, j) == c("0", 
                                                    "1", "2", "3", "4", "5", "6", "7", "8", "9")) == 
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          number1 <- c(number1, as.numeric(substr(formula1[i], 
                                                  k, m)))
        }
        j <- j + 1
      }
      both <- unique(c(element1, element2))
      both = both[order(both)]
      counts <- c()
      for (i in 1:length(both)) {
        if (any(element1 == both[i])) {
          it1 <- c(number1[element1 == both[i]])
        }
        else {
          it1 <- c(0)
        }
        if (any(element2 == both[i])) {
          it2 <- c(number2[element2 == both[i]])
        }
        else {
          it2 <- c(0)
        }
        counts <- c(counts, it1 + it2)
      }
      formula_all <- ""
      for (i in 1:length(both)) {
        if(counts[i]==0){next}
        formula_all <- paste(formula_all, both[i], counts[i], 
                             sep = "")
      }
      formulas <- c(formulas, formula_all)
    }
    return(formulas)
  }
  #My_mergefrom(formula1,formula2)
  
  My_subfrom =function (formula1, formula2) 
  {
    formula2 <- gsub("D", "[2]H", formula2)
    ende2 <- nchar(formula2)
    element2 <- c()
    number2 <- c()
    j <- c(1)
    while (j <= ende2) {
      if (substr(formula2, j, j) == c("[")) {
        b <- j
        while (any(substr(formula2, j, j) == c("]")) != 
               TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula2, j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element2 <- c(element2, substr(formula2, b, m))
      }
      if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                          "4", "5", "6", "7", "8", "9")) != TRUE) {
        k <- j
        while (any(substr(formula2, j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) != 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element2 <- c(element2, substr(formula2, k, m))
      }
      if (any(substr(formula2, j, j) == c("0", "1", "2", "3", 
                                          "4", "5", "6", "7", "8", "9")) == TRUE) {
        k <- j
        while (any(substr(formula2, j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) == 
               TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number2 <- c(number2, as.numeric(substr(formula2, 
                                                k, m)))
      }
      j <- j + 1
    }
    formulas <- c()
    for (i in 1:length(formula1)) {
      formula1[i] <- gsub("D", "[2]H", formula1[i])
      ende1 <- nchar(formula1[i])
      element1 <- c()
      number1 <- c()
      j <- c(1)
      while (j <= ende1) {
        if (substr(formula1[i], j, j) == c("[")) {
          b <- j
          while (any(substr(formula1[i], j, j) == c("]")) != 
                 TRUE) {
            j <- c(j + 1)
          }
          k <- j
          while (any(substr(formula1[i], j, j) == c("0", 
                                                    "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          element1 <- c(element1, substr(formula1[i], 
                                         b, m))
        }
        if (any(substr(formula1[i], j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) != 
            TRUE) {
          k <- j
          while (any(substr(formula1[i], j, j) == c("0", 
                                                    "1", "2", "3", "4", "5", "6", "7", "8", "9")) != 
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          element1 <- c(element1, substr(formula1[i], 
                                         k, m))
        }
        if (any(substr(formula1[i], j, j) == c("0", "1", 
                                               "2", "3", "4", "5", "6", "7", "8", "9")) == 
            TRUE) {
          k <- j
          while (any(substr(formula1[i], j, j) == c("0", 
                                                    "1", "2", "3", "4", "5", "6", "7", "8", "9")) == 
                 TRUE) {
            j <- c(j + 1)
          }
          m <- c(j - 1)
          j <- c(j - 1)
          number1 <- c(number1, as.numeric(substr(formula1[i], 
                                                  k, m)))
        }
        j <- j + 1
      }
      formula_all <- TRUE
      for (i in 1:length(element2)) {
        if (any(element2[i] == element1) == FALSE) {
          return(F)
          formula_all <- paste(element2[i], " from formula 2 not part of formula1", 
                               sep = "")
        }
        else {
          if (number2[i] > number1[element2[i] == element1]) {
            return(F)
            formula_all <- paste("Atom number of ", element2[i], 
                                 " from formula 2 not fully subset of formula1 atom number", 
                                 sep = "")
          }
        }
      }
      number1=number1[order(element1)]
      element1=element1[order(element1)]
      if (formula_all == TRUE) {
        formula_all <- ""
        counts <- c()
        for (i in 1:length(element1)) {
          if (any(element2 == element1[i])) {
            counts <- c(counts, number1[i] - (number2[element2 == 
                                                        element1[i]]))
          }
          else {
            counts <- c(counts, number1[i])
          }
        }
        element1 <- element1[counts > 0]
        counts <- counts[counts > 0]
        for (i in 1:length(counts)) {
          formula_all <- paste(formula_all, element1[i], 
                               counts[i], sep = "")
        }
      }
      formulas <- c(formulas, formula_all)
    }
    return(formulas)
  }
  #My_subfrom(formula1, formula2)
  
}

#Initialize predict_formula from HMDB known formula
{
  data_str = data.frame(id=as.numeric(),
                        formula=as.character(), 
                        steps=as.numeric(), 
                        parent=as.numeric(), 
                        score=as.numeric(), stringsAsFactors = F)
  #sf stands for summary_formula
  merge_node_list = merge_node_list[with(merge_node_list, order(ID)),]
  mnl=merge_node_list
  sf = lapply(1:nrow(merge_node_list),function(i)(data_str))

  
  for(i in (1):1){
    Initial_formula =sf[[i]]
    Initial_formula[1,]= list(i,mnl$MF[i],0,1, 1)
    sf[[i]]=Initial_formula

  }
}

#while loop to Predict formula based on known formula and edgelist 

top_formula_n = 3
New_nodes_in_network = 1
step=0
timer=Sys.time()
while(New_nodes_in_network==1){
  
  all_nodes_df = bind_rows(lapply(sf, head,top_formula_n))
  new_nodes_df = all_nodes_df[all_nodes_df$steps==step 
                              &all_nodes_df$score>0,]
  print(paste("nrow",nrow(all_nodes_df),"in step",step,"elapsed="))
  print((Sys.time()-timer))
  
  if(nrow(new_nodes_df)==0){break}
  
  #Handling head
  {
    edge_list_node1 = edge_list_sub[edge_list_sub$node1 %in% new_nodes_df$id,]
    head_list = edge_list_node1$node1
    
    for (n in 1: nrow(new_nodes_df)){
      head = new_nodes_df$id[n]
      head_formula = new_nodes_df$formula[n]
      
      temp_edge_list=subset(edge_list_node1, edge_list_node1$node1==head)
      if(nrow(temp_edge_list)==0){next}
      #head_score = sf[[head]]$score[min(which(sf[[head]]$formula==head_formula))]
      head_score = sf[[head]]$score[match(head_formula,sf[[head]]$formula)]
      
      if(head_score==0){next}
      
      for(i in 1:nrow(temp_edge_list)){
        tail=temp_edge_list$node2[i]
        temp_fg = temp_edge_list$linktype[i]
        if(temp_fg=="Same"){temp_formula=head_formula}
        #else{temp_formula = get.formula(paste(head_formula,"+", temp_fg))@string}
        else{temp_formula = My_mergefrom(head_formula,temp_fg)}
        temp_steps = step+1
        temp_parent = head
        temp_score = head_score*temp_edge_list$edge_massdif_score[i]
        #Criteria to enter new entry into formula list
        #1. new formula
        temp = sf[[tail]]
        temp_subset=subset(temp, temp$formula==temp_formula)
        if(nrow(temp_subset)!=0){
          #2. much higher scores
          if(temp_score<=(1.2*max(temp_subset$score))){
            next
          }
        }
        #Enter new entry
        temp[nrow(temp)+1, ] = list(tail, temp_formula, temp_steps, temp_parent, temp_score)
        temp = temp[with(temp, order(-score)),]
        sf[[tail]]=temp
      }
    }
  }
  
  #Handling tail
  {
    edge_list_node2 = edge_list_sub[edge_list_sub$node2 %in% new_nodes_df$id,]
    tail_list = edge_list_node2$node2
    
    for (n in 1: nrow(new_nodes_df)){
      tail = new_nodes_df$id[n]
      tail_formula = new_nodes_df$formula[n]
      
      temp_edge_list=subset(edge_list_node2, edge_list_node2$node2==tail)
      if(nrow(temp_edge_list)==0){next}
      
      #tail_score = sf[[tail]]$score[min(which(sf[[tail]]$formula==tail_formula))]
      tail_score = sf[[tail]]$score[match(tail_formula,sf[[tail]]$formula)]
      
      if(tail_score==0){next}
      
      for(i in 1:nrow(temp_edge_list)){
        head=temp_edge_list$node1[i]
        temp_fg = temp_edge_list$linktype[i]
        if(temp_fg=="Same"){temp_formula=tail_formula}
        else{temp_formula = My_subfrom(tail_formula, temp_fg)}
        #else{temp_formula = formula_manipulate(tail_formula, temp_fg, -1)}
        if(is.logical(temp_formula))#function return false if not valid formula
        { temp_score=0
        temp_formula=paste(tail_formula,temp_fg,sep="-")
        }
        else{temp_score = tail_score*temp_edge_list$edge_massdif_score[i]}
        temp_steps = step+1
        temp_parent = tail
        #Criteria to enter new entry into formula list
        #1. new formula
        temp = sf[[head]]
        temp_subset=subset(temp, temp$formula==temp_formula)
        if(nrow(temp_subset)!=0){
          #2. much higher scores
          if(temp_score<=(1.2*max(temp_subset$score))){
            next
          }
        }
        #Enter new entry
        temp[nrow(temp)+1, ] = list(head, temp_formula, temp_steps, temp_parent, temp_score)
        temp = temp[with(temp, order(-score)),]
        sf[[head]]=temp
        
      }
    }
  }
  
  step=step+1
  
}

#merge_edge_list = rbind(edge_list_sub, HMDB_edge_list)
merge_edge_list = edge_list_sub

#Update the mnl from [List]
merge_node_list["Score"]=NA
merge_node_list["Predict_formula"]=NA
for(i in 1:nrow(merge_node_list)){
  merge_node_list$Predict_formula[[i]]=sf[[i]]$formula[1]
  merge_node_list$Score[[i]]=sf[[i]]$score[1]
}

merge_node_list["same"]=(merge_node_list$Predict_formula==merge_node_list$MF)

merge_node_list_false = merge_node_list[!merge_node_list$same,]
merge_node_list_pred = merge_node_list[!is.na(merge_node_list$Predict_formula),]
merge_node_list_same = merge_node_list_pred[merge_node_list_pred$same,]

test = merge_node_list[merge_node_list$degree==0,]

# test=sample_n(All_formula_predict,1000)
# for(i in 1:1000){
#   if(!grepl("-",test$formula[i]))
#   test[i,6]=isvalid.formula(get.formula(test$formula[i]))
#   
# }




All_formula_predict=bind_rows(sf)
which.max(table(All_formula_predict$id))


time["end"]=Sys.time()
time_used=as.numeric()
for (i in 2:ncol(time)){
  time_used[i-1] = round(time[1,i]-time[1,i-1],2)
}
names(time_used)=colnames(time)[1:(ncol(time)-1)]
time_used

if(output_csv){
  write.csv(merge_node_list, "merge_node_list.csv",row.names = F)
  write.csv(All_formula_predict, "All_formula_predict.csv",row.names = F)
  write.csv(merge_edge_list,"merge_edge_list.csv", row.names = F)
}

read_from_csv = 0
if(read_from_csv){
  merge_node_list=read.csv("merge_node_list.csv")
  merge_edge_list=read.csv("merge_edge_list.csv")
  edge_list_sub=read.csv("edge_list_sub.csv")
}

##Network analysis
time["network_analysis"]=Sys.time()
#Generate network graph using node_list and edge_list above

g <- graph_from_data_frame(d = merge_edge_list, vertices = merge_node_list, directed = FALSE)
colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")
V(g)$color = colors[merge_node_list$category+1]
E(g)$color = colors[merge_edge_list$category+1]
merge_node_list["degree"]=degree(g, mode = c("all"))

g_sub = graph_from_data_frame(d = edge_list_sub, vertices = merge_node_list[merge_node_list$ID %in% c(edge_list_sub[,1], edge_list_sub[,2]),], directed = FALSE)
colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")
V(g_sub)$color = colors[vertex.attributes(g_sub)$category+1]
E(g_sub)$color = colors[edge_list_sub$category+1]

g_correct = graph_from_data_frame(d = edge_list_sub[edge_list_sub$node1%in%merge_node_list_same$ID
                                                    &edge_list_sub$node2%in%merge_node_list_same$ID,], 
                                  vertices = merge_node_list_same, directed = FALSE)
colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")
V(g_correct)$color = colors[vertex.attributes(g_correct)$category+1]
E(g_correct)$color = colors[edge_list_sub[edge_list_sub$node1%in%merge_node_list_same$ID
                                          &edge_list_sub$node2%in%merge_node_list_same$ID,"category"]+1]

#Plot TCA related graph
{
  merge_node_list_same[grep("C6H6O6", merge_node_list_same$MF),]
  
  g_intrest <- make_ego_graph(g_sub,1, nodes = "45", mode = c("all"))[[1]]
  
  plot(g_intrest,
       #vertex.color = 'white',
       #vertex.label = vertex.attributes(g_intrest)$Predict_formula,
       #vertex.label = vertex.attributes(g_intrest)$MF,
       #vertex.label = vertex.attributes(g_intrest)$RT,
       #vertex.label = vertex.attributes(g_intrest)$mz,
       vertex.label = vertex.attributes(g_intrest)$ID,
       #vertex.label = vertex.attributes(g_intrest)$originID,
       vertex.label.color = "black",
       vertex.label.cex = 1,
       #edge.color = 'black',
       edge.label = edge.attributes(g_intrest)$linktype,
       vertex.size = 10,
       edge.arrow.size = 0.05
       
  )
}


graph.union()


#Basic graph characteristics, distance, degree, betweeness
{
  #distance
  farthest_vertices(g) 
  #Degree
  g.degree <- degree(g, mode = c("all"))
  table(g.degree)
  hist(g.degree)
  # which.max(g.degree)
  merge_node_list_0degree = merge_node_list[g.degree==0,]
  # node_list$mz[g.degree==0]
  # node_list$mz[71]
  # node_list$MF[71]
  #Betweenness
  g.b <- betweenness(g, directed = T)
  colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")
  V(g)$color <- colors[merge_node_list$category+1]
  plot(g_correct,
       vertex.label = NA,
       #edge.color = 'black',
       vertex.size = sqrt(degree(g_sub, mode = c("all"))),
       edge.arrow.size = 0.05,
       layout = layout_nicely(g_correct))
}

target_mz = 645.09594727
ppm = 10/10^6
merge_node_list[abs(merge_node_list$mz-target_mz+(H_mass-e_mass)*mode)<target_mz*ppm,]
#Analyze the network/subgraph of specific node
{
  interested_node = "174"
  g_intrest <- make_ego_graph(g_sub,1, nodes = 1:3, mode = c("all"))[[1]]
  test=data.frame(vertex.attributes(g_sub))
  #dists = distances(g_intrest, interested_node)
  colors <- c("black", "red", "orange", "blue", "dodgerblue", "cyan")
  #V(g_intrest)$color <- colors[dists+1]
  png(filename=paste("Subnetwork of mz=", target_mz,".png",sep=""),
      width = 2400, height=2400,
      res=300)
  plot(g_intrest,
       #vertex.color = 'white',
       vertex.label = vertex.attributes(g_intrest)$Predict_formula,
       #vertex.label = vertex.attributes(g_intrest)$MF,
       #vertex.label = vertex.attributes(g_intrest)$RT,
       #vertex.label = vertex.attributes(g_intrest)$mz,
       #vertex.label = vertex.attributes(g_intrest)$ID,
       #vertex.label.color = "black",
       vertex.label.cex = 1,
       #edge.color = 'black',
       edge.label = edge.attributes(g_intrest)$linktype,
       vertex.size = 5,
       edge.arrow.size = 0.05,
       main = paste("Subnetwork of mz", interested_node)
  )
  dev.off()
  plot(g_intrest,
       #vertex.color = 'white',
       #vertex.label = vertex.attributes(g_intrest)$Predict_formula,
       #vertex.label = vertex.attributes(g_intrest)$MF,
       #vertex.label = vertex.attributes(g_intrest)$RT,
       #vertex.label = vertex.attributes(g_intrest)$mz,
       vertex.label = vertex.attributes(g_intrest)$ID,
       #vertex.label = vertex.attributes(g_intrest)$originID,
       vertex.label.color = "black",
       vertex.label.cex = 1,
       #edge.color = 'black',
       edge.label = edge.attributes(g_intrest)$linktype,
       vertex.size = 10,
       edge.arrow.size = 0.05,
       main = paste("Subnetwork of mz", interested_node)
  )
}
merge_node_list[merge_node_list$ID==1862,]
test=data.frame(vertex.attributes(g_intrest))


output_network_csv = merge_node_list[vertex.attributes(g_intrest)$compound_name,]
output_network_csv$mz=output_network_csv$mz +(H_mass-e_mass)*mode
write.csv(output_network_csv, paste("network_",target_mz,".csv",sep=""), row.names=F)



merge_node_list[which(merge_node_list$MF=="C6H13O4PS2"),]
merge_node_list[which(merge_node_list$Predict_formula=="C15H4O5"),]

#Analysis of Subnetwork 
{
  clu=components(g_sub)
  #subnetwork criteria 
  subnetwork = igraph::groups(clu)[table(clu$membership)>=2&table(clu$membership)<100]
  g_subnetwork_list = lapply(subnetwork, make_ego_graph, graph=g_sub, order=diameter(g_sub), mode="all")
  for (i in 1:length(subnetwork)){
    plot(g_subnetwork_list[[i]][[1]],
         #vertex.color = 'white',
         vertex.label = vertex.attributes(g_subnetwork_list[[i]][[1]])$MF,
         #vertex.label = vertex.attributes(g_subnetwork_list[[i]][[1]])$mz,
         vertex.label.color = "red",
         vertex.label.cex = 1,
         #edge.color = 'black',
         edge.label = edge.attributes(g_subnetwork_list[[i]][[1]])$linktype,
         vertex.size = 10,
         edge.arrow.size = .05,
         main = paste("Subnetwork",names(subnetwork)[[i]])
    )
  }
  merge_node_list[subnetwork[["233"]],]
}

print(merge_node_list[merge_node_list$ID==9750,])

Match_HMDB = function(mz, data_known, mode){
  H_mass = 1.00782503224
  e_mass = 0.00054857990943
  temp_mass = mz*abs(mode) - (H_mass-e_mass)*mode
  temp_table = data_known[abs(data_known$mz-temp_mass)<0.002,]
  return(temp_table)
}


# #Interactive network w/ threejs
# install.packages("htmlwidgets ")
# library(threejs)
# library(htmlwidgets)
# 
# graphjs(g, main="Network!", bg="gray10", showLabels=F, stroke=F, 
#         curvature=0.1, attraction=0.9, repulsion=0.8, opacity=0.9)


time["end"]=Sys.time()
time_used=as.numeric()
for (i in 2:ncol(time)){
  time_used[i-1] = round(time[1,i]-time[1,i-1],2)
}
names(time_used)=colnames(time)[1:(ncol(time)-1)]
time_used




#Calculate mz difference
# {
# merge_node_list["MF_measured_mz_dif"]=as.numeric()
# merge_node_list["Predict_measured_mz_dif"]=as.numeric()
# for(i in 1: nrow(node_list)){
#   if(!is.na(merge_node_list$`MF`[i])){
#     merge_node_list$`MF_measured_mz_dif`[i]=get.formula(merge_node_list$MF[i])@mass-merge_node_list$mz[i]
#   }
#   if((!is.na(merge_node_list$`Predict_formula`[i]))&(!grepl("Illegal_formula", merge_node_list$Predict_formula[i]))){
#     merge_node_list$`Predict_measured_mz_dif`[i]=get.formula(merge_node_list$Predict_formula[i])@mass-merge_node_list$mz[i]
#   }
# }
# }
# #Plot machine measurement difference
# {
# cat1 = merge_node_list[merge_node_list$category==1,]
# cat2 = merge_node_list[merge_node_list$category==2,]
# pred_same_cat1 = cat1[(cat1$MF==cat1$Predict_formula)&(!is.na(cat1$Predict_formula)),]
# hist(abs(pred_same_cat1$MF_measured_mz_dif)*10^3)
# hist(abs(pred_same_cat1$MF_measured_mz_dif/pred_same_cat1$mz)*10^6)
# }
# 
# #Debug & test
# {
# merge_node_list_pred=merge_node_list[complete.cases(merge_node_list[,c(3,6)]),]
# good_node = merge_node_list_pred[merge_node_list_pred$MF==merge_node_list_pred$Predict_formula,]
# {
#   goodnode_FIT <- fitdist(as.numeric(good_node[,9])*10^3, "norm")    
#   plot(goodnode_FIT)  
#   edge_list_sub["edge_massdif_score"]=dnorm(edge_list_sub[,4]*10^3, edge_mzdif_FIT$estimate[1], edge_mzdif_FIT$estimate[2])
#   edge_list_sub["edge_massdif_score"]=edge_list_sub["edge_massdif_score"]/max(edge_list_sub["edge_massdif_score"])
# }
# 
# data_analysis = merge_node_list[merge_node_list$category!=0&!(is.na(merge_node_list$MF)),]
# data_analysis = data_analysis[with(data_analysis, order(`MF_measured_mz_dif`)),]
# test_logic=data_analysis$MF!=data_analysis$Predict_formula
# test_logic[is.na(test_logic)]=T
# data_analysis_dif = data_analysis[test_logic,]
# 
# sf[[merge_node_list$ID[which(merge_node_list$Predict_formula=="C10H21O7PS2")]]]
# sf[[6578]]
# merge_node_list[merge_node_list$ID==3796,]
# merge_node_list[merge_node_list$ID==3542,]
# merge_node_list[merge_node_list$ID==2563,]
# }
# 
# #write.csv(merge_node_list[merge_node_list$category==1|merge_node_list$category==2,], "neg_unknown_pred.csv")
