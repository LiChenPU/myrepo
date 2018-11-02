#显示中文 
# !diagnostics off
#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
#import library
{
  library(readr)
  #install.packages("rcdk")
  library(rcdk)
  library(igraph)
  #install.packages("fitdistrplus")
  library("fitdistrplus")
  #install.packages("truncnorm")
  library(truncnorm)
  library(mixdist)  
  #install.packages("VennDiagram")
  library(VennDiagram)
  library(dplyr)
  library(tidyr)
  #install.packages("tibble")
  library(tibble)
}

test_time = Sys.time()
##Read basic data information##
time = data.frame(Sys.time())
colnames(time)="read_data"
{
#setwd("C:/Users/lc8/Dropbox/HRMass_ID/")
setwd("C:/Users/Li Chen/Dropbox/HRMass_ID/180629 Triwizard Tournament/JunkMS_ID")

#data1 records select data from raw table.

df_raw = read_csv("hmdb_Yeast-Ecoli-neg-peakpicking_blank.csv")
#df_raw = sample_n(df_raw, 4000, replace = F)
#df_raw = df_raw[1:000,]

#hist(df_raw$medRt)
# df_raw = cbind(id=1:nrow(df_raw), df_raw[with(df_raw, order(medRt)),])
# df_raw = df_raw[,colnames(df_raw)!=".log_p_FDR_corrected"]


sample_names = colnames(df_raw)[8:ncol(df_raw)]
if(length(grep("blank|blk", sample_names, ignore.case = T))!=0){
  sample_names_noblank=sample_names[-grep("blank|blk", sample_names, ignore.case = T)]
} else {
  sample_names_noblank=sample_names
}
sample_names_blank=sample_names[grep("blank|blk", sample_names, ignore.case = T)]
sample_cohort=stri_replace_last_regex(sample_names,'-1|-2|-3|-a|-b|-c|_mean', '',stri_opts_regex(case_insensitive=T))

df_raw=add_column(df_raw, mean_inten=rowMeans(df_raw[,sample_names_noblank]), .after = 7)
df_raw=add_column(df_raw, log10_inten=(log10(df_raw$mean_inten)), .after = 7)
df_raw=add_column(df_raw, artifact=NA, .after = 7)
# df_raw["mean_inten"]mean_inten=rowMeans(df_raw[,sample_names_noblank])
# df_raw["log10_inten"]=log10(df_raw["mean_inten"])
}

# test = df_raw
# 
# test_hasformula = test[!is.na(test$Formula),]
# test_lowblank = test[!test$high_blank,]
# test_significantpeaks = test[test$`_log10_FDR`!=0,]
# 
# test_output=unique(rbind(test_hasformula,test_lowblank,test_significantpeaks))
# 
# df_raw = test_output

df_raw = df_raw[,!colnames(df_raw) %in% sample_names_blank]


##Find adjacent peaks and calculate mass differece and intensity correlation
time_cutoff=0.1
time["edge_list"]=data.frame(Sys.time())
{
{
  df_raw = df_raw[with(df_raw, order(medRt)),]
  i_min = i_max =1
  #edge_ls = data.frame(node1=as.numeric(), node2=as.numeric(), time_dif=as.numeric(), mz_dif=as.numeric(),correlation=as.numeric())
  edge_list = list()
  
  while(i_min!= nrow(df_raw)){
    temp_t = df_raw$medRt[i_min]
    temp_t_end = temp_t+time_cutoff
    
    while(df_raw$medRt[i_max]<temp_t_end){
      if(i_max == nrow(df_raw)){break}
      i_max = i_max+1
    }
    i_max = i_max-1
    
    temp_df_raw = df_raw[i_min:i_max,]
    
    temp_df_raw$time_dif=temp_df_raw$medRt-temp_df_raw$medRt[1]
    temp_df_raw$mz_dif = round(temp_df_raw$medMz-temp_df_raw$medMz[1], digits=5)
    
    temp_x = t(temp_df_raw[1,sample_names_noblank])
    temp_y = t(temp_df_raw[,sample_names_noblank])
    temp_df_raw$correlation = as.numeric(t(cor((temp_x), (temp_y))))
     
    # test_x = matrix(c(1,1,1,1000,1000,1005,1,1,1,100,100,105), ncol=1) 
    # test_y = matrix(c(0,0,0,24,23,22,1,1,1,10,10,15), ncol=1) 
    # cor((test_x), (test_y), method = "pearson")
    
    
    temp_df_raw["node1"]=temp_df_raw$ID[1]
    temp_df_raw["node2"]=temp_df_raw$ID
    
    
    temp_df_raw["mz_node1"] = temp_df_raw$medMz[1]
    temp_df_raw["mz_node2"] = temp_df_raw$medMz
    
    temp_df_raw["log10_inten_node1"] = temp_df_raw$log10_inten[1]
    temp_df_raw["log10_inten_node2"] = temp_df_raw$log10_inten
    temp_df_raw["log10_inten_ratio"] = temp_df_raw["log10_inten_node2"]-temp_df_raw["log10_inten_node1"]
    
    temp_edge_ls=temp_df_raw[,c("node1","node2","time_dif", "mz_dif", "correlation","mz_node1",
                                "mz_node2","log10_inten_node1","log10_inten_node2",
                                "log10_inten_ratio")]
    
    # cor(temp_x, temp_y, method = "kendall")
    # cor(temp_x, temp_y, method = "spearman")
    # cosine(tem#p_x, temp_y)
    temp_edge_ls=temp_edge_ls[temp_edge_ls$correlation>0.8,]
    edge_list[[i_min]]=temp_edge_ls
    i_min = i_min+1
  }
}
edge_ls = bind_rows(edge_list)
edge_ls=edge_ls[edge_ls$node1!=edge_ls$node2,]
rm(edge_list)
}

{
  edge_ls = edge_ls[with(edge_ls, order(mz_dif)),]
  temp_data = edge_ls[edge_ls$mz_dif<0,]
  
  temp_data$mz_node1=temp_data$mz_node1+temp_data$mz_dif
  temp_node=temp_data$node1
  temp_data$node1=temp_data$node2
  temp_data$node2=temp_node
  temp_data$time_dif=-temp_data$time_dif
  temp_data$mz_dif=-temp_data$mz_dif
  temp_data$log10_inten_ratio=-temp_data$log10_inten_ratio
  
  edge_ls[edge_ls$mz_dif<0,] = temp_data 
  rm(temp_data)
}



# write.csv(edge_ls,"edge_ls.csv", row.names = F)


##Select only high correlation peaks
edge_ls = edge_ls[with(edge_ls, order(mz_dif)),]
edge_ls_highcor = edge_ls[edge_ls$correlation>0.6,]

##label mz_dif (edge) with natural abundance peaks, adducts, fragments, oligomers and multi-charge
time["edge_annotation_isotope"]=data.frame(Sys.time())
{
#Define junk dataframe used to identify natural abundance peaks, adducts and fragments
{
  isotope_table=list(c("13C", 13.0033548378-12),
                     c("2H", 2.01410177811-1.00782503224),
                     c("18O", 17.9991610-15.99491461956),
                     c("15N", 15.0001088982-14.0030740048),
                     c("34S", 33.96786690-31.97207100),
                     c("37Cl", 36.96590259-34.96885268),
                     #c("53Cr", 52.9406494-51.9405075),
                     #c("65Cu", 64.9277895-62.9295975),
                     #c("54Cr", 53.9388804-51.9405075),
                     c("41K", 40.96182576-38.96370668)
  )
  natural_abundance <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("category", "type", "mass")
  colnames(natural_abundance) <- x
  for(i in 1:length(isotope_table)){
    natural_abundance[nrow(natural_abundance)+1,]=c("natural_abundance", isotope_table[[i]][1],isotope_table[[i]][2])
  }
  
  
  adduct = data.frame(type=c("H3PO4", 
                             "C2H4O2",
                             "H2SO4",
                             "H4O4Si",
                             "C16H32O2",
                             "C18H36O2",
                             "C2HF3O2",
                             "CrO3",
                             "HCl",
                             "HCOOH"
  ),                           
  stringsAsFactors=F)
  adduct["category"]="adduct"
  adduct["mass"]=NA
  for (i in 1:nrow(adduct)){
    adduct$mass[i]=get.formula(adduct[i,1])@mass
  }
  fragment = data.frame(type=c("CO2", 
                               "H2O", 
                               "NH3"
  ),                           
  stringsAsFactors=F)
  fragment["category"]="fragment"
  fragment["mass"]=NA
  
  for (i in 1:nrow(fragment)){
    fragment$mass[i]=get.formula(fragment[i,1])@mass
  }
  H_replacement=data.frame(type=c("K", 
                                  "Na"
  ),                           
  stringsAsFactors=F)
  H_replacement["category"]="H_replacement"
  H_replacement["mass"]=NA
  for (i in 1:nrow(H_replacement)){
    H_replacement$mass[i]=get.formula(H_replacement[i,1])@mass-1.00782503224
  }
  
  junk_df = rbind.data.frame(natural_abundance, adduct, fragment, H_replacement)
  junk_df$category=as.factor(junk_df$category)
  junk_df$mass=as.numeric(junk_df$mass)
  
}
summary(junk_df)
junk_df=junk_df[with(junk_df, order(mass)),]

#Statsitically find the mass range for natural abundance peaks, adducts and fragments mass
{
  test_13C_cor=edge_ls_highcor[(edge_ls_highcor$mz_dif<(1.0033548+0.001))&(edge_ls_highcor$mz_dif>(1.0033548-0.001)),]
  initial_FIT <- fitdist(as.numeric(test_13C_cor$mz_dif)*10^3, "norm")    
  hist(test_13C_cor$mz_dif)
  plot(initial_FIT)  
  initial_FIT$estimate[1]
  initial_FIT$estimate[2]
}

#Match mass_dif to natural abundance peaks, adducts and fragments mass
search_ms_cutoff = 10 * initial_FIT$estimate[2]/1000
{
i=j=1
temp_ls = list()
temp_df = edge_ls_highcor[1,]
temp_df["category"]=NA
temp_df["type"]=NA
temp_df = temp_df[0,]

while (i <= nrow(edge_ls_highcor)){
  i=i+1
  if(edge_ls_highcor$mz_dif[i]<(junk_df$mass[j]-search_ms_cutoff)){
    next
  }
  if(edge_ls_highcor$mz_dif[i]>(junk_df$mass[j]+search_ms_cutoff)){
    temp_ls[[j]]=temp_df
    temp_df = temp_df[0,]
    j=j+1
    if(j>nrow(junk_df)){break}
    search_ms_cutoff = 10 * initial_FIT$estimate[2]/1000 
    while(edge_ls_highcor$mz_dif[i]>(junk_df$mass[j]-search_ms_cutoff)){i=i-1}
    next
  }
  temp_df[(nrow(temp_df)+1),]=c(edge_ls_highcor[i,],as.character(junk_df$category[j]), junk_df$type[j])
}
}

temp_df_isotope = bind_rows(temp_ls)

junk_summary=table(temp_df_isotope$type)
junk_summary

}


##Handle multiple-charge and oligomer species

{
#Oligomers. Note that it is indistinguishable between 2-charge-parent pair and parent-dimer pair
  time["edge_annotation_oligomer"]=data.frame(Sys.time())
  test_time = Sys.time()
  {
H_mass = 1.00782503224
e_mass = 0.00054857990943
mode=-1
ppm=2/10^6
temp_df_oligo = temp_df[0,]

temp_edge_ls = data.frame(ratio=edge_ls_highcor$mz_dif/edge_ls_highcor$mz_node1)
temp_edge_ls["rounding"]=round(temp_edge_ls[,1],digit=0)
temp_edge_ls["dif"]=temp_edge_ls["rounding"]-temp_edge_ls[,1]
temp_edge_ls = temp_edge_ls[(temp_edge_ls$rounding!=0)
                            &(abs(temp_edge_ls$dif)<0.15),]

for(i in 1:nrow(temp_edge_ls)){
  temp_data = edge_ls_highcor[rownames(temp_edge_ls)[i],]
  temp_mz1 = df_raw$medMz[which(df_raw$ID==temp_data$node1)]
  temp_mz2 = temp_mz1 + temp_data$mz_dif
  for(j in 2:10){
    if(abs(temp_mz1*j-(H_mass-e_mass)*mode*(j-1)-temp_mz2)<temp_mz2*ppm){
      temp_df_oligo[(nrow(temp_df_oligo)+1),]=c(temp_data,"oligomer",paste(j,"-mer",sep=""))
    }
  }
}
  }
  rm(temp_edge_ls)
  temp_df_oligo
  nrow(temp_df_oligo)
  test_time = Sys.time()-test_time
  

  
#Multiple-charge, based on isotope features.
time["edge_annotation_multicharge"]=data.frame(Sys.time())
{
temp_ls_multi = list()
temp_df_multi = temp_df[0,]
for(n in 2:3){
  i=j=1
  while (i <= nrow(edge_ls_highcor)){
    i=i+1
    if(edge_ls_highcor$mz_dif[i]<(junk_df$mass[j]-search_ms_cutoff)/n){
      next
    }
    if(edge_ls_highcor$mz_dif[i]>(junk_df$mass[j]+search_ms_cutoff)/n){
      j=j+1
      if(j>nrow(junk_df)){break}
      while(edge_ls_highcor$mz_dif[i]>(junk_df$mass[j]-search_ms_cutoff)/n){i=i-1}
      next
    }
    temp_df_multi[(nrow(temp_df_multi)+1),]=c(edge_ls_highcor[i,],"multi_charge",paste(junk_df$type[j],"/",n,sep=""))
  }
  temp_ls_multi[[n]] = temp_df_multi
  temp_df_multi = temp_df_multi[0,]
}
}

temp_df_multi = bind_rows(temp_ls_multi)
temp_df_multi
nrow(temp_df_multi)

}

##Annotation summary
edge_ls_annotate=rbind(temp_df_isotope,temp_df_oligo,temp_df_multi)



##Put likelihood on annotation
time["edge_annotation_score"]=data.frame(Sys.time())

#mz_dif_score
{
edge_ls_annotate["mz_dif_score"]=0
junk_df["mass_group"]=0
n=1
junk_df$mass_group[1]=1
for(i in 2:nrow(junk_df)){
  if(junk_df$mass[i]-junk_df$mass[i-1]>search_ms_cutoff*2){n=n+1}
  junk_df$mass_group[i]=n
}

i=1
for (i in 1:max(junk_df$mass_group)){
  temp_junk_df = junk_df[junk_df$mass_group==i,]
  temp_data=edge_ls_annotate[edge_ls_annotate$type %in% temp_junk_df$type,]
  temp_data=unique(temp_data[,1:5])
  
  if(nrow(temp_junk_df)!=1){
    breaks <- 30
    his <- hist(temp_data$mz_dif*1000, breaks=breaks)  
    df <- data.frame(mid=his$mids, cou=his$counts)  
    guemea <- temp_junk_df$mass*1000
    guesig <- rep(initial_FIT$estimate[2], nrow(temp_junk_df))*1000
    guedis <- "norm"  
    emsteps = 20
    (fitpro <- mix(as.mixdata(df), mixparam(mu=guemea, sigma=guesig), 
                   constr=mixconstr(conmu ="MFX", fixmu=rep(T,nrow(temp_junk_df)),
                                    consigma = "SEQ"), 
                   dist=guedis,emsteps = emsteps))
    junk_df[junk_df$mass_group==i,"pi"]=fitpro$parameters$pi
    junk_df[junk_df$mass_group==i,"mu"]=fitpro$parameters$mu
    junk_df[junk_df$mass_group==i,"sigma"]=fitpro$parameters$sigma
  }
  else{
    temp_FIT <- fitdist(as.numeric(temp_data$mz_dif)*10^3, "norm")    
    #hist(temp_data$mz_dif)
    plot(temp_FIT)  
    # temp_data["mz_dif_score"]=dnorm(temp_data$mz_dif*1000, temp_FIT$estimate[1], temp_FIT$estimate[2])
    # temp_data["mz_dif_score"]=temp_data["mz_dif_score"]/max(temp_data["mz_dif_score"])
    junk_df[junk_df$mass_group==i,"pi"]=1
    temp_junk_df["pi"]=1
    junk_df[junk_df$mass_group==i,"mu"]=temp_FIT$estimate[1]
    junk_df[junk_df$mass_group==i,"sigma"]=temp_FIT$estimate[2]
  }
  # temp_data["mz_dif_score"]=dnorm(temp_data$mz_dif*1000, temp_FIT$estimate[1], temp_FIT$estimate[2])
  # temp_data["mz_dif_score"]=temp_data["mz_dif_score"]/max(temp_data["mz_dif_score"])
}
for (i in 1:nrow(junk_df)){
  temp_data=edge_ls_annotate[edge_ls_annotate$type==junk_df$type[i],]
  edge_ls_annotate[edge_ls_annotate$type==junk_df$type[i],"mz_dif_score"]=
    dnorm(temp_data$mz_dif*1000, junk_df$mu[i], junk_df$sigma[i])*junk_df$pi[i]
}

for (i in 1:max(junk_df$mass_group)){
  temp_junk_df = junk_df[junk_df$mass_group==i,]
  temp_data=edge_ls_annotate[edge_ls_annotate$type %in% temp_junk_df$type,]
  temp_data["mz_dif_score"]=temp_data["mz_dif_score"]/max(temp_data["mz_dif_score"])
  edge_ls_annotate[edge_ls_annotate$type %in% temp_junk_df$type,]=temp_data
}


}

#Statisically find the time_dif distribution for annotated peaks
{
time_distribute <- fitdist(as.numeric(edge_ls_annotate$time_dif)*10^3, "norm", method="mme")
plot(time_distribute)
hist(as.numeric(edge_ls_annotate$time_dif)*10^3)
#hist((temp_summary$time_dif)*10^3)
edge_ls_annotate["time_dif_score"]=dnorm(as.numeric(edge_ls_annotate$time_dif)*10^3, time_distribute$estimate[1],time_distribute$estimate[2])
edge_ls_annotate["time_dif_score"]=edge_ls_annotate["time_dif_score"]/max(edge_ls_annotate["time_dif_score"])
}

#Switch the direction of fragment category
{
  temp_data = edge_ls_annotate[edge_ls_annotate$category=="fragment",]
  
  temp_node=temp_data$node1
  temp_data$node1=temp_data$node2
  temp_data$node2=temp_node
  temp_data$time_dif=-temp_data$time_dif
  temp_data$mz_dif=-temp_data$mz_dif
  temp_data$log10_inten_ratio=-temp_data$log10_inten_ratio
  temp_data$type = paste("-",temp_data$type,sep="")
  
  edge_ls_annotate[edge_ls_annotate$category=="fragment",] = temp_data 
}



edge_ls_annotate$category=as.factor(edge_ls_annotate$category)
edge_ls_annotate["color"]=as.numeric(factor(edge_ls_annotate$category))
color_palette=c("blue", "blue", "black", "green", "black", "green")

edge_ls_annotate["Confidence"]=edge_ls_annotate$correlation+edge_ls_annotate$mz_dif_score+edge_ls_annotate$time_dif_score

edge_ls_confidence = edge_ls_annotate[edge_ls_annotate$Confidence>1.8,]

edge_ls_confidence = edge_ls_confidence[with(edge_ls_confidence,order(-Confidence)),]
edge_ls_confidence_unique = edge_ls_confidence[!duplicated(edge_ls_confidence[,c("node1","node2")]),]

write.csv(edge_ls_annotate,"edge_ls_annotate.csv", row.names = F)
write.csv(edge_ls_confidence_unique,"edge_ls_confidence_unique.csv", row.names = F)



##label df_raw and output 
df_raw$artifact=NA
for(i in 1:nrow(edge_ls_confidence_unique)){
  temp_node = edge_ls_confidence_unique$node2[i]
  temp_text = paste(edge_ls_confidence_unique$type[i],
                    "(",edge_ls_confidence_unique$Confidence[i],") of ",
                    edge_ls_confidence_unique$node1[i],
                    sep="")
  
  if(is.na(df_raw$artifact[df_raw$ID==temp_node])){df_raw$artifact[df_raw$ID==temp_node]=temp_text}
  else{df_raw$artifact[df_raw$ID==temp_node] = paste(df_raw$artifact[df_raw$ID==temp_node],temp_text, sep="; ")}
}

df_raw=df_raw[,colnames(df_raw)!="name"]
df_raw["annotation"]=NA
for(i in 1:nrow(df_raw)){
  if(!is.na(df_raw$Formula[i])){
    df_raw$annotation[i]=paste(df_raw$Formula[i],df_raw$Metabolite[i],df_raw$medRt[i],sep="_")
  } else if(!is.na(df_raw$artifact[i])){
    #df_raw$annotation[i]=paste(df_raw$artifact[i],df_raw$medRt[i],sep="_")
    df_raw$annotation[i]=paste("artifact",df_raw$medMz[i],df_raw$medRt[i],sep="_")
  } else {
    df_raw$annotation[i]=paste(df_raw$medMz[i],df_raw$medRt[i],sep="_")
  }
}

write.csv(df_raw,"hmdb_adduct.csv", row.names = F)







##Network analysis
time["network_analysis"]=Sys.time()

#Initiate network
{
g_annotate <- graph_from_data_frame(d = edge_ls_annotate, vertices = df_raw, directed = TRUE)
g_confidence = graph_from_data_frame(d = edge_ls_confidence, vertices = df_raw, directed = TRUE)
g_confidence_unique = graph_from_data_frame(d = edge_ls_confidence_unique, vertices = df_raw, directed = TRUE)

colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")

# g_highcor = graph_from_data_frame(d = edge_ls_highcor, vertices = df_raw, directed = FALSE)
# colors <- c("white", "red", "orange", "blue", "dodgerblue", "cyan")

## Understand "multi_charge" species belong to other category
# test_multi_charge["Connected"]=0
# for(i in 1:nrow(test_multi_charge)){
#   temp_connect = vertex.disjoint.paths(g_annotate, 
#                                      source = test_multi_charge$node1[i],
#                                      target = test_multi_charge$node2[i]
#                                      )
#   if(temp_connect){test_multi_charge$Connected[i]=1}
# }
# test_multi_charge=test_multi_charge[test_multi_charge$Connected==0,]
# 
node_ls=df_raw[,1:7]
node_ls[df_raw$ID==483,]
query_mzdif= 0.5
query_mz = 359.1
query=node_ls[df_raw$medMz<(query_mz+query_mzdif)&df_raw$medMz>(query_mz-query_mzdif),]
node_ls[grep("glucose", df_raw$Metabolite, ignore.case = T),]

}

#Whole network analysis
{
  
  #Generate network graph using node_list and edge_list above
  
  #str(edge_ls_highcor)
  #str(edge_ls_annotate)
  
  #Basic graph characteristics, distance, degree, betweeness
  basic_graph=function(g)
  {
    #distance
    farthest_vertices(g)
    #Degree
    g.degree <- degree(g, mode = c("all"))
    table(g.degree)
    # which.max(g.degree)
    node_ls_0degree = node_ls[g.degree==0,]
    # node_list$Exact_Mass[g.degree==0]
    # node_list$Exact_Mass[71]
    # node_list$MF[71]
    #Betweenness
    g.b <- betweenness(g, directed = T)
    V(g)$color <- colors[2]
    plot(g,
         vertex.label = NA,
         #edge.color = 'black',
         vertex.size = sqrt(g.b)+1,
         edge.arrow.size = 0.05,
         layout = layout_nicely(g))
  }
  basic_graph(g_annotate)
  #basic_graph(g_highcor)
  g.degree <- degree(g_annotate, mode = c("all"))
  dgree_table=table(g.degree)
  #which(g.degree==18)
}

#Analyze the network/subgraph of specific node
# edge_ls_Cr = edge_ls_annotate[grep("Cr", edge_ls_annotate$type),]
# g_Cr = graph_from_data_frame(d = edge_ls_Cr, vertices = df_raw, directed = TRUE)
# Cr_adduct = edge_ls_Cr[edge_ls_Cr$category=="adduct",]
# Cr_isotope = edge_ls_Cr[edge_ls_Cr$category=="natural_abundance",]
# filter = Cr_adduct[Cr_adduct$node2 %in% intersect(Cr_adduct$node2,Cr_isotope$node1),]
# g_cr = graph_from_data_frame(d = edge_ls_Cr, vertices = df_raw, directed = F)
# write.csv(filter,"Cr_list.csv", row.names = F)
# 
# for(i in intersect(Cr_adduct$node2,Cr_isotope$node1)){
# subgraph_specific_node(paste(i),g_cr)
# }

subgraph_specific_node = function(interested_node, g)
{
  #Glucose
  interested_node = interested_node
  g.degree <- degree(g, mode = c("all"))
  g_intrest <- make_ego_graph(g, 
                              max(g.degree), 
                              #1,
                              nodes = interested_node, 
                              mode = c("out"))[[1]]
  dists = distances(g_intrest, interested_node)
  colors <- c("black", "red", "orange", "blue", "dodgerblue", "cyan")
  #V(g_intrest)$color <- colors[dists+1]
  # png(filename=paste("Subnetwork of node ", interested_node,".png",sep=""),
  #     width = 2400, height=2400,
  #     res=300)
  # 
  plot(g_intrest,
       #vertex.color = 'white',
       vertex.label = vertex.attributes(g_intrest)$medMz,
       #vertex.label = vertex.attributes(g_intrest)$medRt,
       vertex.label.color = "red",
       vertex.label.cex = 1,
       #vertex.label.dist = 2,
       vertex.size = log(vertex.attributes(g_intrest)$mean_inten)*2-12,
       edge.width = edge.attributes(g_intrest)$Confidence*2-2,
       edge.color = color_palette[edge.attributes(g_intrest)$color],
       #edge.label = edge.attributes(g_intrest)$mz_dif,
       edge.label = edge.attributes(g_intrest)$type,
       edge.label.color = "red",
       edge.label.cex = 1,
       edge.arrow.size = 0.5,
       edge.arrow.width = 1,
       main = paste("Subnetwork of node", interested_node),
       layout = layout_nicely(g_intrest)
       #layout = layout_as_tree(g_intrest)
  )
  # dev.off()
}

#Glucose subgraph
subgraph_specific_node("10", g_confidence_unique)


time["end"]=Sys.time()
time_used=as.numeric()
for (i in 2:ncol(time)){
  time_used[i-1] = round(time[1,i]-time[1,i-1],2)
}
names(time_used)=colnames(time)[1:(ncol(time)-1)]
time_used


# #Plotting bar graph of different score distribution
# {
# hist(edge_ls_annotate$mz_dif_score)
# hist(edge_ls_annotate$correlation)
# hist(edge_ls_annotate$time_dif_score)
# hist(edge_ls_annotate$log10_inten_ratio)
# hist(edge_ls_annotate$Confidence)
# count=0
# for (junk_category in unique(junk_df$category)){
#   count=count+1
#   hist(edge_ls_annotate$mz_dif_score[edge_ls_annotate$category==junk_category],
#         xlab="mz_dif_score",
#         main=junk_category
#   )
#   hist(edge_ls_annotate$time_dif[edge_ls_annotate$category==junk_category],
#        xlab="time_dif_score",
#        main=junk_category
#   )
#   hist(edge_ls_annotate$correlation[edge_ls_annotate$category==junk_category],
#        xlab="intensity_correlation",
#        main=junk_category
#   )
#   hist(edge_ls_annotate$log10_inten_ratio[edge_ls_annotate$category==junk_category],
#        xlab="log10_inten_ratio",
#        main=junk_category
#   )
# }
# }
# #Venn diagram evaluation
# {
# mz_dif_score_0.2=edge_ls_annotate[edge_ls_annotate$mz_dif_score>0.2,]
# correlation_0.8=edge_ls_annotate[edge_ls_annotate$correlation>0.8,]
# time_dif_score_0.2=edge_ls_annotate[edge_ls_annotate$time_dif_score>0.2,]
# 
# n1=nrow(mz_dif_score_0.2)
# n2=nrow(correlation_0.8)
# n3=nrow(time_dif_score_0.2)
# n12=nrow(inner_join(mz_dif_score_0.2, correlation_0.8))
# n23=nrow(inner_join(correlation_0.8, time_dif_score_0.2))
# n13=nrow(inner_join(mz_dif_score_0.2, time_dif_score_0.2))
# n123=nrow(inner_join(mz_dif_score_0.2, inner_join(correlation_0.8,time_dif_score_0.2)))
# grid.newpage()
# venn.plot <- draw.triple.venn(n1, n2, n3,
#                               n12, n23, n13, n123, c("mz_dif_score=0.2", "correlation=0.8", "time_dif_score=0.2"));
# grid.draw(venn.plot)
# }


##Debug
{
  # edge_ls_no_annotate=edge_ls_highcor[edge_ls_highcor$annotate==0,]
  # 
  # #Is two node connected?
  # edge_ls_no_annotate["Connected"]=0
  # for(i in 1:nrow(edge_ls_no_annotate)){
  #   temp_connect = vertex.disjoint.paths(g_annotate, 
  #                                        source = edge_ls_no_annotate$node1[i],
  #                                        target = edge_ls_no_annotate$node2[i])
  #   if(temp_connect){edge_ls_no_annotate$Connected[i]=1}
  # }
  # edge_ls_no_connect=edge_ls_no_annotate[edge_ls_no_annotate$Connected==0,]
  # 
  # 
  # 
  # 
  # test = fitdist(as.numeric(temp_summary$time_dif)*10^3,
  #                "pois",
  #                fix.arg=list(a=0),
  #                start = list(mean = mean(as.numeric(temp_summary$time_dif)*10^3),
  #                             sd = sd(as.numeric(temp_summary$time_dif)*10^3)),
  #                optim.method="L-BFGS-B",
  #                lower=c(0, 0)
  # )
}