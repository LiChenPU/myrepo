#显示中文 
# !diagnostics off
#Sys.setlocale(category = "LC_ALL", locale = "Chinese")
library(readr)
library(tidyr)
library(ggplot2)
#install.packages("stringi")
library(stringi)
#install.packages("matrixStats")
library(matrixStats)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

filenames = c("Yeast-Ecoli-neg-peakpicking_blank.csv")
hmdb_df = read_csv("hmdb_unique.csv")
known_mz = read_csv("known_mz.csv")

mode = -1
normalized_to_col_median = F

time=Sys.time()


#Consider two groups have the same m/z if difference is below threshold
ms_dif_ppm=5
ms_dif_ppm = ms_dif_ppm/10^6
#Consider two groups have the same RT if difference is below threshold
rt_dif_min=0.2
#Signal below detection_limit will be replaced as detection_limit
detection_limit=2500



#Read data
{
  filename=filenames[1]
  raw <- read_csv(filename)
  #Remove rows that has incomplete intensity data
  raw = raw[complete.cases(raw[, 14:ncol(raw)]),]
  
  raw_rn = nrow(raw)
  sample_names=colnames(raw)[15:ncol(raw)]
  if(length(grep("blank|blk", sample_names, ignore.case = T))!=0){
    sample_names_noblank=sample_names[-grep("blank|blk", sample_names, ignore.case = T)]
  } else {
    sample_names_noblank=sample_names
  }
  sample_names_blank=sample_names[grep("blank|blk", sample_names, ignore.case = T)]
  sample_cohort=stri_replace_last_regex(sample_names,'-1|-2|-3|-a|-b|-c|_mean', '',stri_opts_regex(case_insensitive=T))
  
}

#Data manipulation to increase replicate number
{
  data_replicate_number = 3
  if(data_replicate_number <3)
  {
    for (i in names(table(sample_cohort))){
      raw[paste(i,"_mean",sep="")]=rowMedians(raw[, 14+which(sample_cohort==i)])
    }
  }
  sample_names=colnames(raw)[15:ncol(raw)]
  if(length(grep("blank|blk", sample_names, ignore.case = T))!=0){
    sample_names_noblank=sample_names[-grep("blank|blk", sample_names, ignore.case = T)]
  } else {
    sample_names_noblank=sample_names
  }
  sample_names_blank=sample_names[grep("blank|blk", sample_names, ignore.case = T)]
  sample_cohort=stri_replace_last_regex(sample_names,'-1|-2|-3|-a|-b|-c|_mean', '',stri_opts_regex(case_insensitive=T))
  
}

###### Don't need to change from here #####
### Merge similar m/z and RT, prepare files for MetaboAnalyst ####
##Group MS groups
{
  s = raw[with(raw, order(medMz, medRt)),]
  s["merge_group"]=NA
  mgMS_count = 1 
  i_max=i_min=1
  
  #
  while (i_min <= nrow(raw)){
    while(s$medMz[i_max]-s$medMz[i_min]< (s$medMz[i_min]*ms_dif_ppm)){
      i_max = i_max+1
      if(i_max>nrow(raw)){
        break
      }
    }
    if(i_max-i_min == 1){
      s$merge_group[i_min] = 0
      i_min = i_min + 1
      next
    }
    while(i_min < i_max){
      s$merge_group[i_min] = mgMS_count
      i_min=i_min+1
    }  
    mgMS_count=mgMS_count+1
  }
}

##Group RT similar groups based on MS groups
{
  s2 = s[with(s, order(merge_group, medRt)),]
  s2["RTmerge_group"]=NA
  
  mgRT_count=1
  
  j_min=length(which(s2$merge_group == 0))
  
  for (n in 1:j_min){
    s2$RTmerge_group[n] = 0
  }
  
  j_min=j_min+1
  if(j_min <=nrow(raw)){s2$RTmerge_group[j_min] = 0}
  j_max=j_min+1
  
  while (j_max <= nrow(raw)){
    while (s2$merge_group[j_min] != s2$merge_group[j_max] |
           s2$medRt[j_max] - s2$medRt[j_min] >= rt_dif_min){
      s2$RTmerge_group[j_max] = 0
      j_max = j_max+1
      j_min = j_min+1
      if(j_max > nrow(raw)){break}
    }
    if(j_max > nrow(raw)){break}
    
    while (s2$medRt[j_max] - s2$medRt[j_max-1] < rt_dif_min){
      j_max = j_max+1
      if(j_max>nrow(raw) | s2$merge_group[j_min] != s2$merge_group[j_max]) {break}
    }
    j_max = j_max-1
    
    for(j in j_min:j_max){
      s2$RTmerge_group[j] = mgRT_count
    }
    
    mgRT_count = mgRT_count+1
    j_min=j_max
    j_max=j_max+1
  }
}

##Flag groups for deletion & combine signal
{
  s3 = s2[with(s2, order(RTmerge_group)),]
  s3[["flag"]]=NA
  k_max=k_min=length(which(s3$RTmerge_group==0))
  for(k in 1:k_min){
    s3$flag[k] = T
  }
  
  while (k_max <= nrow(raw)){
    k_min = k_max
    s3$flag[k_min] = T
    while (s3$RTmerge_group[k_min] == s3$RTmerge_group[k_max]){
      k_max = k_max+1
      if(k_max > nrow(raw)){break}
      s3$flag[k_max] = F
    }
    
    if(k_min == k_max-1){
      next
    }
    
    s3$medMz[k_min]=mean(s3$medMz[k_min:(k_max-1)])
    s3$medRt[k_min]=mean(s3$medRt[k_min:(k_max-1)])
    s3$goodPeakCount[k_min]=max(s3$goodPeakCount[k_min:(k_max-1)])
    for (n in 14:ncol(raw)){
      s3[k_min,n]=max(s3[k_min:(k_max-1),n])
    }
  }
}

#intermediate files, replace below detection number to random small number
{
  s4 = s3[s3$flag, 1:ncol(raw)]
  out_filename = paste("merge_", filename, sep="")
  # hist(result$goodPeakCount)
  # hist(result$maxQuality)
  s4 = s3[s3$flag, 1:ncol(raw)]
  s4 = s4[,-c(1:2,4, 7:14)]
  colnames(s4)[1]="ID"
  s4[,4:ncol(s4)][s4[,4:ncol(s4)]<detection_limit]=sample(1:detection_limit, 
                                                          size=sum(s4[,4:ncol(s4)]<detection_limit), 
                                                          replace=T)
}
  

##Perform row-mean and col-median normalization
{
H_mass = 1.00782503224
e_mass = 0.00054857990943
known_mz["adjustmz"] = 0
known_mz["adjustmz"] = known_mz$mz + (H_mass-e_mass)*mode
known_mz=known_mz[with(known_mz, order(adjustmz)),]

s5=s4
s5 = s5[with(s5, order(medMz)),]

s5["Formula"]=as.character()

s5["Metabolite"]=as.character()
i_min=1
j=1
ppm=5/10^6
while(j<=nrow(s5)){
  while(known_mz$adjustmz[i_min+1] < (s5$medMz[j]*(1-ppm))){
    i_min=i_min+1
  }
  temp_metabolite=as.character()
  temp_formula=as.character()
  k=1
  while(known_mz$adjustmz[i_min+k]<(s5$medMz[j]*(1+ppm))){
    if(length(temp_metabolite)==0){
      temp_metabolite=known_mz$Name[i_min+k]
      temp_formula=known_mz$MF[i_min+k]
    }
    else{
      temp_metabolite=paste(temp_metabolite,";",known_mz$Name[i_min+k])
      temp_formula=paste(temp_formula,";",known_mz$MF[i_min+k])
    }
    k=k+1
  }
  if(length(temp_metabolite)!=0){
    s5$Metabolite[j]=temp_metabolite
    s5$Formula[j]=temp_formula
  }
  j=j+1
}
s6=s5[complete.cases(s5),]
s6[,4:(ncol(s6)-2)]=s6[,4:(ncol(s6)-2)]/rowMedians(as.matrix(s6[,4:(ncol(s6)-2)]))
s6[,4:(ncol(s6)-2)]=s6[,4:(ncol(s6)-2)]/rowMeans(as.matrix(s6[,4:(ncol(s6)-2)]))
col_median = data.frame(col_median=apply(s6[,4:(ncol(s6)-2)],2,median))
 
col_median
write.csv(col_median, paste("col_median_",filename,sep=""), row.names = T)
}


## Annontate base on HMDB mz
{
hmdb_df["adjustmz"] = hmdb_df$Exact_Mass + (H_mass-e_mass)*mode
hmdb_df=hmdb_df[with(hmdb_df, order(adjustmz)),]

s5["Formula"]=as.character()
s5["Metabolite"]=as.character()

i_min=1
j=1
ppm=5/10^6
while(j<=nrow(s5)){
  while(hmdb_df$adjustmz[i_min+1] < (s5$medMz[j]*(1-ppm))){
    i_min=i_min+1
  }
  temp_metabolite=as.character()
  temp_formula=as.character()
  k=1
  while(hmdb_df$adjustmz[i_min+k]<(s5$medMz[j]*(1+ppm))){
    if(length(temp_metabolite)==0){
      temp_metabolite=hmdb_df$Name[i_min+k]
      temp_formula=hmdb_df$MF[i_min+k]
    }
    else{
      temp_metabolite=paste(temp_metabolite,";",hmdb_df$Name[i_min+k])
      temp_formula=paste(temp_formula,";",hmdb_df$MF[i_min+k])
    }
    k=k+1
  }
  if(length(temp_metabolite)!=0){
    s5$Metabolite[j]=temp_metabolite
    s5$Formula[j]=temp_formula
  }
  j=j+1
}
}



##Remove high blank samples
s7 = s5
s7["high_blank"]=NA
if(length(sample_names_blank)>0){
  s7["high_blank"]= rowMeans(s7[,sample_names_noblank]) < 2*rowMeans(s7[,sample_names_blank])
  #test = s7[s7$high_blank&is.na(s7$Formula),]
}




##Output files for MetaboAnalyst
{
  merge_output = s4[, !colnames(s4) %in% c("medRt","medMz")]
  #merge_output = merge_output[,-grep("blank|blk", colnames(merge_output), ignore.case = T)]
  if(normalized_to_col_median){
    merge_output=cbind(merge_output[,1],sweep(merge_output[,-1],2,col_median$col_median,FUN = "/"))
  }
  
  MA_output = rbind(c("Cohort", sample_cohort), merge_output)
  if(length(sample_names_blank)>0){
    MA_output = MA_output[,-grep("blank|blk", colnames(MA_output), ignore.case = T)]
  }
  
  # out_filename = paste("merge_", filename, sep="")
  # write.csv(merge_output, file=out_filename,row.names = F)
  MetaboAnalyst_filename = paste("MA_", filename, sep="")
  write.csv(MA_output, file=MetaboAnalyst_filename, row.names=F)
}

{
  library(MetaboAnalystR)
  mSet <- InitDataObjects("pktable", "stat", FALSE)
  mSet<-Read.TextData(mSet, MetaboAnalyst_filename, "colu", "disc");
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet);
  mSet<-FilterVariable(mSet, "iqr", "F", 25)
  mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", "a11", ratio=FALSE, ratioNum=20)
  mSet<-ANOVA.Anal(mSet, F, 0.05, "fisher")
  #mSet<-PlotHeatMap(mSet, "heatmap_0_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview", T, T, NA, T, F)
  #mSet<-PlotSubHeatMap(mSet, "heatmap_1_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 25, "overview", F, T, T, F)
  #mSet<-PlotSubHeatMap(mSet, "heatmap_1_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 50, "overview", T, T, T, F)

  ANOVA_file = "anova_posthoc.csv"
  ANOVA_raw <- read_csv(ANOVA_file)
  
  ANOVA_FDR = ANOVA_raw[,c("X1","FDR")]
  ANOVA_FDR$FDR=-log10(ANOVA_FDR$FDR)
  colnames(ANOVA_FDR)=c("ID","_log10_FDR")
}

##Output
{
  hmdb_match_output=s7[with(s7, order(ID)),]
  merge_output=merge_output[with(merge_output,order(ID)),]
  
  hmdb_match_output[sample_names_noblank]=merge_output[sample_names_noblank]
  
  hmdb_match_output=merge(hmdb_match_output,ANOVA_FDR,by="ID",all=T)
  hmdb_match_output$`_log10_FDR`[is.na(hmdb_match_output$`_log10_FDR`)]=0
  
  refcols = c("ID","medMz","medRt","Formula","Metabolite","_log10_FDR","high_blank")
  hmdb_match_output = hmdb_match_output[,c(refcols, setdiff(colnames(hmdb_match_output), refcols))]
}

if(normalized_to_col_median){
  write.csv(hmdb_match_output, file=paste("hmdb_","normalized_col_median_",filename,sep=""), row.names=F)
}else{
  write.csv(hmdb_match_output, file=paste("hmdb_","without_normalization_",filename,sep=""), row.names=F)
}


fn <-MetaboAnalyst_filename
if (file.exists(fn)) file.remove(fn)

fn <-ANOVA_file
if (file.exists(fn)) file.remove(fn)


# 
# hmdb_match_output = read_csv("hmdb_Yeast-Ecoli-neg-peakpicking_blank.csv")
# test = hmdb_match_output
# 
# test_hasformula = hmdb_match_output[!is.na(hmdb_match_output$Formula),]
# test_lowblank = hmdb_match_output[!hmdb_match_output$high_blank,]
# test_significantpeaks = hmdb_match_output[hmdb_match_output$`_log10_FDR`!=0,]
# 
# test_output=unique(rbind(test_hasformula,test_lowblank,test_significantpeaks))





# 
# ##Include alternative names and SMILE 
# 
# hmdb_full_df = read_csv("hmdb_full.csv")
# 
# hmdb_SMILE = hmdb_full_df[,c("MF","Name","HMDB_ID","SMILES") ]
# 
# s6 = s5[!is.na(s5$Formula),]
# 
# ls = list()
# for (i in 1:nrow(s6)){
#   temp_df = hmdb_SMILE[hmdb_SMILE$MF==s6$Formula[i],]
#   temp_data = s6[rep(i,nrow(temp_df)),1:(ncol(s6)-2)]
#   ls[[i]] = cbind(temp_df, temp_data)
# }
# 
# s7 = bind_rows(ls)
# 
# 
# write.csv(s7, file=paste("ful_hmdb_SMILE_",filename,sep=""), row.names=F)



time = Sys.time() - time
